#include "gcre.h"

#include <chrono>
#include <thread>
#include <atomic>
#include <mutex>

using namespace std;

JoinExec::JoinExec(const int num_cases, const int num_ctrls) : num_cases(num_cases), num_ctrls(num_ctrls), width_ul(vector_width_ul(num_cases, num_ctrls)) {

  // create case mask from case/control ranges
  case_mask = new uint64_t[width_ul];
  for(int k = 0; k < width_ul; k++)
    case_mask[k] = bit_zero_ul;
  for(int k = 0; k < num_cases; k++)
    case_mask[k/64] |= bit_one_ul << k % 64;

  printf("init exec: %d x %d (width: %d ul)\n", num_cases, num_ctrls, width_ul);
}

// TODO size check and precompute max table
// create a copy of the vt for simplicity
void JoinExec::setValueTable(vec2d_d table){
  printf("setting value table: %lu x %lu\n", table.size(), table.size() > 0 ? table.front().size() : 0);
  value_table = table;
}

// data comes in as perm x patient: 0=flipped, 1=unflipped
// here converted to 1=case, 0=control
// intermediate mask is inverted, to keep 0-padding of the vector tail from becoming mostly cases
void JoinExec::setPermutedCases(const vec2d_i& data) {

  printf("loading permuted case masks: %lu x %lu (iterations: %d)\n", data.size(), data.size() > 0 ? data.front().size() : 0, iterations);

  if(iterations < data.size())
    printf("  ** WARN more permuted cases than iterations, input will be truncated\n");

  perm_case_mask = new uint64_t[iterations * width_ul];
  for(int k = 0; k < iterations * width_ul; k++)
    perm_case_mask[k] = bit_zero_ul;

  for(int r = 0; r < min( (size_t) iterations, data.size()); r++) {
    uint64_t perm[width_ul];
    for(int k = 0; k < width_ul; k++)
      perm[k] = bit_zero_ul;
    for(int c = 0; c < data[r].size(); c++) {
      if(data[r][c] != 1)
        perm[c/64] |= bit_one_ul << c % 64;
    }
    for(int k = 0; k < width_ul; k++) {
      perm_case_mask[r * width_ul + k] = case_mask[k] ^ perm[k];
    }
  }

  if(iterations > data.size()) {

    printf("  ** WARN not enough permuted cases, some will be reused to match iterations - hope this is for testing!\n");

    for(int r = data.size(); r < iterations; r++) {
      int c = r % data.size();
      for(int k = 0; k < width_ul; k++) 
        perm_case_mask[r * width_ul + k] = perm_case_mask[c * width_ul + k];
    }
  }

}

// TODO width based on method needs
unique_ptr<PathSet> JoinExec::createPathSet(int size) const {
  // std::make_unique is C++14
  return unique_ptr<PathSet>(new PathSet(size, width_ul, width_ul * 2));
}

joined_res JoinExec::join(const vector<uid_ref>& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const {

  // TODO
  // int iters = 8 * ceil(max(1, conf.iterations) / 8.0);
  const int iters = iterations;
  const int vlen = width_ul;
  // printf("adjusting width: %d -> %d\n", conf.num_cases + conf.num_controls, vlen * 64);
  // printf("adjusting iterations: %d -> %d\n\n", conf.iterations, iters);

  // priority queue for the indices of the top K paths in the data
  // add dummy score to avoid empty check
  priority_queue<Score> global_scores;
  global_scores.push(Score());

  atomic<unsigned> uid_idx(0);
  mutex g_mutex;

  // C++14 capture init syntax would be nice,
  // const auto exec = this;
  auto worker = [this, &uids, &uid_idx, &g_mutex, iters, &paths0, &paths1, &paths_res, &global_scores](int tid) {

    bool keep_paths = paths_res.size != 0;
    int flipped_pivot_length = paths1.size;

    priority_queue<Score> scores;
    scores.push(Score());

    double perm_score[iters];
    double perm_flips[iters];
    for(int k = 0; k < iters; k++){
      perm_score[k] = 0;
      perm_flips[k] = 0;
    }

    // allocate as single block for storage in path set
    uint64_t joined[width_ul * 2];
    uint64_t* joined_pos = joined;
    uint64_t* joined_neg = joined + width_ul;

    uint64_t joined_tp[width_ul];
    uint64_t joined_tn[width_ul];

    int p_case_pos[iters];
    int p_ctrl_pos[iters];

    int idx = -1;
    while((idx = uid_idx.fetch_add(1)) < uids.size()) {

      const auto& uid = uids[idx];

      if(uid.count == 0)
        continue;

      // start of result path block for this uid
      int path_idx = uid.path_idx;

      // printf("[thread-%d] idx: %d, src: %d, trg: %d, count: %d, loc: %d --", tid, idx, uid.src, uid.trg, uid.count, uid.location);

      const uint64_t* path_pos0 = paths0[idx];
      const uint64_t* path_neg0 = path_pos0 + width_ul;

      for(int loc_idx = 0, loc = uid.location; loc < uid.location + uid.count; loc_idx++, loc++) {

        auto sign = uid.signs[loc_idx];
        const uint64_t* path_pos1 = (sign ? paths1[loc] : paths1[loc] + width_ul);
        const uint64_t* path_neg1 = (sign ? paths1[loc] + width_ul : paths1[loc]);

        int case_pos = 0;
        int case_neg = 0;
        int ctrl_pos = 0;
        int ctrl_neg = 0;
        int total_pos = 0;
        int total_neg = 0;

        uint64_t bit_pos, bit_neg, bit_con, true_pos, true_neg, mask;

        // zero out join path holder
        memset(joined_pos, 0, width_ul * sizeof(uint64_t));
        memset(joined_neg, 0, width_ul * sizeof(uint64_t));

        // for(int r = 0; r < iters; r++){
        //   p_case_pos[r] = 0;
        //   p_ctrl_pos[r] = 0;
        // }

        for(int k = 0; k < width_ul; k++) {

          bit_pos = path_pos0[k] | path_pos1[k];
          bit_neg = path_neg0[k] | path_neg1[k];

          if(bit_pos != 0 || bit_neg != 0) {

            bit_con = bit_pos & bit_neg;
            true_pos = bit_pos & ~bit_con;
            true_neg = bit_neg & ~bit_con;
            mask = case_mask[k];

            total_pos += __builtin_popcountl(true_pos);
            total_neg += __builtin_popcountl(true_neg);
            case_pos += __builtin_popcountl(true_pos &  mask);
            case_neg += __builtin_popcountl(true_neg & ~mask);
            ctrl_pos += __builtin_popcountl(true_neg &  mask);
            ctrl_neg += __builtin_popcountl(true_pos & ~mask);

            // TODO permute mask access
/*
            if(true_pos != 0) {
              for(int r = 0; r < iters; r++) {
                // uint64_t* p_mask = permute_mask + r * vlen;
                // p_case_pos[r] += __builtin_popcountl(true_pos & p_mask[k]);
              }
            }

            if(true_neg != 0) {
              for(int r = 0; r < iters; r++){
                // uint64_t* p_mask = permute_mask + r * vlen;
                // p_ctrl_pos[r] += __builtin_popcountl(true_neg & p_mask[k]);
              }
            }
*/

            joined_pos[k] = bit_pos;
            joined_neg[k] = bit_neg;

          }

        }

        if(keep_paths)
          paths_res.set(path_idx, joined);
        path_idx += 1;

        int cases = case_pos + case_neg;
        int ctrls = ctrl_pos + ctrl_neg;

        double score = value_table[cases][ctrls];
        double flips = value_table[ctrls][cases];

        if(score > scores.top().score)
          scores.push(Score(score, idx, loc));
        if(flips > scores.top().score)
          scores.push(Score(flips, idx, loc + flipped_pivot_length));
        while(scores.size() > top_k)
          scores.pop();

        // for(int r = 0; r < conf.iterations; r++){
        //   int p_cases = p_case_pos[r] + (total_neg - p_ctrl_pos[r]);
        //   int p_ctrls = p_ctrl_pos[r] + (total_pos - p_case_pos[r]);

        //   double p_score = value_table[p_cases][p_ctrls];
        //   if(p_score > perm_score[r])
        //     perm_score[r] = p_score;
        //   double p_flips = value_table[p_ctrls][p_cases];
        //   if(p_flips > perm_flips[r])
        //     perm_flips[r] = p_flips;
        // }

      }

    }

    // synchronize final section to merge results
    lock_guard<mutex> lock(g_mutex);

    while(!scores.empty()){
      auto score = scores.top();
      if(score.score > global_scores.top().score)
        global_scores.push(score);
      scores.pop();
    }

    printf("[%d] done\n", tid);
  };

  if(nthreads > 0) {
    vector<thread> pool(nthreads);
    for(int tid = 0; tid < pool.size(); tid++)
      pool[tid] = thread(worker, tid + 1);
    for(auto& t : pool)
      t.join();
  } else {
    // '0' is main thread
    worker(0);
  }

  while(global_scores.size() > top_k)
    global_scores.pop();

  joined_res res;
  res.permuted_scores.resize(iterations, 0);
  for(int k = 0; k < iterations; k++)
    res.permuted_scores[k] = 0; // max(perm_score[k], perm_flips[k]);
  res.scores.clear();
  while(!global_scores.empty()){
    res.scores.push_back(global_scores.top());
    global_scores.pop();
  }

  return res;
}

// TODO method to generate case permutations
  // printf("generating simplified case masks for permutation (it = %d)\n", conf.iterations);
  // srand((unsigned) time(0));

  // vec2d_u16 permute_cases;
  // for(int k = 0; k < conf.iterations; k++){
  //   set<int> cases;
  //   permute_cases.push_back(vector<uint16_t>());
  //   while(cases.size() < conf.num_cases)
  //     cases.insert(rand() % (conf.num_cases + conf.num_controls));
  //   for(auto c : cases)
  //     permute_cases.back().push_back(c);
  // }
  // printf("permuted_cases: %lu (%lu)\n", permute_cases.size(), permute_cases.back().size());
