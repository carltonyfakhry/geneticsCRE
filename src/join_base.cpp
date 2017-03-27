#include "gcre.h"

#include <thread>
#include <atomic>
#include <chrono>

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
// intermediate mask is inverted, so keep 0-padding of the vector tail from becoming mostly cases
void JoinExec::setPermutedCases(vec2d_i& data) {

  printf("loading permuted case masks: %lu x %lu (iterations: %d)\n", data.size(), data.size() > 0 ? data.front().size() : 0, iterations);

  perm_case_mask = new uint64_t[iterations * width_ul];
  for(int k = 0; k < iterations * width_ul; k++)
    perm_case_mask[k] = bit_zero_ul;

  for(int r = 0; r < data.size(); r++) {
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

}

unique_ptr<PathSet> JoinExec::createPathSet(int size) const {
  // std::make_unique is C++14
  return unique_ptr<PathSet>(new PathSet_BlockM2(size, num_cases + num_ctrls));
}

joined_res JoinExec::join(const vector<uid_ref>& uids, const vector<int>& join_gene_signs, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const {


  // TODO
  // int iters = 8 * ceil(max(1, conf.iterations) / 8.0);
  const int iters = iterations;
  const int vlen = width_ul;
  // printf("adjusting width: %d -> %d\n", conf.num_cases + conf.num_controls, vlen * 64);
  // printf("adjusting iterations: %d -> %d\n\n", conf.iterations, iters);


  bool keep_paths = paths_res.size != 0;
  int flipped_pivot_length = paths1.size;
  int path_idx = 0;

  // priority queue for the indices of the top K paths in the data
  // add dummy score to avoid empty check
  priority_queue<Score> scores;
  scores.push(Score());

  atomic<unsigned> uid_idx(0);

  // C++14 capture init syntax would be nice
  // const auto exec = this;
  auto worker = [this, &uids, &uid_idx](int tid) {

    // uint64_t joined_pos[vlen];
    // uint64_t joined_neg[vlen];
    // uint64_t joined_tps[vlen];
    // uint64_t joined_tns[vlen];

    int idx = -1;
    while((idx = uid_idx.fetch_add(1)) < uids.size()) {
      const auto& uid = uids[idx];
      // printf("[thread-%d] idx: %d, src: %d, trg: %d, count: %d, loc: %d --", tid, idx, uid.src, uid.trg, uid.count, uid.location);
      if(uid.count == 0)
        continue;

      // for(int loc = uid.location; loc < uid.location + uid.count; loc++)
        // printf(" %d", loc);
      // printf("\n");
      // this_thread::sleep_for(chrono::milliseconds(10));
    }
    printf("[%d] done.\n", tid);
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

  // exit(0);
  return joined_res();


/*
  for(int i = 0; i < uids.size(); i++){

    // ** thread partition
    // int i
    // uid_ref& uid
    // -- sign stuff
    // uint64_t* left (1 x width_ul)
    // uint64_t* right (uid.count x width_ul)
    // uint64_t* to store (if keep paths)


    const uid_ref& uid = uids[i];
    printf("[%06d] src: %d, trg: %d, count: %d, loc: %d\n", i, uid.src, uid.trg, uid.count, uid.location);
    if(uid.count == 0)
      continue;

    // uint64_t* path_pos0 = paths0->pos + i * vlen;
    // uint64_t* path_neg0 = paths0->neg + i * vlen;

    for(int j = uid.location; j < (uid.location + uid.count); j++){

      // printf("[%06d] src: %d, trg: %d, count: %d, loc: %d --> %d\n", i, uid.src, uid.trg, uid.count, uid.location, j);


      int sign = 0;
      // if(conf.path_length > 3)
      //   sign = join_gene_signs[i];
      // else if(conf.path_length < 3)
      //   sign = join_gene_signs[j];
      // else if(conf.path_length == 3)
      //   sign = (join_gene_signs[i] + join_gene_signs[j] == 0) ? -1 : 1;

      // uint64_t* path_pos1 = (sign == 1 ? paths1->pos : paths1->neg) + j * vlen;
      // uint64_t* path_neg1 = (sign == 1 ? paths1->neg : paths1->pos) + j * vlen;


      if(keep_paths){
        // memcpy(pathsr->pos + path_idx * vlen, joined_pos, vlen * sizeof(uint64_t));
        // memcpy(pathsr->neg + path_idx * vlen, joined_neg, vlen * sizeof(uint64_t));
      }
      path_idx += 1;

      // int cases = case_pos + case_neg;
      // int ctrls = ctrl_pos + ctrl_neg;


    }

  }
*/
  return joined_res();
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
