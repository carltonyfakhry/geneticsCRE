#include <thread>
#include <atomic>
#include <mutex>

#include "gcre.h"

using namespace std;

// TODO uid and loc index integer sizes

class JoinMethod {

public:

  priority_queue<Score> scores;
  double* perm_scores;

  JoinMethod(const JoinExec* exec, const int flip_pivot_len, double* p_perm_scores, uint64_t* p_joined_block, int* p_perm_count_block) :

  exec(exec),
  flip_pivot_len(flip_pivot_len),
  case_mask(exec->case_mask),
  perm_case_mask(exec->perm_case_mask),
  value_table(exec->value_table),
  value_table_max(exec->value_table_max) {

    scores.push(Score());

    perm_scores = p_perm_scores;
    for(int k = 0; k < exec->iterations; k++)
      perm_scores[k] = 0;

    joined_block = p_joined_block;
    perm_count_block = p_perm_count_block;

  }

  // run at end of worker thread to collect results
  void drain_scores(const JoinExec* exec) {

    while(!scores.empty()) {
      auto score = scores.top();
      if(score.score > exec->scores.top().score)
        exec->scores.push(score);
      scores.pop();
    }

    for(int r = 0; r < exec->iterations; r++) {
      if(perm_scores[r] > exec->perm_scores[r])
        exec->perm_scores[r] = perm_scores[r];
    }

  }

protected:

  const JoinExec* exec;
  const int flip_pivot_len;

  const uint64_t* case_mask;
  const uint64_t* perm_case_mask;

  const vec2d_d& value_table;
  const vec2d_d& value_table_max;

  uint64_t* joined_block;
  int* perm_count_block;

};

class JoinMethod1 : public JoinMethod {

public:

  JoinMethod1(const JoinExec* exec, const int flip_pivot_len, double* p_perm_scores, uint64_t* p_joined_block, int* p_perm_count_block) : JoinMethod(exec, flip_pivot_len, p_perm_scores, p_joined_block, p_perm_count_block) {
  }

  void score_permute_cpu(int idx, int loc, const uint64_t* path0, const uint64_t* path1);

};

class JoinMethod2 : public JoinMethod {

public:
  JoinMethod2(const JoinExec* exec, const int flip_pivot_len, double* p_perm_scores, uint64_t* p_joined_block, int* p_perm_count_block) : JoinMethod(exec, flip_pivot_len, p_perm_scores, p_joined_block, p_perm_count_block) {

    joined_pos = joined_block;
    joined_neg = joined_block + exec->width_ul;

    perm_case_pos = perm_count_block;
    perm_ctrl_pos = perm_count_block + exec->iterations;

  }

  void score_permute_cpu(int idx, int loc, const uint64_t* path0_pos, const uint64_t* path0_neg, const uint64_t* path1_pos, const uint64_t* path1_neg);

protected:

  uint64_t* joined_pos;
  uint64_t* joined_neg;
  int* perm_case_pos;
  int* perm_ctrl_pos;

};

JoinExec::JoinExec(string method_name, int num_cases, int num_ctrls, int iters) :

method(to_method(method_name)),
num_cases(num_cases),
num_ctrls(num_ctrls),
width_ul(vector_width_ul(num_cases, num_ctrls)),
iterations(iter_size(iters)) {

  check_true(num_cases > 0 && num_ctrls > 0 && iters > 0);

  width_vec = width_ul / gs_vec_width_ul;

  // create case mask from case/control ranges
  case_mask = new uint64_t[width_ul];
  for(int k = 0; k < width_ul; k++)
    case_mask[k] = bit_zero_ul;
  for(int k = 0; k < num_cases; k++)
    case_mask[k/64] |= bit_one_ul << k % 64;

  printf("[init exec] method: '%d' | data: %d x %d (width: %d ul)\n", method, num_cases, num_ctrls, width_ul);
}

// create a copy of the vt for simplicity
void JoinExec::setValueTable(vec2d_d table){
  // TODO size check
  printf("setting value table: %lu x %lu\n", table.size(), table.size() > 0 ? table.front().size() : 0);

  // indices can be flipped and counts could include all 0 to all cases/controls (so max index == size)
  int max_idx = max(num_cases, num_ctrls);
  check_index(max_idx, table.size());
  for(const auto& row : table)
    check_index(max_idx, row.size());

  value_table = table;
  // precompute max including flipped indices
  value_table_max = vec2d_d(value_table.size(), vec_d(value_table.front().size()));
  for(int r = 0; r < value_table.size(); r++) {
    for(int c = 0; c < value_table[r].size(); c++) {
      value_table_max[r][c] = max(value_table[r][c],value_table[c][r]);
    }
  }

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
    check_equal(num_cases + num_ctrls, data[r].size());
    for(int c = 0; c < data[r].size(); c++) {
      if(data[r][c] != 1)
        perm[c/64] |= bit_one_ul << c % 64;
    }
    int count_cases = 0;
    for(int k = 0; k < width_ul; k++) {
      uint64_t mask_k = case_mask[k] ^ perm[k];
      perm_case_mask[r * width_ul + k] = mask_k;
      count_cases += __builtin_popcountl(mask_k);
    }
    if(count_cases != num_cases)
      printf("  ** WARN perm mask did not set correct number of cases: %d != %d (mask %04d)\n", count_cases, num_cases, r);
  }

  // TODO correct padding behavior

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
unique_ptr<PathSet> JoinExec::createPathSet(st_pathset_size size) const {
  // std::make_unique is C++14
  // also the vlen size parameter is a ridiculous hack
  return unique_ptr<PathSet>(new PathSet(size, width_ul, width_ul * (int) method));
}

// TODO overload the method based on JoinMethod implementation
// not entirely sure what dynamic dispatch will to performance-wise, so doing everything statically for now
// JoinExec itself is very super not thread-safe
joined_res JoinExec::join(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const {

  // why no .clear() C++?
  while(!scores.empty())
    scores.pop();
  scores.push(Score());

  check_equal(uids.size(), paths0.size);
  check_true(paths_res.size == 0 || paths_res.size == uids.count_total_paths());
  for(const auto& uid : uids.uids)
    check_index(uid.location + uid.count - 1, paths1.size);

  double perm_scores[iterations];
  for(int k = 0; k < iterations; k++)
    perm_scores[k] = 0;
  this->perm_scores = perm_scores;

  if(method == Method::method1)
    return join_method1(uids, paths0, paths1, paths_res);
  if(method == Method::method2)
    return join_method2(uids, paths0, paths1, paths_res);
  throw std::logic_error("bad method");

}

joined_res JoinExec::format_result() const {

  while(scores.size() > top_k)
    scores.pop();

  joined_res res;
  res.permuted_scores.resize(iterations, 0);
  for(int k = 0; k < iterations; k++)
    res.permuted_scores[k] = perm_scores[k];
  res.scores.clear();
  while(!scores.empty()){
    res.scores.push_back(scores.top());
    scores.pop();
  }
  return res;

}

constexpr int prog_min_size = 4000;

static inline void prog_print(mutex& g_mutex, const int prog_size, size_t idx) {
  if(prog_size >= prog_min_size && idx % (prog_size / 10) == 0) {
    lock_guard<mutex> lock(g_mutex);
    printf(" %1.0lf%%", (idx * 100.0) / prog_size);
    std::cout << std::flush;
  }
}

joined_res JoinExec::join_method1(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const {

  // TODO not really size_t, should be limited to max total_paths size
  atomic<size_t> uid_idx(0);
  mutex g_mutex;

  if(uids.size() >= prog_min_size)
    printf("\nprogress:");

  // C++14 capture init syntax would be nice,
  // const auto exec = this;
  auto worker = [this, &uid_idx, &uids, &g_mutex, &paths0, &paths1, &paths_res](int tid) {

    const int iters = this->iterations;
    const int width_ul = this->width_ul;

    const bool keep_paths = paths_res.size != 0;

    uint64_t joined_block[width_ul];
    double perm_scores_block[iters];
    int perm_count_block[iters];

    // this mess actually makes a significant performance difference
    // don't know if there is a different way to stack-allocate a dynamic array in an instance field
    // or it may be a thread access thing... either way, this works
    JoinMethod1 method(this, paths1.size, perm_scores_block, joined_block, perm_count_block);

    const size_t prog_size = uids.size();

    size_t idx = -1;
    while((idx = uid_idx.fetch_add(1)) < uids.size()) {

      const auto& uid = uids[idx];

      prog_print(g_mutex, prog_size, idx);
      // if(prog_size >= prog_min_size && idx % (prog_size / 10) == 0) {
      //   lock_guard<mutex> lock(g_mutex);
      //   printf(" %1.0lf%%", (idx * 100.0) / prog_size);
      //   std::cout << std::flush;
      // }

      if(uid.count > 0) {

        // TODO long
        // start of result path block for this uid
        int path_idx = uid.path_idx;
        const uint64_t* path0 = paths0[idx];

        for(int loc_idx = 0, loc = uid.location; loc < uid.location + uid.count; loc_idx++, loc++) {
          method.score_permute_cpu(idx, loc, path0, paths1[loc]);
          if(keep_paths) {
            paths_res.set(path_idx, joined_block);
            path_idx += 1;
          }
        }

      }

    }

    // synchronize final section to merge results
    lock_guard<mutex> lock(g_mutex);
    method.drain_scores(this);
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

  if(uids.size() >= prog_min_size)
    printf(" - done!\n");

  return format_result();
}

// TODO check where path set sizes exceed max int
joined_res JoinExec::join_method2(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const {

  atomic<unsigned> uid_idx(0);
  mutex g_mutex;

  if(uids.size() >= prog_min_size)
    printf("\nprogress:");

  // C++14 capture init syntax would be nice,
  // const auto exec = this;
  auto worker = [this, &uid_idx, &uids, &g_mutex, &paths0, &paths1, &paths_res](int tid) {

    const int iters = this->iterations;
    const int width_ul = this->width_ul;

    const bool keep_paths = paths_res.size != 0;
    // TODO clean these up
    const int width_full_ul = width_ul * 2;
    const int p_block_size = iters * 2;

    double perm_scores_block[iters];
    uint64_t joined_block[width_ul * 2];
    int perm_pos_block[p_block_size];

    // this mess actually makes a significant performance difference
    // don't know if there is a different way to stack-allocate a dynamic array in an instance field
    // or it may be a thread access thing... either way, this works
    JoinMethod2 method(this, paths1.size, perm_scores_block, joined_block, perm_pos_block);

    const size_t prog_size = uids.size();

    int idx = -1;
    while((idx = uid_idx.fetch_add(1)) < uids.size()) {

      const auto& uid = uids[idx];

      prog_print(g_mutex, prog_size, idx);

      if(uid.count > 0) {

        // start of result path block for this uid
        int path_idx = uid.path_idx;

        const uint64_t* path_pos0 = paths0[idx];
        const uint64_t* path_neg0 = path_pos0 + width_ul;

        for(int loc_idx = 0, loc = uid.location; loc < uid.location + uid.count; loc_idx++, loc++) {
          const auto do_flip = uids.need_flip(idx, loc);
          const uint64_t* path_pos1 = (do_flip ? paths1[loc] : paths1[loc] + width_ul);
          const uint64_t* path_neg1 = (do_flip ? paths1[loc] + width_ul : paths1[loc]);
          method.score_permute_cpu(idx, loc, path_pos0, path_neg0, path_pos1, path_neg1);
          if(keep_paths) {
            paths_res.set(path_idx, joined_block);
            path_idx += 1;
          }
        }

      }

    }

    // synchronize final section to merge results
    lock_guard<mutex> lock(g_mutex);
    method.drain_scores(this);
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

  if(uids.size() >= prog_min_size)
    printf(" - done!\n");

  return format_result();
}

void JoinMethod1::score_permute_cpu(int idx, int loc, const uint64_t* path0, const uint64_t* path1) {

  // avoid de-refs like the plague
  const int top_k = exec->top_k;
  const int iters = exec->iterations;
  const int width_ul = exec->width_ul;
  const int flip_pivot_len = this->flip_pivot_len;

  int cases = 0;
  int ctrls = 0;

  // zero out join path holder
  memset(joined_block, 0, width_ul * sizeof(uint64_t));
  memset(perm_count_block, 0, iters * sizeof(int));

  uint64_t joined = 0;

  for(int k = 0; k < width_ul; k++) {

    joined = path0[k] | path1[k];

    if(joined != 0) {

      joined_block[k] = joined;

      cases += __builtin_popcountl(joined &  case_mask[k]);
      ctrls += __builtin_popcountl(joined & ~case_mask[k]);

      const uint64_t* p_mask = perm_case_mask + k;
      for(int r = 0; r < iters; r++)
        perm_count_block[r] += __builtin_popcountl(joined & p_mask[r * width_ul]);

    }

  }

  double score = value_table[cases][ctrls];
  double flips = value_table[ctrls][cases];

  if(score > scores.top().score)
    scores.push(Score(score, idx, loc));
  if(flips > scores.top().score)
    scores.push(Score(flips, idx, loc + flip_pivot_len));
  while(scores.size() > top_k)
    scores.pop();

  int total = cases + ctrls;
  for(int r = 0; r < iters; r++){
    int p_cases = perm_count_block[r];
    double p_score = value_table_max[p_cases][total - p_cases];
    if(p_score > perm_scores[r])
      perm_scores[r] = p_score;
  }

}

void JoinMethod2::score_permute_cpu(int idx, int loc, const uint64_t* path_pos0, const uint64_t* path_neg0, const uint64_t* path_pos1, const uint64_t* path_neg1) {

  // avoid de-refs like the plague
  const int top_k = exec->top_k;
  const int iters = exec->iterations;
  const int width_ul = exec->width_ul;
  const int flip_pivot_len = this->flip_pivot_len;

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

  for(int r = 0; r < iters; r++){
    perm_case_pos[r] = 0;
    perm_ctrl_pos[r] = 0;
  }

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

      if(true_pos != 0) {
        for(int r = 0; r < iters; r++) {
          const uint64_t* p_mask = perm_case_mask + r * width_ul;
          perm_case_pos[r] += __builtin_popcountl(true_pos & p_mask[k]);
        }
      }

      if(true_neg != 0) {
        for(int r = 0; r < iters; r++){
          const uint64_t* p_mask = perm_case_mask + r * width_ul;
          perm_ctrl_pos[r] += __builtin_popcountl(true_neg & p_mask[k]);
        }
      }

      joined_pos[k] = bit_pos;
      joined_neg[k] = bit_neg;

    }

  }

  int cases = case_pos + case_neg;
  int ctrls = ctrl_pos + ctrl_neg;

  double score = value_table[cases][ctrls];
  double flips = value_table[ctrls][cases];

  if(score > scores.top().score)
    scores.push(Score(score, idx, loc));
  if(flips > scores.top().score)
    scores.push(Score(flips, idx, loc + flip_pivot_len));
  while(scores.size() > top_k)
    scores.pop();

  for(int r = 0; r < iters; r++){
    int perm_cases = perm_case_pos[r] + (total_neg - perm_ctrl_pos[r]);
    int perm_ctrls = perm_ctrl_pos[r] + (total_pos - perm_case_pos[r]);

    double p_score = value_table_max[perm_cases][perm_ctrls];
    if(p_score > perm_scores[r])
      perm_scores[r] = p_score;
  }

}
