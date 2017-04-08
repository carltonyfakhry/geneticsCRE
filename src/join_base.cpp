
#include "gcre.h"

using namespace std;

class JoinMethod {

public:

  priority_queue<Score> scores;
  float* perm_scores;

  JoinMethod(const JoinExec* exec, const int flip_pivot_len, float* p_perm_scores, uint64_t* p_joined_block, int* p_perm_count_block) :

  exec(exec),
  top_k(exec->top_k),
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
  const int top_k;
  const int flip_pivot_len;

  const uint64_t* case_mask;
  const uint64_t* perm_case_mask;

  const vec2d_d& value_table;
  const vec2d_d& value_table_max;

  uint64_t* joined_block;
  int* perm_count_block;

  void keep_score(int idx, int loc, int cases, int ctrls) {

    double score = value_table[cases][ctrls];
    if(score > scores.top().score)
      scores.push(Score(score, idx, loc));
    
    double flips = value_table[ctrls][cases];
    if(flips > scores.top().score)
      scores.push(Score(flips, idx, loc + flip_pivot_len));
    while(scores.size() > top_k)
      scores.pop();

  }

};

class JoinMethod1 : public JoinMethod {

public:

  JoinMethod1(const JoinExec* exec, const int flip_pivot_len, float* p_perm_scores, uint64_t* p_joined_block, int* p_perm_count_block) : JoinMethod(exec, flip_pivot_len, p_perm_scores, p_joined_block, p_perm_count_block) {}

  void score_permute_cpu(int idx, int loc, const uint64_t* path0, const uint64_t* path1);
  void score_permute_sse2(int idx, int loc, const uint64_t* path0, const uint64_t* path1);
  void score_permute_sse4(int idx, int loc, const uint64_t* path0, const uint64_t* path1);
  void score_permute_avx2(int idx, int loc, const uint64_t* path0, const uint64_t* path1);

};

class JoinMethod2 : public JoinMethod {

public:

  JoinMethod2(const JoinExec* exec, const int flip_pivot_len, float* p_perm_scores, uint64_t* p_joined_block, int* p_perm_count_block) : JoinMethod(exec, flip_pivot_len, p_perm_scores, p_joined_block, p_perm_count_block) {

    joined_pos = joined_block;
    joined_neg = joined_block + exec->width_ul;

    perm_case_pos = perm_count_block;
    perm_ctrl_pos = perm_count_block + exec->iterations;

  }

  void score_permute_cpu(int idx, int loc, const uint64_t* path0_pos, const uint64_t* path0_neg, const uint64_t* path1_pos, const uint64_t* path1_neg);
  void score_permute_sse2(int idx, int loc, const uint64_t* path0_pos, const uint64_t* path0_neg, const uint64_t* path1_pos, const uint64_t* path1_neg);
  void score_permute_sse4(int idx, int loc, const uint64_t* path0_pos, const uint64_t* path0_neg, const uint64_t* path1_pos, const uint64_t* path1_neg);
  void score_permute_avx2(int idx, int loc, const uint64_t* path0_pos, const uint64_t* path0_neg, const uint64_t* path1_pos, const uint64_t* path1_neg);

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
width_ul(vector_width_cast(num_cases, num_ctrls, 64)),
width_dw(vector_width_cast(num_cases, num_ctrls, 32)),
width_qw(vector_width_cast(num_cases, num_ctrls, 64)),
width_dq(vector_width_cast(num_cases, num_ctrls, 128)),
width_qq(vector_width_cast(num_cases, num_ctrls, 256)),
width_oq(vector_width_cast(num_cases, num_ctrls, 512)),
iters_requested(iters),
iterations(iter_size_dw(iters)) {

  check_true(num_cases > 0 && num_ctrls > 0 && iters > 0);
  width_vec = width_ul / gs_vec_width_qw;

  // create case mask from case/control ranges
  case_mask = new uint64_t[width_ul];
  for(int k = 0; k < width_ul; k++)
    case_mask[k] = bit_zero_ul;
  for(int k = 0; k < num_cases; k++)
    case_mask[k/64] |= bit_one_ul << k % 64;

  printf("adjusting input to vector width: %d -> %d\n", num_cases + num_ctrls, vector_width(num_cases, num_ctrls));
  printf("adjusting iterations to vector width: %d -> %d\n", iters_requested, iterations);
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

  // when iterations > iters_requested, the tails are zeroed out
  perm_case_mask = new uint64_t[iterations * width_ul];
  for(int k = 0; k < iterations * width_ul; k++)
    perm_case_mask[k] = bit_zero_ul;

  for(int r = 0; r < min( (size_t) iters_requested, data.size()); r++) {
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
      // group by slice
      perm_case_mask[k * iterations + r] = mask_k;
      // group by iteration
      // perm_case_mask[r * width_ul + k] = mask_k;
      count_cases += __builtin_popcountl(mask_k);
    }
    if(count_cases != num_cases)
      printf("  ** WARN perm mask did not set correct number of cases: %d != %d (mask %04d)\n", count_cases, num_cases, r);
  }

  if(iters_requested > data.size()) {
    printf("  ** WARN not enough permuted cases, some will be reused to match iterations - hope this is for testing!\n");
    for(int r = data.size(); r < iters_requested; r++) {
      int c = r % data.size();
      for(int k = 0; k < width_ul; k++) 
        perm_case_mask[k * iterations + r] = perm_case_mask[k * iterations + c];
    }
  }

}

// TODO width based on method needs
unique_ptr<PathSet> JoinExec::createPathSet(int size) const {
  // std::make_unique is C++14
  // also the size parameter is a ridiculous hack
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

  float perm_scores[iterations];

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
  res.permuted_scores.resize(iters_requested, 0);
  for(int k = 0; k < iters_requested; k++)
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

  atomic<size_t> uid_idx(0);
  mutex g_mutex;

  // C++14 capture init syntax would be nice,
  // const auto exec = this;
  auto worker = [this, &uid_idx, &uids, &g_mutex, &paths0, &paths1, &paths_res](int tid) {

    const int iters = this->iterations;
    const int width_ul = this->width_ul;

    const bool keep_paths = paths_res.size != 0;

    // TODO method implementation should request size

    uint64_t joined_block[width_ul] __attribute__ ((aligned (gs_align_size)));
    float perm_scores_block[iters] __attribute__ ((aligned (gs_align_size)));
    int perm_count_block[iters] __attribute__ ((aligned (gs_align_size)));

    // this mess actually makes a significant performance difference
    // don't know if there is a different way to stack-allocate a dynamic array in an instance field
    // or it may be a thread access thing... either way, this works
    JoinMethod1 method(this, paths1.size, perm_scores_block, joined_block, perm_count_block);

    const size_t prog_size = uids.size();

    size_t idx = -1;
    while((idx = uid_idx.fetch_add(1)) < uids.size()) {

      const auto& uid = uids[idx];

      prog_print(g_mutex, prog_size, idx);

      if(uid.count > 0) {

        // TODO long
        // start of result path block for this uid
        int path_idx = uid.path_idx;
        const uint64_t* path0 = paths0[idx];

        for(int loc_idx = 0, loc = uid.location; loc < uid.location + uid.count; loc_idx++, loc++) {
          method.SCORE_METHOD_NAME(idx, loc, path0, paths1[loc]);
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

  return execute(worker, uids.size() >= prog_min_size);
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

    float perm_scores_block[iters] __attribute__ ((aligned (gs_align_size)));
    uint64_t joined_block[width_ul * 2] __attribute__ ((aligned (gs_align_size)));
    int perm_pos_block[p_block_size] __attribute__ ((aligned (gs_align_size)));

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
          method.SCORE_METHOD_NAME(idx, loc, path_pos0, path_neg0, path_pos1, path_neg1);
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

  return execute(worker, uids.size() >= prog_min_size);
}

// .cxx extensions so that these don't get picked up by Rcpps' makefile (also that they shouldn't be compiled independently)

#ifdef COMPILE_CPU
#include "impl_cpu.cxx"
#endif

#ifdef COMPILE_SSE2
#include "impl_sse2b.cxx"
#endif

#ifdef COMPILE_SSE4
#include "impl_sse4.cxx"
#endif

#ifdef COMPILE_AVX2
#include "impl_avx2.cxx"
#endif
