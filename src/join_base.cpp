
#include "gcre.h"
#include "methods.h"
#include <Rcpp.h>

// 'size' is the input element count, 'width' is bits needed for each element
// returns count that fits evenly into current vector size
static inline int pad_vector_size(int size, int width) {
  int bits = size * width;
  int nvec = (int) ceil(bits / (double) gs_vec_width);
  int vlen = nvec * gs_vec_width / width;
  return vlen;
}

static int vector_width(int num_cases, int num_ctrls) {
  return gs_vec_width * (int) ceil((num_cases + num_ctrls) / (double) gs_vec_width);
}

static int vector_width_cast(int num_cases, int num_ctrls, int width) {
  if(width > gs_vec_width)
    return 0;
  return vector_width(num_cases, num_ctrls) / width;
}

constexpr int prog_min_size = 4000;

static inline void prog_print(mutex& g_mutex, const int prog_size, st_uids_size idx) {
  if(prog_size >= prog_min_size && idx % (prog_size / 10) == 0) {
    lock_guard<mutex> lock(g_mutex);
    printf(" %1.0lf%%", (idx * 100.0) / prog_size);
    // Rcpp::Rcout << std::flush;
  }
}

// TODO uid and loc index integer sizes

JoinExec::JoinExec(string method_name, int num_cases, int num_ctrls, int iters) :

method(to_method(method_name)),
num_cases(num_cases),
num_ctrls(num_ctrls),
width_ul(vector_width_cast(num_cases, num_ctrls, 64)),
// need to store 32-bit ints and single-precision floats per iteration
iters_requested(iters),
iterations(pad_vector_size(iters, 32)) {

  check_true(num_cases > 0 && num_ctrls > 0 && iters >= 0);

  // create case mask from case/control ranges
  case_mask = new uint64_t[width_ul];
  for(int k = 0; k < width_ul; k++)
    case_mask[k] = bit_zero_ul;
  for(int k = 0; k < num_cases; k++)
    case_mask[k/64] |= bit_one_ul << k % 64;

  // printf("adjusting input to vector width: %d -> %d\n", num_cases + num_ctrls, vector_width(num_cases, num_ctrls));
  // printf("adjusting iterations to vector width: %d -> %d\n", iters_requested, iterations);
  // printf("[init exec] method: '%d' | data: %d x %d (width: %d ul)\n", method, num_cases, num_ctrls, width_ul);
}

// TODO floats for lookup
void JoinExec::setValueTable(const vec2d_d& table){

  // printf("received value table: %lu x %lu\n", table.size(), table.size() > 0 ? table.front().size() : 0);

  // indices can be flipped and counts could include all 0 to all cases/controls (so max index == size)
  size_t table_size = num_cases + num_ctrls + 1;

  // if(table.size() < table_size || table.size() == 0 || table.front().size() < table_size)
  //   printf("  ** WARN: value table smaller than %lu x %lu, will be padded\n", table_size, table_size);

  value_table = vec2d_d(table_size, vec_d(table_size, -1.0));

  for(int r = 0; r < min(table_size, table.size()); r++) {
    for(int c = 0; c < min(table_size, table[r].size()); c++) {
      value_table[r][c] = table[r][c];
    }
  }

}

// data comes in as perm x patient: 0=flipped, 1=unflipped
// here converted to 1=case, 0=control
// intermediate mask is inverted, to keep 0-padding of the vector tail from becoming mostly cases
void JoinExec::setPermutedCases(const vec2d_i& data) {

  // printf("loading permuted case masks: %lu x %lu (iterations: %d)\n", data.size(), data.size() > 0 ? data.front().size() : 0, iterations);

  if(iters_requested < data.size())
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
      perm_case_mask[k * iterations + r] = mask_k;
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

// create instance of join implementation based on specified method type
// TJoinMethod JoinExec::createMethod(const UidRelSet& uids, const int flip_pivot_len, float* p_perm_scores) const {
TJoinMethod JoinExec::createMethod(const UidRelSet& uids, float* p_perm_scores) const {
  if(method == Method::method1)
    return TJoinMethod(new JoinMethod1(this, uids, p_perm_scores));
  if(method == Method::method2)
    return TJoinMethod(new JoinMethod2(this, uids, p_perm_scores));
    // return TJoinMethod(new JoinMethod2(this, uids, flip_pivot_len, p_perm_scores));
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

// TODO width based on method needs
TPathSet JoinExec::createPathSet(st_pathset_size size) const {
  // std::make_unique is C++14
  // also the vlen size parameter is a ridiculous hack
  return TPathSet(new PathSet(size, width_ul, width_ul * (int) method));
}

template<typename Worker>
joined_res JoinExec::execute(const Worker& worker, bool show_progress) const {

  if(show_progress)
    printf("progress:");

  // use '0' to identify main thread
  if(max(0, nthreads) == 0) {
    worker(0);
  } else {
    vector<thread> pool(nthreads);
    for(int tid = 0; tid < pool.size(); tid++)
      pool[tid] = thread(worker, tid + 1);
    for(auto& th : pool)
      th.join();
  }

  if(show_progress)
    printf(" - done!\n\n");
    // printf(" - done!\n");

  return format_result();
}


// JoinExec is very super not thread-safe
joined_res JoinExec::join(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const {

  // why no .clear() C++?
  while(!scores.empty())
    scores.pop();
  scores.push(Score());

  check_equal(uids.size(), paths0.size);
  check_true(paths_res.size == 0 || paths_res.size == uids.count_total_paths());
  for(const auto& uid : uids.uids)
    check_index(uid.location + uid.count - 1, paths1.size);

  // again, funky allocation to get a stack array
  float perm_scores[iterations] ALIGNED;
  for(int k = 0; k < iterations; k++)
    perm_scores[k] = 0;
  this->perm_scores = perm_scores;

  atomic<st_uids_size> uid_idx(0);
  mutex g_mutex;

  // C++14 capture init syntax would be nice,
  // const auto exec = this;
  auto worker = [this, &uid_idx, &uids, &g_mutex, &paths0, &paths1, &paths_res](int tid) {

    const int iters = this->iterations;
    const int width_ul = this->width_ul;
    const bool keep_paths = paths_res.size != 0;

    // this mess actually makes a significant performance difference
    // don't know if there is a different way to stack-allocate a dynamic array in an instance field
    // or it may be a thread access thing... either way, this works
    float perm_scores_block[iters] ALIGNED;
    // TJoinMethod method = createMethod(uids, paths1.size, perm_scores_block);
    TJoinMethod method = createMethod(uids, perm_scores_block);

    uint64_t path_res[paths_res.vlen] ALIGNED;

    const st_uids_size prog_size = uids.size();
    st_uids_size idx = -1;
    while((idx = uid_idx.fetch_add(1)) < uids.size()) {

      const auto& uid = uids[idx];

      prog_print(g_mutex, prog_size, idx);

      if(uid.count > 0) {

        // start of result path block for this uid
        st_path_count path_idx = uid.path_idx;
        const uint64_t* path0 = paths0[idx];

        for(int loc_idx = 0, loc = uid.location; loc < uid.location + uid.count; loc_idx++, loc++) {
          if(keep_paths)
            memset(path_res, 0, paths_res.vlen * 8);
          method->score_permute(idx, loc, path0, paths1[loc], path_res, keep_paths);
          if(keep_paths) {
            paths_res.set(path_idx, path_res);
            path_idx += 1;
          }
        }

      }

    }

    // synchronize final section to merge results
    lock_guard<mutex> lock(g_mutex);
    method->merge_scores();
  };

  // return execute(worker, uids.size() >= prog_min_size);
  return execute(worker, true);

}
