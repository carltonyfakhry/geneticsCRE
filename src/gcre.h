#ifndef GCRE_H
#define GCRE_H

#include <cstring>
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <memory>
#include <vector>
#include <queue>
#include <iostream>

#include <thread>
#include <atomic>
#include <mutex>


// will not compile without sse2 (it's a clang thing), but cpu-only can be enabled manually for testing
#ifdef COMPILE_CPU
constexpr int gs_vec_width = 64;
const std::string gs_instr_label = "x86-64";
#elif defined __AVX512__
constexpr int gs_vec_width = 512;
const std::string gs_instr_label = "AVX-512";
#elif defined __AVX2__
constexpr int gs_vec_width = 256;
const std::string gs_instr_label = "AVX2";
#elif defined __SSE4_2__
constexpr int gs_vec_width = 256;
const std::string gs_instr_label = "SSE4";
#else
constexpr int gs_vec_width = 128;
const std::string gs_instr_label = "SSE2";
#endif

constexpr int gs_align_size = gs_vec_width / 8;
#define ALIGNED __attribute__ ((aligned (gs_align_size)))

#include "gcre_types.h"

// 'size' is the input element count, 'width' is bits needed for each element
// returns count that fits evenly into current vector size
inline int pad_vector_size(int size, int width) {
  int bits = size * width;
  int nvec = (int) ceil(bits / (double) gs_vec_width);
  int vlen = nvec * gs_vec_width / width;
  return vlen;
}

using namespace std;

inline void check_true(bool condition) {
  if(!condition)
    throw std::logic_error("assertion");
}

inline void check_equal(size_t one, size_t two) {
  if(one != two)
    throw std::logic_error("assertion");
}

inline void check_index(long value, size_t size) {
  if(value < 0 || value >= size)
    throw std::out_of_range("assertion");
}

inline void check_range(long value, long min, long max) {
  if(value < min || value > max)
    throw std::out_of_range("assertion");
}


// TODO need to check memory use and performance for larger lengths
class UidRelSet {

public:

  const int path_length;
  const vector<uid_ref> uids;
  const vector<int> signs;

  UidRelSet(int path_length, vector<uid_ref> uids, vector<int> signs) : path_length(path_length), uids(uids), signs(signs) {
  }

  size_t size() const {
    return uids.size();
  }

  // TODO will this copy?
  inline const uid_ref& operator[](int idx) const {
    check_index(idx, uids.size());
    return uids[idx];
  }

  // TODO this is probably slow - need a short-circuit for path_length != 3
  inline bool need_flip(int idx, int loc) const {
    // TODO check correctness (path_length >5?)
    int sign = 0;
    if(path_length > 3)
      sign = signs[idx];
    else if(path_length < 3)
      sign = signs[loc];
    else if(path_length == 3)
      sign = (signs[idx] + signs[loc] == 0) ? -1 : 1;
    return sign == 1;
  }

  uint64_t count_total_paths() const {
    uint64_t total = 0;
    for(const auto& uid : uids)
      total += uid.count;
    return total;
  }

};

class PathSet {

public:

  const int size;
  const int width_ul;
  const int vlen;

  // TODO abstract interface for other storage types
  PathSet(const int size, const int width_ul, const int vlen) : size(size), width_ul(width_ul), vlen(vlen) {
    const size_t min_block_size = 16 * 1024 * 1024;
    const size_t block_size = size * vlen * sizeof(uint64_t);
    if(block_size > min_block_size)
      printf("allocate block for path-set (%5d x %d ul) %'15lu bytes: ", size, vlen, block_size);
    block = unique_ptr<uint64_t[]>((uint64_t*) aligned_alloc(gs_align_size, sizeof(uint64_t) * size * vlen));
    for(int k = 0; k < size * vlen; k++)
      block[k] = bit_zero_ul;
    if(block_size > min_block_size)
      printf("[%p]\n", block.get());
  }

  inline const uint64_t* operator[](int idx) const {
    check_index(idx, size);
    return block.get() + idx * vlen;
  }

  inline void set(int idx, const uint64_t* data) {
    check_index(idx, size);
    memcpy(block.get() + idx * vlen, data, vlen * sizeof(uint64_t));
  }

  // TODO check input size dimensions
  // positives are loaded into the first half of the record, so m1/m2 loading is the same
  void load(const vec2d_i& data) {

    printf("loading data to path-set (%lu x %lu): ", data.size(), data.size() > 0 ? data.front().size() : 0);

    check_true(size == data.size());
    int count_in = 0;
    for(int r = 0; r < data.size(); r++) {
      check_index(data[r].size(), width_ul * 64);
      for(int c = 0; c < data[r].size(); c++) {
        if(data[r][c] != 0) {
          count_in += 1;
          block[r * vlen + c/64] |= bit_one_ul << c % 64;
        }
      }
    }

    int count_set = 0;
    for(int k = 0; k < size * vlen; k++)
      count_set += __builtin_popcountl(block[k]);

    printf("%d (of %d)\n", count_set, count_in);
    check_true(count_set == count_in);
  }

  // create new path set from provided indices
  unique_ptr<PathSet> select(const vector<int>& indices) const {
    unique_ptr<PathSet> pset(new PathSet(indices.size(), width_ul, vlen));
    for(int k = 0; k < pset->size; k++) {
      check_index(indices[k], size);
      memcpy(pset->block.get() + vlen * k, block.get() + vlen * indices[k], vlen * sizeof(uint64_t));
    }
    int count_set = 0;
    for(int k = 0; k < pset->size * vlen; k++)
      count_set += __builtin_popcountl(pset->block[k]);
    // printf("copied path set by index selection, indices: %lu; total set: %d\n", indices.size(), count_set);
    return pset;
  }

protected:

  unique_ptr<uint64_t[]> block;

};

class JoinExec {

  friend class JoinMethod;
  friend class JoinMethod1;
  friend class JoinMethod2;

public:

  const Method method;
  const int num_cases;
  const int num_ctrls;
  const int width_ul;
  const int iterations;
  const int iters_requested;
  
  int top_k = 12;
  int nthreads = 0;

  int width_vec;

  static inline Method to_method(string name) {
    if(name == "method1")
      return Method::method1;
    if(name == "method2")
      return Method::method2;
    cout << "unknown method: " << name << endl;
    exit(1);
  }

  JoinExec(string method_name, int num_cases, int num_ctrls, int iters);

  ~JoinExec(){
    if(case_mask)
      delete[] case_mask;
    if(perm_case_mask)
      delete[] perm_case_mask;
  }

  void print_vector_info() {
    printf("\n");
    printf("########################\n");
    printf("  instruction set: %s\n", gs_instr_label.c_str());
    printf("     vector width: %u\n", gs_vec_width);
    printf("        alignment: %u\n", gs_align_size);
    printf("########################\n");
    printf("\n");
  }

  void setValueTable(vec2d_d table);

  void setPermutedCases(const vec2d_i& perm_cases);

  unique_ptr<PathSet> createPathSet(int size) const;
  
  joined_res join(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const;

protected:

  mutable priority_queue<Score> scores;
  mutable float* perm_scores;

  joined_res format_result() const;
  joined_res join_method1(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const;
  joined_res join_method2(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const;

  static int vector_width(int num_cases, int num_ctrls) {
    return gs_vec_width * (int) ceil((num_cases + num_ctrls) / (double) gs_vec_width);
  }

  // width larger than current will not be set
  static int vector_width_cast(int num_cases, int num_ctrls, int width) {
    if(width > gs_vec_width)
      return 0;
    return vector_width(num_cases, num_ctrls) / width;
  }

  vec2d_d value_table;
  vec2d_d value_table_max;
  uint64_t* case_mask = nullptr;
  uint64_t* perm_case_mask = nullptr;

  // TODO constness isn't really true here
  template<typename Worker>
  joined_res execute(const Worker& worker, bool show_progress) const {

    if(show_progress)
      printf("\nprogress:");

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
      printf(" - done!\n");

    return format_result();
  }


};

#endif
