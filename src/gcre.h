#ifndef GCRE_H
#define GCRE_H

#include <cstring>
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <limits>
#include <chrono>
#include <memory>
#include <vector>
#include <queue>
#include <iostream>

#include <unistd.h>

// will not compile without sse2, but cpu can be enable manually for testing
#ifdef COMPILE_CPU

#define SCORE_METHOD_NAME score_permute_cpu
constexpr int gs_vec_width = 64;
const std::string gs_impl_label = "CPU";

#else

#if defined __AVX512__

#define COMPILE_EVEX
#define SCORE_METHOD_NAME score_permute_evex
constexpr int gs_vec_width = 512;
const std::string gs_impl_label = "AVX-512";

#elif defined __AVX2__

#define COMPILE_AVX2
#define SCORE_METHOD_NAME score_permute_avx2
constexpr int gs_vec_width = 256;
const std::string gs_impl_label = "AVX2";

// SSE4 also requires AVX (don't think that affects many CPUs)
#elif defined __AVX__

#define COMPILE_SSE4
#define SCORE_METHOD_NAME score_permute_sse4
constexpr int gs_vec_width = 256;
const std::string gs_impl_label = "SSE4";

#else

#define COMPILE_SSE2
#define SCORE_METHOD_NAME score_permute_sse2
constexpr int gs_vec_width = 128;
const std::string gs_impl_label = "SSE2";

#endif
#endif

constexpr int gs_align_size = gs_vec_width / 8;

constexpr int gs_vec_width_b  = gs_vec_width / 8;
constexpr int gs_vec_width_dw = gs_vec_width / 32;
constexpr int gs_vec_width_qw = gs_vec_width / 64;
constexpr int gs_vec_width_dq = gs_vec_width / 128;
constexpr int gs_vec_width_qq = gs_vec_width / 256;
constexpr int gs_vec_width_oq = gs_vec_width / 512;

using namespace std;

const uint64_t bit_zero_ul = 0;
const uint64_t bit_one_ul = 1;

// TODO apparently aliases are now the thing?
typedef std::vector<int> vec_i;
typedef std::vector<double> vec_d;
typedef std::vector<uint64_t> vec_u64;

typedef std::vector<std::vector<double>> vec2d_d;
typedef std::vector<std::vector<int>> vec2d_i;
typedef std::vector<std::vector<int8_t>> vec2d_i8;
typedef std::vector<std::vector<uint16_t>> vec2d_u16;
typedef std::vector<std::vector<uint64_t>> vec2d_u64;

enum class Method { method1 = 1, method2 = 2 };

class Score {
public:
  double score = -numeric_limits<double>::infinity();
  int src = -1;
  int trg = -1;
  inline Score() {}
  inline Score(double score, int src, int trg) : score(score), src(src), trg(trg) {}
  // reverse sort for priority queue
  friend bool operator<(Score a, Score b) { return a.score > b.score; }
};

struct uid_ref {
  int src;
  int trg;
  int count;
  int location;
  int path_idx;
};

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

struct joined_res {
  vector<Score> scores;
  vec_d permuted_scores;
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
    return block.get() + idx * vlen;
  }

  inline void set(int idx, const uint64_t* data) {
    memcpy(block.get() + idx * vlen, data, vlen * sizeof(uint64_t));
  }

  // TODO check input size dimensions
  // positives are loaded into the first half of the record, so m1/m2 loading is the same
  void load(const vec2d_i& data) {

    printf("loading data to path-set (%lu x %lu): ", data.size(), data.size() > 0 ? data.front().size() : 0);

    int count_in = 0;
    for(int r = 0; r < data.size(); r++) {
      for(int c = 0; c < data[r].size(); c++) {
        if(data[r][c] != 0) {
          block[r * vlen + c/64] |= bit_one_ul << c % 64;
          count_in += 1;
        }
      }
    }

    int count_set = 0;
    for(int k = 0; k < size * vlen; k++)
      count_set += __builtin_popcountl(block[k]);

    printf("%d (of %d)\n", count_set, count_in);
  }

  // create new path set from provided indices
  unique_ptr<PathSet> select(const vector<int>& indices) const {
    unique_ptr<PathSet> pset(new PathSet(indices.size(), width_ul, vlen));
    for(int k = 0; k < pset->size; k++) {
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

public:

  const Method method;
  const int num_cases;
  const int num_ctrls;
  const int width_ul;
  const int width_dw;
  const int width_qw;
  const int width_dq;
  const int width_qq;
  const int width_oq;
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
    printf("  vector sizes with: %s\n", gs_impl_label.c_str());
    printf("########################\n");
    printf("\n    aligned to: %u\n", gs_align_size);
    printf("\n  vector width: %u\n", gs_vec_width);
    printf("    bytes: %u\n", gs_vec_width_b);
    printf("    dword: %u\n", gs_vec_width_dw);
    printf("    qword: %u\n", gs_vec_width_qw);
    printf("    dquad: %u\n", gs_vec_width_dq);
    printf("    qquad: %u\n", gs_vec_width_qq);
    printf("    oquad: %u\n", gs_vec_width_oq);
    printf("\n   input width: %d (%d)\n", width_ul * 64, num_cases + num_ctrls);
    printf("    dword: %u\n", width_dw);
    printf("    qword: %u\n", width_qw);
    printf("    dquad: %u\n", width_dq);
    printf("    qquad: %u\n", width_qq);
    printf("    oquad: %u\n", width_oq);
    printf("\n########################\n");
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

  static int iter_size_dw(int iters) {
    return gs_vec_width_dw * (int) ceil(iters / (double) gs_vec_width_dw);
  }

  vec2d_d value_table;
  vec2d_d value_table_max;
  uint64_t* case_mask = nullptr;
  uint64_t* perm_case_mask = nullptr;

};

class Timer {

public:

  static void print_header() {
    printf("\nTIME:PID IMPL METHOD WIDTH LENGTH PATHS PERMS THREADS MS\n\n");
  }

  Timer(const JoinExec& exec, int path_length, uint64_t total_paths) : exec(exec), path_length(path_length), total_paths(total_paths) {}

  ~Timer(){
    auto time = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - start);
    printf("\nTIME:%d %s %d %d %d %lu %d %d %lu\n\n", getpid(), gs_impl_label.c_str(), exec.method, exec.width_ul * 64, path_length, total_paths, exec.iterations, exec.nthreads, (unsigned long) time.count());
  }

private:

  const chrono::system_clock::time_point start = chrono::system_clock::now();
  const JoinExec& exec;
  const int path_length;
  const uint64_t total_paths;

};

#endif
