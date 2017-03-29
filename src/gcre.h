#ifndef GCRE_H
#define GCRE_H

#include <cstring>
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <queue>
#include <iostream>

// #if defined __AVX2__
// #define SCORE_METHOD_NAME score_permute_avx2
// #define SCORE_IMPL_LABEL "AVX2"
// constexpr int gs_vec_width    = 256;
// #elif defined __AVX__
// #define SCORE_METHOD_NAME score_permute_sse4
// #define SCORE_IMPL_LABEL "SSE4"
// constexpr int gs_vec_width    = 256;
// #else
// #define SCORE_METHOD_NAME score_permute_sse2
// #define SCORE_IMPL_LABEL "SSE2"
// constexpr int gs_vec_width    = 128;
// #endif

#define SCORE_METHOD_NAME score_permute_cpu
#define SCORE_IMPL_LABEL "CPU"
constexpr int gs_vec_width    = 64;

constexpr int gs_vec_width_b  = gs_vec_width / 8;
constexpr int gs_vec_width_32 = gs_vec_width / 32;
constexpr int gs_vec_width_ul = gs_vec_width / 64;

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

  unsigned count_total_paths() const {
    unsigned total = 0;
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
    block = unique_ptr<uint64_t[]>(new uint64_t[size * vlen]);
    for(int k = 0; k < size * vlen; k++)
      block[k] = bit_zero_ul;
    printf("[%p] allocated block for path set: %d x %d ul, total width: %d (%'1lu bytes)\n", this, size, width_ul, vlen, size * vlen * sizeof(uint64_t));
  }

  inline const uint64_t* operator[](int idx) const {
    return block.get() + idx * vlen;
  }

  inline void set(int idx, const uint64_t* data) {
    memcpy(block.get() + idx * vlen, data, vlen * sizeof(uint64_t));
  }

  // positives are loaded into the first half of the record, so m1/m2 loading is the same
  void load(const vec2d_i& data) {

    printf("loading data to path set: %lu x %lu\n", data.size(), data.size() > 0 ? data.front().size() : 0);

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
    printf("total set: %d (input: %d)\n", count_set, count_in);
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
    printf("copied path set by index selection, indices: %lu; total set: %d\n", indices.size(), count_set);
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
  const int iterations;
  
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

  void setValueTable(vec2d_d table);

  void setPermutedCases(const vec2d_i& perm_cases);

  unique_ptr<PathSet> createPathSet(int size) const;
  
  joined_res join(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const;

protected:

  mutable priority_queue<Score> scores;
  mutable double* perm_scores;

  joined_res format_result() const;
  joined_res join_method1(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const;
  joined_res join_method2(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const;

  static int vector_width_ul(int num_cases, int num_ctrls) {
    int width = (int) ceil((num_cases + num_ctrls) / (double) gs_vec_width);
    printf("adjusting input to vector width: %d -> %d\n", num_cases + num_ctrls, width * gs_vec_width);
    return width * gs_vec_width_ul;
  }

  static int iter_size(int iters) {
    int width = gs_vec_width_32 * (int) ceil(iters / (double) gs_vec_width_32);
    printf("adjusting iterations to vector width: %d -> %d\n", iters, width);
    return width;
  }

  vec2d_d value_table;
  vec2d_d value_table_max;
  uint64_t* case_mask = nullptr;
  uint64_t* perm_case_mask = nullptr;

};

#endif
