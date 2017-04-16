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

constexpr int gs_vec_width    = 64;

constexpr int gs_vec_width_b  = gs_vec_width / 8;
constexpr int gs_vec_width_32 = gs_vec_width / 32;
constexpr int gs_vec_width_ul = gs_vec_width / 64;

using namespace std;

const uint64_t bit_zero_ul = 0;
const uint64_t bit_one_ul = 1;

// types that can hold max allowed values for various entities
using st_pathset_size = uint32_t;
using st_total_paths = uint64_t;

// TODO apparently aliases are now the thing?
typedef std::vector<int> vec_i;
typedef std::vector<double> vec_d;
typedef std::vector<uint64_t> vec_u64;

typedef std::vector<std::vector<double>> vec2d_d;
typedef std::vector<std::vector<int>> vec2d_i;
typedef std::vector<std::vector<int8_t>> vec2d_i8;
typedef std::vector<std::vector<uint16_t>> vec2d_u16;
typedef std::vector<std::vector<uint64_t>> vec2d_u64;

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

  st_total_paths count_total_paths() const {
    st_total_paths total = 0;
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

  const st_pathset_size size;
  const uint16_t width_ul;
  const uint16_t vlen;

  // TODO abstract interface for other storage types
  PathSet(st_pathset_size size, const int width_ul, const uint16_t vlen) : size(size), width_ul(width_ul), vlen(vlen) {
    const size_t min_block_size = 16 * 1024 * 1024;
    const size_t block_size = size * vlen * sizeof(uint64_t);
    if(block_size > min_block_size)
      printf("allocate block for path-set (%5d x %d ul) %'15lu bytes: ", size, vlen, block_size);
    block = unique_ptr<uint64_t[]>(new uint64_t[size * vlen]);
    for(int k = 0; k < size * vlen; k++)
      block[k] = bit_zero_ul;
    if(block_size > min_block_size)
      printf("[%p]\n", block.get());
  }

  inline const uint64_t* operator[](st_pathset_size idx) const {
    check_index(idx, size);
    return block.get() + idx * vlen;
  }

  inline void set(st_pathset_size idx, const uint64_t* data) {
    check_index(idx, size);
    memcpy(block.get() + idx * vlen, data, vlen * sizeof(uint64_t));
  }

  // TODO check input size dimensions
  // positives are loaded into the first half of the record, so m1/m2 loading is the same
  void load(const vec2d_i& data) {

    printf("loading data to path-set (%lu x %lu): ", data.size(), data.size() > 0 ? data.front().size() : 0);

    check_true(size == data.size());
    long count_in = 0;
    for(auto r = 0; r < data.size(); r++) {
      check_index(data[r].size(), width_ul * 64);
      for(auto c = 0; c < data[r].size(); c++) {
        if(data[r][c] != 0) {
          count_in += 1;
          block[r * vlen + c/64] |= bit_one_ul << c % 64;
        }
      }
    }

    long count_set = 0;
    for(auto k = 0; k < size * vlen; k++)
      count_set += __builtin_popcountl(block[k]);
    printf("%ld (of %ld)\n", count_set, count_in);
    check_true(count_set == count_in);
  }

  // create new path set from provided indices
  // this is only done for the initial path sets, and the index values come from R, so int instead of st_pathset_size
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

  unique_ptr<PathSet> createPathSet(st_pathset_size size) const;
  
  joined_res join(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const;

protected:

  mutable priority_queue<Score> scores;
  mutable double* perm_scores;

  joined_res format_result() const;
  joined_res join_method1(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const;
  joined_res join_method2(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const;

  static int vector_width_ul(int num_cases, int num_ctrls) {
    check_true(num_cases > 0 && num_ctrls > 0);
    int width = (int) ceil((num_cases + num_ctrls) / (double) gs_vec_width);
    printf("adjusting input to vector width: %d -> %d\n", num_cases + num_ctrls, width * gs_vec_width);
    return width * gs_vec_width_ul;
  }

  static int iter_size(int iters) {
    check_true(iters > 0);
    int width = gs_vec_width_32 * (int) ceil(iters / (double) gs_vec_width_32);
    printf("adjusting iterations to vector width: %d -> %d\n", iters, width);
    return width;
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

  Timer(const JoinExec& exec, int path_length, st_total_paths total_paths) : exec(exec), path_length(path_length), total_paths(total_paths) {}

  ~Timer(){
    auto time = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - start);
    printf("\nTIME:%d CPU %d %d %d %lu %d %d %lu\n\n", getpid(), exec.method, exec.width_ul * 64, path_length, total_paths, exec.iterations, exec.nthreads, (unsigned long) time.count());
  }

private:

  const chrono::system_clock::time_point start = chrono::system_clock::now();
  const JoinExec& exec;
  const int path_length;
  const st_total_paths total_paths;

};

#endif
