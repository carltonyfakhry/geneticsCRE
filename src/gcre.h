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
#define gs_align_size 8
#elif defined __AVX512F__
constexpr const int gs_vec_width = 512;
const std::string gs_instr_label = "AVX-512";
#define gs_align_size 64
#elif defined __AVX2__
constexpr int gs_vec_width = 256;
const std::string gs_instr_label = "AVX2";
#define gs_align_size 32
#elif defined __SSE4_2__
constexpr int gs_vec_width = 256;
const std::string gs_instr_label = "SSE4";
#define gs_align_size 32
#else
constexpr int gs_vec_width = 128;
const std::string gs_instr_label = "SSE2";
#define gs_align_size 16
#endif

//constexpr int gs_align_size = gs_vec_width / 8;
#define ALIGNED __attribute__ ((aligned (gs_align_size)))

#include "gcre_types.h"
#include "gcre_paths.h"

using namespace std;

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

  st_path_count count_total_paths() const {
    st_path_count total = 0;
    for(const auto& uid : uids)
      total += uid.count;
    return total;
  }

};

class JoinMethod {

public:

  virtual void score_permute(int idx, int loc, const uint64_t* path0, const uint64_t* path1, uint64_t* path_res, bool keep_paths) = 0;
  virtual void merge_scores() = 0;

};

using TJoinMethod = unique_ptr<JoinMethod>;

class JoinExec {

  // TODO
  friend class JoinMethod;
  friend class JoinMethod_Base;
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
    // if(name == "method2")
    else
      return Method::method2;
    // cout << "unknown method: " << name << endl;
    // exit(1);
  }

  JoinExec(string method_name, int num_cases, int num_ctrls, int iters);

  ~JoinExec(){
    if(case_mask)
      delete[] case_mask;
    if(perm_case_mask)
      delete[] perm_case_mask;
  }

  void print_vector_info() {
    // printf("\n");
    // printf("########################\n");
    // printf("  instruction set: %s\n", gs_instr_label.c_str());
    // printf("     vector width: %u\n", gs_vec_width);
    // printf("        alignment: %u\n", gs_align_size);
    // printf("########################\n");
    // printf("\n");
  }

  void setValueTable(const vec2d_d& table);

  void setPermutedCases(const vec2d_i& perm_cases);

  // TJoinMethod createMethod(const UidRelSet& uids, const int flip_pivot_len, float* p_perm_scores) const;
  TJoinMethod createMethod(const UidRelSet& uids, float* p_perm_scores) const;

  TPathSet createPathSet(st_pathset_size size) const;

  joined_res join(const UidRelSet& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const;

protected:

  mutable priority_queue<Score> scores;
  mutable float* perm_scores;

  joined_res format_result() const;

  vec2d_d value_table;
  uint64_t* case_mask = nullptr;
  uint64_t* perm_case_mask = nullptr;

  // TODO constness isn't really true here
  template<typename Worker>
  joined_res execute(const Worker& worker, bool show_progress) const;

};

#endif
