#ifndef GCRE_H
#define GCRE_H

#include <inttypes.h>
#include <stdio.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <cmath>
#include <limits>
#include <vector>
#include <set>
#include <queue>
#include <cstdio>
#include <iostream>
#include <fstream>

// these are types visible to R; keeping that minimal

using namespace std;

const uint64_t bit_zero_ul = 0;
const uint64_t bit_one_ul = 1;

typedef std::vector<int> vec_i;
typedef std::vector<double> vec_d;
typedef std::vector<uint64_t> vec_u64;

typedef std::vector<std::vector<double>> vec2d_d;
typedef std::vector<std::vector<int>> vec2d_i;
typedef std::vector<std::vector<int8_t>> vec2d_i8;
typedef std::vector<std::vector<uint16_t>> vec2d_u16;
typedef std::vector<std::vector<uint64_t>> vec2d_u64;

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

struct paths_type {};

struct paths_base : paths_type {
  int size = 0;
  int width_ul = 0;
  int num_cases = 0;
};

struct join_config {
  int num_cases = 0;
  int num_controls = 0;
  int top_k = 0;
  int path_length = 0;
  int iterations = 0;
  int nthreads = 0;
};

struct uid_ref {
  int src;
  int trg;
  int count;
  int location;
};

struct joined_res {
  vector<Score> scores;
  vec_d permuted_scores;
};

class JoinMethod {
public:
  virtual joined_res join(uid_ref& uid) const = 0;
};

// class JoinMethod1Native : public JoinMethod {
// public:
//   virtual joined_res join(join_config& conf, vector<uid_ref>& uids, vector<int>& join_gene_signs, vec2d_d& value_table, vec2d_u16& permute_cases, paths_type* p_paths0, paths_type* p_paths1, paths_type* p_paths_res, uint64_t total_paths) const;
// };

// class JoinMethod2Native : public JoinMethod {
// public:
//   virtual joined_res join(join_config& conf, vector<uid_ref>& uids, vector<int>& join_gene_signs, vec2d_d& value_table, vec2d_u16& permute_cases, paths_type* p_paths0, paths_type* p_paths1, paths_type* p_paths_res, uint64_t total_paths) const;
// };

class PathSet {

public:

  const int size;
  const int width_ul;
  const int vlen;

  PathSet(const int size, const int width_ul, const int vlen) : size(size), width_ul(width_ul), vlen(vlen) {};
  virtual ~PathSet() {};

  virtual const uint64_t* operator[](int idx) const = 0;
  virtual void load(const vec2d_i& data) = 0;

  // create new path set from provided indices
  virtual unique_ptr<PathSet> select(const vector<int>& indices) const = 0;

protected:

  uint64_t* block = nullptr;

};

class PathSet_BlockM1 : public PathSet {
public:
  PathSet_BlockM1(const int size, const int width_ul);
  virtual ~PathSet_BlockM1();
  virtual const uint64_t* operator[](int idx) const;
  virtual void load(const vec2d_i& data);
  virtual unique_ptr<PathSet> select(const vector<int>& indices) const;
};

class PathSet_BlockM2 : public PathSet {
public:
  PathSet_BlockM2(const int size, const int width_ul);
  virtual ~PathSet_BlockM2();
  virtual const uint64_t* operator[](int idx) const;
  virtual void load(const vec2d_i& data);
  virtual unique_ptr<PathSet> select(const vector<int>& indices) const;
};

class JoinExec {

public:

  const int num_cases;
  const int num_ctrls;
  const int width_ul;
  int top_k = 12;
  int iterations = 0;
  int nthreads = 0;

  // TODO shouldn't really be here
  static inline unsigned count_total_paths(const vector<uid_ref>& uids) {
    unsigned total = 0;
    for(const auto& uid : uids)
      total += uid.count;
    return total;
  }

  JoinExec(const int num_cases, const int num_ctrls);

  ~JoinExec(){
    if(case_mask)
      delete[] case_mask;
    if(perm_case_mask)
      delete[] perm_case_mask;
  }

  void setValueTable(vec2d_d table);

  void setPermutedCases(vec2d_i& perm_cases);

  unique_ptr<PathSet> createPathSet(int size) const;
  
  joined_res join(const vector<uid_ref>& uids, const vector<int>& join_gene_signs, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const;

protected:

  // pad to uint64_t
  static inline int vector_width_ul(int num_cases, int num_ctrls) {
    return (int) ceil((num_cases + num_ctrls) / 64.0);
  }

  vec2d_d value_table;
  uint64_t* case_mask = nullptr;
  uint64_t* perm_case_mask = nullptr;


};

#endif
