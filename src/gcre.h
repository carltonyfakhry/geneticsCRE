#ifndef GCRE_H
#define GCRE_H

#include <inttypes.h>
#include <cmath>
#include <limits>
#include <vector>
#include <queue>
#include <iostream>
#include <fstream>

// these are types visible to R; keeping that minimal
#include "geneticsCRE_types.h"

using namespace std;

const uint64_t ONE_UL = 1;

typedef std::vector<uint64_t> vec_u64;
typedef std::vector<std::vector<uint64_t>> vec2d_u64;
typedef std::vector<double> vec_d;
typedef std::vector<std::vector<double>> vec2d_d;
typedef std::vector<std::vector<int>> vec2d_i;

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

struct paths_base : paths_type {
  int size = 0;
  int width_ul = 0;
  int num_cases = 0;
};

struct paths_vec : paths_base {
  vec2d_u64 pos;
  vec2d_u64 neg;
  vec2d_u64 con;
  void resize(int size, int width_ul, int num_cases);
};

struct paths_block : paths_base {
  struct path {
    uint64_t hash = 0;
    short lengh = 0;
    short* index = NULL;
    uint64_t* pos;
    uint64_t* neg;
    uint64_t* con;
  };
};

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

struct joined_res {
  vector<Score> scores;
  vec_d permuted_scores;
};

// TODO must be a better way to define overrides
class JoinMethod {
public:
  virtual paths_type* createPathSet() const = 0;
  virtual paths_type* createPathSet(vec2d_i& data, int num_cases, int num_controls) const = 0;
  virtual joined_res join(join_config& conf, vector<uid_ref>& uids, vector<int>& join_gene_signs, vec2d_d& value_table, vec2d_i& cases, paths_vec* paths0, paths_vec* paths1, paths_vec* paths_res, uint64_t total_paths) const = 0;
};

class JoinMethod2Vector : public JoinMethod {
public:
  virtual paths_type* createPathSet() const;
  virtual paths_type* createPathSet(vec2d_i& data, int num_cases, int num_controls) const;
  virtual joined_res join(join_config& conf, vector<uid_ref>& uids, vector<int>& join_gene_signs, vec2d_d& value_table, vec2d_i& cases, paths_vec* paths0, paths_vec* paths1, paths_vec* paths_res, uint64_t total_paths) const;
};

class JoinMethod2Native : public JoinMethod {
};

#endif
