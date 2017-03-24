#ifndef GCRE_H
#define GCRE_H

#include <inttypes.h>
#include <stdio.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
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

const uint64_t ZERO_UL = 0;
const uint64_t ONE_UL = 1;

typedef std::vector<int> vec_i;
typedef std::vector<double> vec_d;
typedef std::vector<uint64_t> vec_u64;

typedef std::vector<std::vector<double>> vec2d_d;
typedef std::vector<std::vector<int>> vec2d_i;
typedef std::vector<std::vector<int8_t>> vec2d_i8;
typedef std::vector<std::vector<uint16_t>> vec2d_u16;
typedef std::vector<std::vector<uint64_t>> vec2d_u64;


struct uid_ref {
  int src;
  int trg;
  int count;
  int location;
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

#include "geneticsCRE_types.h"

class PathSet {

public:

  const int size;
  const int width_ul;

  PathSet(const int size, const int width_ul) : size(size), width_ul(width_ul) {};
  virtual ~PathSet() {};

  virtual const uint64_t* operator[](int idx) const = 0;
  virtual void load(vec_i data) = 0;

protected:

  uint64_t* block = nullptr;

};

class PathSet_BlockM1 : public PathSet {
public:
  PathSet_BlockM1(const int size, const int width_ul);
  virtual ~PathSet_BlockM1();
  virtual const uint64_t* operator[](int idx) const;
  virtual void load(vec_i data);
};

class PathSet_BlockM2 : public PathSet {
public:
  PathSet_BlockM2(const int size, const int width_ul);
  virtual ~PathSet_BlockM2();
  virtual const uint64_t* operator[](int idx) const;
  virtual void load(vec_i data);
};

class JoinExec {

public:

  const int num_cases;
  const int num_ctrls;
  vec2d_d value_table;
  vec2d_u16 permute_cases;
  int top_k = 12;
  int iterations = 0;
  int nthreads = 0;

  JoinExec(const int num_cases, const int num_ctrls) : num_cases(num_cases), num_ctrls(num_ctrls) {}
  ~JoinExec();

  PathSet& createPathSet(int size) const;
  joined_res join(int path_length, const vector<uid_ref>& uids, const vector<int>& join_gene_signs, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const;

private:

  mutable vector<PathSet*> paths;

};

#endif
