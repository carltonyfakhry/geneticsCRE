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

struct uid_ref {
  int src;
  int trg;
  int count;
  int location;
  int path_idx;
  vector<bool> signs;
};

struct joined_res {
  vector<Score> scores;
  vec_d permuted_scores;
};

class JoinMethod {
public:
  virtual joined_res join(uid_ref& uid) const = 0;
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

  void setPermutedCases(const vec2d_i& perm_cases);

  unique_ptr<PathSet> createPathSet(int size) const;
  
  joined_res join(const vector<uid_ref>& uids, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const;

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
