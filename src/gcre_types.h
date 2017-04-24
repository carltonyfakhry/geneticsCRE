#ifndef GCRE_TYPES_H
#define GCRE_TYPES_H

// low-level types needed in parts other than scoring

#include <cstdint>
#include <limits>
#include <vector>

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

// types that can hold max allowed values for various entities
using st_path_count = uint64_t;
using st_uids_size = uint32_t;
using st_pathset_size = uint32_t;

enum class Method { method1 = 1, method2 = 2 };

class Score {
public:
  double score = -std::numeric_limits<double>::infinity();
  int src = -1;
  int trg = -1;
  inline Score() {}
  inline Score(double score, int src, int trg) : score(score), src(src), trg(trg) {}
  // reverse sort for priority queue
  friend bool operator<(Score a, Score b) { return a.score > b.score; }
};

struct joined_res {
  std::vector<Score> scores;
  vec_d permuted_scores;
};

struct uid_ref {
  int src;
  int trg;
  int count;
  st_pathset_size location;
  st_path_count path_idx;
};

#endif