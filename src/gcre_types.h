#ifndef GCRE_TYPES_H
#define GCRE_TYPES_H

// low-level types needed in parts other than scoring

#include <cstdint>
#include <limits>
#include <vector>

const uint64_t bit_zero_ul = 0;
const uint64_t bit_one_ul = 1;

// TODO apparently aliases are now the thing?
using vec_i   = std::vector<int>;
using vec_d   = std::vector<double>;
using vec_u64 = std::vector<uint64_t>;

using vec2d_d   = std::vector<std::vector<double>>;
using vec2d_i   = std::vector<std::vector<int>>;
using vec2d_i8  = std::vector<std::vector<int8_t>>;
using vec2d_u16 = std::vector<std::vector<uint16_t>>;
using vec2d_u64 = std::vector<std::vector<uint64_t>>;

// types that can hold max allowed values for various entities
using st_path_count   = uint64_t;
using st_uids_size    = uint32_t;
using st_pathset_size = uint32_t;

enum class Method { method1 = 1, method2 = 2 };

class Score {
public:
  double score = -std::numeric_limits<double>::infinity();
  int src = -1;
  int trg = -1;
  int cases = 0;
  int ctrls = 0;
  inline Score() {}
  inline Score(double score, int src, int trg, int cases, int ctrls) : score(score), src(src), trg(trg), cases(cases), ctrls(ctrls) {}
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

#endif