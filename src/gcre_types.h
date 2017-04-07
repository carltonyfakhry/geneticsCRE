#ifndef GCRE_TYPES_H
#define GCRE_TYPES_H

#include <cstdint>
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

struct uid_ref {
  int src;
  int trg;
  int count;
  int location;
  int path_idx;
};

#endif