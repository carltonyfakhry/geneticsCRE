#include "gcre.h"

using namespace std;

PathSet_BlockM1::PathSet_BlockM1(const int size, const int width_ul) : PathSet(size, width_ul) {
  block = new uint64_t[size * width_ul] {};
  printf("[%p] allocated block for m1 paths: %d x %d (%lu bytes)\n", this, size, width_ul, size * width_ul * sizeof(uint64_t));
}

PathSet_BlockM1::~PathSet_BlockM1(){
  printf("[%p] m1 block cleanup\n", this);
  delete[] block;
}

const uint64_t* PathSet_BlockM1::operator[](int idx) const {
  return block + idx * width_ul;
}

void PathSet_BlockM1::load(vec_i data) {}

PathSet_BlockM2::PathSet_BlockM2(const int size, const int width_ul) : PathSet(size, width_ul) {
  block = new uint64_t[size * width_ul * 2] {};
  printf("[%p:%p] allocated block for m2 paths: %d x %d (%lu bytes)\n", this, block, size, width_ul, size * width_ul * 2 * sizeof(uint64_t));
}

PathSet_BlockM2::~PathSet_BlockM2(){
  printf("[%p:%p] m2 block cleanup\n", this, block);
  delete[] block;
}

const uint64_t* PathSet_BlockM2::operator[](int idx) const {
  return block + idx * width_ul * 2;
}

void PathSet_BlockM2::load(vec_i data) {}

unique_ptr<PathSet> JoinExec::createPathSet(int size) const {
  // std::make_unique is C++14
  return unique_ptr<PathSet>(new PathSet_BlockM2(size, num_cases + num_ctrls));
}

joined_res JoinExec::join(int path_length, const vector<uid_ref>& uids, const vector<int>& join_gene_signs, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const {
  return joined_res();
}
