#include "gcre.h"

using namespace std;

PathSet_BlockM1::PathSet_BlockM1(const int size, const int width_ul) : PathSet(size, width_ul, width_ul) {
  block = new uint64_t[size * vlen] {};
  printf("[%p] allocated block for m1 paths: %d x %d (%lu bytes)\n", this, size, width_ul, size * vlen * sizeof(uint64_t));
}

PathSet_BlockM1::~PathSet_BlockM1(){
  printf("[%p] m1 block cleanup\n", this);
  delete[] block;
}

const uint64_t* PathSet_BlockM1::operator[](int idx) const {
  return block + idx * vlen;
}

void PathSet_BlockM1::load(const vec2d_i& data) {}

PathSet_BlockM2::PathSet_BlockM2(const int size, const int width_ul) : PathSet(size, width_ul, width_ul * 2) {
  block = new uint64_t[size * vlen] {};
  printf("[%p:%p] allocated block for m2 paths: %d x %d (%lu bytes)\n", this, block, size, width_ul, size * vlen * sizeof(uint64_t));
}

PathSet_BlockM2::~PathSet_BlockM2(){
  printf("[%p:%p] m2 block cleanup\n", this, block);
  delete[] block;
}

const uint64_t* PathSet_BlockM2::operator[](int idx) const {
  return block + idx * vlen;
}

// interleave pos/neg uint64 blocks
void PathSet_BlockM2::load(const vec2d_i& data) {

  printf("loading data to path set: %lu x %lu\n", data.size(), data.size() > 0 ? data.front().size() : 0);

  int count_in = 0;
  for(int r = 0; r < data.size(); r++) {
    for(int c = 0; c < data[r].size(); c++) {
      if(data[r][c] != 0) {
        block[r * vlen + c/64*2] |= bit_one_ul << c % 64;
        count_in += 1;
      }
    }
  }

  int count_set = 0;
  for(int k = 0; k < size * vlen; k++)
    count_set += __builtin_popcountl(block[k]);
  printf("total set: %d (input: %d)\n", count_set, count_in);
}


// std::vector<std::vector<uint64_t> > matchData(const std::vector<std::vector<uint64_t> > &parseddata, IntegerVector data_inds){
//   std::vector<std::vector<uint64_t> > matcheddata(data_inds.size(), std::vector<uint64_t>(parseddata[0].size(),0));
//   for(int i = 0; i < data_inds.size(); i++){
//     int data_index = data_inds[i];
//     matcheddata[i] = parseddata[data_index];
//   }
//   return matcheddata;
// }
// }



JoinExec::JoinExec(const int num_cases, const int num_ctrls) : num_cases(num_cases), num_ctrls(num_ctrls), width_ul(vector_width_ul(num_cases, num_ctrls)) {

  // create case mask from case/control ranges
  case_mask = new uint64_t[width_ul];
  for(int k = 0; k < width_ul; k++)
    case_mask[k] = bit_zero_ul;
  for(int k = 0; k < num_cases; k++)
    case_mask[k/64] |= bit_one_ul << k % 64;

  printf("init exec: %d x %d (width: %d ul)\n", num_cases, num_ctrls, width_ul);
}

// TODO size check and precompute max table
// create a copy of the vt for simplicity
void JoinExec::setValueTable(vec2d_d table){
  printf("setting value table: %lu x %lu\n", table.size(), table.size() > 0 ? table.front().size() : 0);
  value_table = table;
}


unique_ptr<PathSet> JoinExec::createPathSet(int size) const {
  // std::make_unique is C++14
  return unique_ptr<PathSet>(new PathSet_BlockM2(size, num_cases + num_ctrls));
}

joined_res JoinExec::join(int path_length, const vector<uid_ref>& uids, const vector<int>& join_gene_signs, const PathSet& paths0, const PathSet& paths1, PathSet& paths_res) const {
  return joined_res();
}
