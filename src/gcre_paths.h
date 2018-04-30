  #ifndef GCRE_PATHS_H
#define GCRE_PATHS_H

using namespace std;

// TODO consolidate class forward declares
class PathSet;
using TPathSet = unique_ptr<PathSet>;

class PathSet {

public:

  const st_pathset_size size;
  const uint16_t width_ul;
  const uint16_t vlen;

  // TODO abstract interface for other storage types
  PathSet(st_pathset_size size, const int width_ul, const uint16_t vlen) : size(size), width_ul(width_ul), vlen(vlen) {
    const size_t min_block_size = 16 * 1024 * 1024;
    const size_t block_size = size * vlen * sizeof(uint64_t);
    // if(block_size > min_block_size)
      // printf("allocate block for path-set (%5d x %d ul) %'15lu bytes: ", size, vlen, block_size);
    // block = unique_ptr<uint64_t[]>((uint64_t*) aligned_alloc(gs_align_size, sizeof(uint64_t) * size * vlen));
    try{
      block = unique_ptr<uint64_t[]>(new uint64_t[size * vlen]);
    }catch(std::bad_alloc &mess){
      printf("bad_alloc caught when dynammically allocating array: most likely cause is not enough memory available!\n");
    }

    // Touch allocated array
    try{
      uint64_t touch = block[size*vlen-1];
    }catch(std::exception &e){
      printf("Segmentation fault caught: most likely cause is not enough memory available to allocate array!\n");
    }

    for(int k = 0; k < size * vlen; k++)
      block[k] = bit_zero_ul;
    // if(block_size > min_block_size)
      // printf("[%p]\n", block.get());
  }

  inline const uint64_t* operator[](st_pathset_size idx) const {
    check_index(idx, size);
    return block.get() + idx * vlen;
  }

  inline void set(st_pathset_size idx, const uint64_t* data) {
    check_index(idx, size);
    memcpy(block.get() + idx * vlen, data, vlen * sizeof(uint64_t));
  }

  // TODO check input size dimensions
  // positives are loaded into the first half of the record, so m1/m2 loading is the same
  void load(const vec2d_i& data) {

    // printf("loading data to path-set (%lu x %lu): ", data.size(), data.size() > 0 ? data.front().size() : 0);

    check_true(size == data.size());
    long count_in = 0;
    for(auto r = 0; r < data.size(); r++) {
      check_index(data[r].size(), width_ul * 64);
      for(auto c = 0; c < data[r].size(); c++) {
        if(data[r][c] != 0) {
          count_in += 1;
          block[r * vlen + c/64] |= bit_one_ul << c % 64;
        }
      }
    }

    long count_set = 0;
    for(auto k = 0; k < size * vlen; k++)
      count_set += __builtin_popcountl(block[k]);

    // printf("%ld (of %ld)\n", count_set, count_in);
    check_true(count_set == count_in);
  }

  // create new path set from provided indices
  // this is only done for the initial path sets, and the index values come from R, so int instead of st_pathset_size
  TPathSet select(const vector<int>& indices) const {
    TPathSet pset(new PathSet(indices.size(), width_ul, vlen));
    for(int k = 0; k < pset->size; k++) {
      check_index(indices[k], size);
      memcpy(pset->block.get() + vlen * k, block.get() + vlen * indices[k], vlen * sizeof(uint64_t));
    }
    int count_set = 0;
    for(int k = 0; k < pset->size * vlen; k++)
      count_set += __builtin_popcountl(pset->block[k]);
    return pset;
  }

protected:

  unique_ptr<uint64_t[]> block;

};

#endif
