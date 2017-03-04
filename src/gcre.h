#ifndef GCRE_H
#define GCRE_H

#include <vector>
#include <inttypes.h>

// these are types visible to R; keeping that minimal
#include "geneticsCRE_types.h"

using namespace std;

typedef std::vector<uint64_t> vec_u64;
typedef std::vector<std::vector<uint64_t>> vec2d_u64;
typedef std::vector<double> vec_d;
typedef std::vector<std::vector<double>> vec2d_d;

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

struct score_t {
  int src_uid;
  int trg_uid;
  double score;
};

class JoinMethod {
public:
  virtual paths_type* createPathSet() const = 0;
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

struct joined_res {
  vec2d_u64 ids;
  vec_d scores;
  vec_d permuted_scores;
};

// // struct path_set {
// //   int len;
// //   path_set(int len);
// //   ~path_set();
// // };

// path_set::path_set(int len) : len(len) {
//   printf("i am created: %d, %d\n", len, len);
// }

// path_set::~path_set(){
//   printf("i am destructed\n");
// }

joined_res* join_method2_new(join_config& conf, vector<uid_ref>& uids,
  vector<int>& join_gene_signs, vec2d_d& value_table, vec2d_u64& case_mask,
  paths_vec* paths0, paths_vec* paths1, paths_vec* paths_res, uint64_t total_paths);

// joined_res* join_method2(vector<int>& src_uids, vector<int>& trg_uids, Rcpp::List& uids_CountLoc, vector<int>& join_gene_signs,
//   vec2d_d& value_table, int num_cases, int num_controls, int top_k,
//   int iterations, vec2d_u64& case_mask, int path_length, int nthreads,
//   paths_vec* paths0, paths_vec* paths1, paths_vec* paths_res, int total_paths);

#endif
