#ifndef GCRE_H
#define GCRE_H

#include <inttypes.h>

// TODO isolate to join indices
#include <Rcpp.h>

// these are types visible to R; keeping that minimal
#include "geneticsCRE_types.h"

using namespace std;

typedef std::vector<uint64_t> vec_u64;
typedef std::vector<std::vector<uint64_t>> vec2d_u64;
typedef std::vector<double> vec_d;
typedef std::vector<std::vector<double>> vec2d_d;

struct score_t {
  int src_uid;
  int trg_uid;
  double score;
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

int getTotalPaths(Rcpp::IntegerVector trguids, Rcpp::List uids_CountLoc);
int getTotalCountsCountLoc(Rcpp::List uids_CountLoc);

joined_res* join_method2(vector<int>& src_uids, vector<int>& trg_uids, Rcpp::List& uids_CountLoc, vector<int>& join_gene_signs,
  vec2d_d& value_table, int nCases, int nControls, int K,
  int iterations, vec2d_u64& case_mask, int pathLength, int nthreads,
  paths_vec* paths0, paths_vec* paths1, paths_vec* paths_res, int total_paths);

#endif
