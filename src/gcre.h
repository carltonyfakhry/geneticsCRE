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
  int size;
  int width_ul;
  int num_cases;
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
    short* index;
    uint64_t* pos;
    uint64_t* neg;
    uint64_t* con;
  };
};

struct joined_res {
  vec2d_u64 ids;
  vec_d scores;
  vec_d permuted_scores;
  paths_type* paths = NULL;
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

vec2d_u64 parseCaseORControl(Rcpp::IntegerMatrix CaseORControl, int nCases, int nControls);
void parsePaths(Rcpp::IntegerMatrix data, int nCases, int nControls, string file_path);
void StorePaths(vec2d_u64 &paths, string file_path);
vec2d_u64 readPaths(string file_path);
int getTotalPaths(Rcpp::IntegerVector trguids, Rcpp::List uids_CountLoc);
vec2d_u64 getZeroMatrix(int dim1, int dim2);
int getTotalCountsCountLoc(Rcpp::List uids_CountLoc);

joined_res* join_method2(vector<int> src_uids, vector<int> trg_uids, Rcpp::List uids_CountLoc, vector<int> join_gene_signs,
  vec2d_d value_table, int nCases, int nControls, int K,
  int iterations, Rcpp::IntegerMatrix CaseORControl, int pathLength, int nthreads,
  bool keep_joined,
  string pos_path1, string neg_path1, string conflict_path1,
  string pos_path2, string neg_path2, string conflict_path2,
  string dest_path_pos, string dest_path_neg, string dest_path_conflict, int total_paths);

Rcpp::List JoinIndicesMethod1(Rcpp::IntegerVector srcuid, Rcpp::IntegerVector trguids2, Rcpp::List uids_CountLoc, Rcpp::IntegerVector joining_gene_sign,
  Rcpp::NumericMatrix ValueTable, int nCases, int nControls, int K,
  int iterations, Rcpp::IntegerMatrix CaseORControl, int pathLength, int nthreads, string pos_path1,
  string pos_path2, string dest_path_pos);

Rcpp::List JoinIndicesMethod2(Rcpp::IntegerVector srcuid, Rcpp::IntegerVector trguids2, Rcpp::List uids_CountLoc, Rcpp::IntegerVector joining_gene_sign,
  Rcpp::NumericMatrix ValueTable, int nCases, int nControls, int K,
  int iterations, Rcpp::IntegerMatrix CaseORControl, int pathLength, int nthreads, string pos_path1,
  string neg_path1, string conflict_path1, string pos_path2, string neg_path2, string conflict_path2,
  string dest_path_pos, string dest_path_neg, string dest_path_conflict);

#endif
