#ifndef GCRE_H
#define GCRE_H

// TODO isolate to join indices
#include <Rcpp.h>

using namespace std;

typedef std::vector<std::vector<double>> vecd2d; 

vector<vector<uint64_t>> parseCaseORControl(Rcpp::IntegerMatrix CaseORControl, int nCases, int nControls);
void parsePaths(Rcpp::IntegerMatrix data, int nCases, int nControls, string file_path);
void StorePaths(vector<vector<uint64_t>> &paths, string file_path);
vector<vector<uint64_t>> readPaths(string file_path);
int getTotalPaths(Rcpp::IntegerVector trguids, Rcpp::List uids_CountLoc);
vector<vector<uint64_t>> getZeroMatrix(int dim1, int dim2);
int getTotalCountsCountLoc(Rcpp::List uids_CountLoc);

Rcpp::List join_method2(Rcpp::IntegerVector srcuid, Rcpp::IntegerVector trguids2, Rcpp::List uids_CountLoc, Rcpp::IntegerVector joining_gene_sign,
  vecd2d value_table, int nCases, int nControls, int K,
  int iterations, Rcpp::IntegerMatrix CaseORControl, int pathLength, int nthreads, string pos_path1,
  string neg_path1, string conflict_path1, string pos_path2, string neg_path2, string conflict_path2,
  string dest_path_pos, string dest_path_neg, string dest_path_conflict);

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
