#include <Rcpp.h>
#include "gcre.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// conversions from R

vector<int> copy_rvector(IntegerVector& r_vec) {
  vector<int> vec(r_vec.size());
  for(int k = 0; k < r_vec.size(); k += 1)
    vec[k] = r_vec[k];
  return vec;
}

vec2d_d copy_rmatrix(NumericMatrix& matrix) {
  vec2d_d table(matrix.nrow(), vector<double>(matrix.ncol()));
  for(int r = 0; r < matrix.nrow(); r += 1){
    for(int c = 0; c < matrix.ncol(); c += 1)
      table[r][c] = matrix(r,c);
  }
  printf("copied value matrix: %d x %d\n", matrix.nrow(), matrix.ncol());
  return table;
}

vec2d_i copy_rmatrix(IntegerMatrix& matrix) {
  vec2d_i table(matrix.nrow(), vector<int>(matrix.ncol()));
  for(int r = 0; r < matrix.nrow(); r += 1){
    for(int c = 0; c < matrix.ncol(); c += 1)
      table[r][c] = matrix(r,c);
  }
  printf("copied value matrix: %d x %d\n", matrix.nrow(), matrix.ncol());
  return table;
}

JoinMethod& create_method(string name) {
  if(name == "method2") {
    static JoinMethod2Vector method;
    return method;
  } else {
    stop("unknown method: " + name);
  }
}

// API exposed to R

// [[Rcpp::export]]
XPtr<paths_type> createPathSet(IntegerMatrix& r_data, int num_cases, int num_controls, std::string method_name){

  if(r_data.ncol() != num_cases + num_cases)
    stop("data width does not match case/control counts: " + to_string(r_data.ncol()));

  JoinMethod& method = create_method(method_name);
  vec2d_i data = copy_rmatrix(r_data);
  return XPtr<paths_type>(method.createPathSet(data, num_cases, num_controls), true) ;
}

// [[Rcpp::export]]
XPtr<paths_type> createEmptyPathSet(std::string method){
  return XPtr<paths_type>(create_method(method).createPathSet(), true) ;
}

/**
 *
 * Create list indexed by the names of gene uids. Each entry in the list
 * is a vector with the first entry being the count of paths starting
 * with the corresponding gene uids and the second entry be the location
 * of the first path in the relations data frame.
 *
 */
// [[Rcpp::export]]
List getMatchingList(IntegerVector& uids, IntegerVector& counts, IntegerVector& location){
  List uids_2countsloc = List();
  for(unsigned int i = 0; i < uids.size(); i++){
    string uids_str = to_string(uids[i]);
    IntegerVector temp = IntegerVector(2);
    int count = counts[i];
    temp[0] = count;
    temp[1] = location[i];
    uids_2countsloc[uids_str] = temp;
  }
  return uids_2countsloc;
}

// [[Rcpp::export]]
List JoinIndices(IntegerVector r_src_uids, IntegerVector r_trg_uids, List r_uid_count_locs, IntegerVector r_join_gene_signs,
  NumericMatrix r_value_table, int num_cases, int num_controls, int top_k,
  int iterations, IntegerMatrix r_cases, std::string method_name, int path_length, int nthreads,
  SEXP xp_paths0, SEXP xp_paths1, SEXP xp_paths_res){

  printf("\n################################\n");
  printf("  method      : %s\n", method_name.c_str());
  printf("  path_length : %d\n", path_length);
  printf("  uid size    : %d\n", r_trg_uids.size());
  printf("################################\n\n");

  if(r_trg_uids.size() != r_src_uids.size())
    stop("uid size mismatch");

  // convert some common datatypes
  vec2d_i cases = copy_rmatrix(r_cases);
  vec2d_d value_table = copy_rmatrix(r_value_table);
  vector<int> join_gene_signs = copy_rvector(r_join_gene_signs);

  // join parameters
  join_config conf;
  conf.num_cases = num_cases;
  conf.num_controls = num_controls;
  conf.top_k = top_k;
  conf.path_length = path_length;
  conf.iterations = iterations;
  conf.nthreads = nthreads;

  // combobulate uid and location/count data
  uint64_t total_paths = 0;
  vector<uid_ref> uids(r_trg_uids.size(), uid_ref());
  for(int k = 0; k < r_trg_uids.size(); k++){
    uid_ref& uid = uids[k];
    uid.trg = r_trg_uids[k];
    uid.src = r_src_uids[k];
    IntegerVector count_loc = r_uid_count_locs[to_string(uid.trg)];
    uid.count = count_loc[0];
    uid.location = count_loc[1];
    total_paths += uid.count;
  }

  // pointers to method-specific bitset data
  paths_vec* paths1 = (paths_vec*) XPtr<paths_type>(xp_paths1).get();

  paths_vec* paths0 = NULL;
  if(!Rf_isNull(xp_paths0))
    paths0 = (paths_vec*) XPtr<paths_type>(xp_paths0).get();

  paths_vec* paths_res = NULL;
  if(!Rf_isNull(xp_paths_res))
    paths_res = (paths_vec*) XPtr<paths_type>(xp_paths_res).get();

  // get appropriate method type
  JoinMethod& method = create_method(method_name);

  joined_res res = method.join(conf, uids, join_gene_signs, value_table, cases, paths0, paths1, paths_res, total_paths);

  NumericVector permuted_scores(res.permuted_scores.size());
  for(int k = 0; k < res.permuted_scores.size(); k++)
    permuted_scores[k] = res.permuted_scores[k];

  NumericVector scores(res.scores.size());
  IntegerMatrix ids(res.scores.size(), 2);
  for(int k = 0; k < res.scores.size(); k++){
    Score& score = res.scores[k];
    scores[k] = score.score;
    // add 1 because it will be used as an index in R
    ids(k,0) = score.src + 1;
    ids(k,1) = score.trg + 1;
  }

  return List::create(Named("scores") = scores, Named("ids") = ids, Named("TestScores") = permuted_scores);
}