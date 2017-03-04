#include <fstream>
#include "gcre.h"
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

const uint64_t ONE_UL = 1;

// need consistent padding for different alignments
static int vector_width_ul(int count){
  return (int) ceil(count / 64.0);
}

/**
 *
 * Get the total number of paths to be processed.
 *
 */
int getTotalPaths(IntegerVector trguids, List uids_CountLoc){
  int total_paths = 0;
  for(int i = 0; i < trguids.size(); i++){
    int uid = trguids[i];
    string geneuid = to_string(uid);
    IntegerVector temp = as<IntegerVector>(uids_CountLoc[geneuid]);
    total_paths += temp[0];
  }
  return total_paths;
}

/**
 *
 * Get the total counts in uids_CountLoc.
 *
 */
int getTotalCountsCountLoc(List uids_CountLoc){
  int total_paths = 0;
  for(List::iterator it = uids_CountLoc.begin(); it != uids_CountLoc.end(); ++it){
    IntegerVector temp = as<IntegerVector>(*it);
    total_paths += temp[0];
  }
  return total_paths;
}



// conversions from R

vec2d_u64& create_case_mask(IntegerMatrix& r_cases, int num_cases, int num_controls){

  int vlenc = vector_width_ul(num_cases);
  int vlent = vector_width_ul(num_controls);

  vec2d_u64& mask = *new vec2d_u64(r_cases.nrow(), vec_u64(vlenc + vlent, 0));
  for(int i = 0; i < r_cases.nrow(); i++){
    for(int j = 0; j < num_cases; j++){
      int index = j / 64;
      if(r_cases(i,j) == 1)
        mask[i][index] |= ONE_UL << j % 64;
    }
    for(int j = num_cases; j < num_cases + num_controls; j++){
      int index = vlenc + (j-num_cases) / 64;
      if(r_cases(i,j) == 1)
        mask[i][index] |= ONE_UL << (j-num_cases) % 64;
    }

  }

  return mask;
}

vector<int>& copy_rvector(IntegerVector& r_vec)
{
  vector<int>& vec = *new vector<int>(r_vec.size());
  for(int k = 0; k < r_vec.size(); k += 1)
    vec[k] = r_vec[k];
  return vec;
}

vec2d_d& copy_rmatrix(NumericMatrix& matrix)
{
  vec2d_d& table = *new vec2d_d(matrix.nrow(), vector<double>(matrix.ncol()));
  for(int r = 0; r < matrix.nrow(); r += 1){
    for(int c = 0; c < matrix.ncol(); c += 1)
      table[r][c] = matrix(r,c);
  }
  printf("copied value matrix: %d x %d\n", matrix.nrow(), matrix.ncol());
  return table;
}

/**
 *
 * This function parses paths into their 64 bit representations.
 *
 */
// [[Rcpp::export]]
XPtr<paths_type> createPathSet(IntegerMatrix& data, int num_cases, int num_controls, std::string method){

  if(data.ncol() != num_cases + num_cases)
    stop("data width does not match case/control counts: " + to_string(data.ncol()));

  int vlenc = vector_width_ul(num_cases);
  int vlent = vector_width_ul(num_controls);

  // allocate on heap
  paths_vec* paths = new paths_vec;

  paths->size = data.nrow();
  paths->width_ul = vlenc + vlent;
  paths->num_cases = num_cases;
  paths->pos.resize(paths->size, vec_u64(paths->width_ul, 0));
  paths->neg.resize(paths->size, vec_u64(paths->width_ul, 0));
  paths->con.resize(paths->size, vec_u64(paths->width_ul, 0));

  for(int r = 0; r < data.nrow(); r++){
    for(int c = 0; c < data.ncol(); c++){
      if(data(r,c) != 0){
        if(c < num_cases)
          paths->pos[r][c/64] |= ONE_UL << c % 64;
        else
          paths->pos[r][vlenc + (c-num_cases)/64] |= ONE_UL << (c-num_cases) % 64;
      }
    }
  }

  return XPtr<paths_type>(paths, true) ;
}

// [[Rcpp::export]]
XPtr<paths_type> createEmptyPathSet(std::string method){
  paths_vec* paths = new paths_vec;
  paths->size = 0;
  paths->width_ul = 0;
  return XPtr<paths_type>(paths, true) ;
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
    Rcpp::IntegerVector temp = IntegerVector(2);
    int count = counts[i];
    temp[0] = count;
    temp[1] = location[i];
    uids_2countsloc[uids_str] = temp;
  }
  return uids_2countsloc;
}

// [[Rcpp::export]]
List JoinIndices(IntegerVector r_src_uids, IntegerVector r_trg_uids, List uids_CountLoc, IntegerVector r_join_gene_signs,
  NumericMatrix r_value_table, int num_cases, int num_controls, int K,
  int iterations, IntegerMatrix r_case_mask, std::string method, int pathLength, int nthreads,
  SEXP xp_paths0, SEXP xp_paths1, SEXP xp_paths_res){

  vec2d_u64& case_mask = create_case_mask(r_case_mask, num_cases, num_controls);
  vec2d_d& value_table = copy_rmatrix(r_value_table);
  vector<int>& src_uids = copy_rvector(r_src_uids);
  vector<int>& trg_uids = copy_rvector(r_trg_uids);
  vector<int>& join_gene_signs = copy_rvector(r_join_gene_signs);

  paths_vec* paths1 = (paths_vec*) XPtr<paths_type>(xp_paths1).get();

  // TODO fix leak
  paths_vec* paths0 = NULL;
  if(!Rf_isNull(xp_paths0)){
    paths0 = (paths_vec*) XPtr<paths_type>(xp_paths0).get();
  }else{
    paths0 = new paths_vec;
    paths0->size = trg_uids.size();
    paths0->width_ul = 50;
    paths0->num_cases = num_cases;
    paths0->pos.resize(paths0->size, vec_u64(paths0->width_ul, 0));
    paths0->neg.resize(paths0->size, vec_u64(paths0->width_ul, 0));
    paths0->con.resize(paths0->size, vec_u64(paths0->width_ul, 0));
    printf("  ** resized zero matrix: %d x %d\n", paths0->size, paths0->width_ul);
  }

  paths_vec* paths_res = NULL;
  if(!Rf_isNull(xp_paths_res))
    paths_res = (paths_vec*) XPtr<paths_type>(xp_paths_res).get();

  printf("################\n");
  // printf("srcuid size : %d\n", srcuid.size());
  // printf("trguid size : %d\n", trguids2.size());
  // printf("uids cl size: %d\n", uids_CountLoc.size());
  // printf("join genes  : %d\n", joining_gene_sign.size());
  // printf("value table : %d x %d\n", ValueTable.rows(), ValueTable.cols());
  // printf("caseorcon   : %d x %d\n", CaseORControl.rows(), CaseORControl.cols());
  // printf("num_cases      : %d\n", num_cases);
  // printf("num_controls   : %d\n", num_controls);
  // printf("K           : %d\n", K);
  // printf("iterations  : %d\n", iterations);
  printf("method      : %s\n", method.c_str());
  printf("pathlen     : %d\n", pathLength);
  printf("nthread     : %d\n", nthreads);
  printf("\n");
  printf("################\n\n");

  if(method == "method2") {

    joined_res* res = join_method2_new(src_uids, trg_uids, uids_CountLoc, join_gene_signs,
      value_table, num_cases, num_controls, K,
      iterations, case_mask, pathLength, nthreads,
      paths0, paths1, paths_res,
      getTotalPaths(r_trg_uids, uids_CountLoc));

    NumericVector scores(res->scores.size());
    for(int k = 0; k < res->scores.size(); k++)
      scores[k] = res->scores[k];

    NumericVector permuted_scores(res->permuted_scores.size());
    for(int k = 0; k < res->permuted_scores.size(); k++)
      permuted_scores[k] = res->permuted_scores[k];

    IntegerMatrix ids(res->ids.size(),2);
    for(int k = 0; k < res->ids.size(); k++){
      ids(k,0) = res->ids[k][0];
      ids(k,1) = res->ids[k][1];
    }

    return List::create(Named("scores") = scores, Named("ids") = ids, Named("TestScores") = permuted_scores);

  } else {
    stop("unknown method: " + method);
  }
}