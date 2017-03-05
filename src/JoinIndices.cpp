#include <Rcpp.h>
#include "gcre.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

const uint64_t ONE_UL = 1;

// need consistent padding for different alignments
static int vector_width_ul(int count){
  return (int) ceil(count / 64.0);
}

// conversions from R

vec2d_u64 create_case_mask(IntegerMatrix& r_cases, int num_cases, int num_controls){

  int vlenc = vector_width_ul(num_cases);
  int vlent = vector_width_ul(num_controls);

  vec2d_u64 mask(r_cases.nrow(), vec_u64(vlenc + vlent, 0));
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

vector<int> copy_rvector(IntegerVector& r_vec)
{
  vector<int> vec(r_vec.size());
  for(int k = 0; k < r_vec.size(); k += 1)
    vec[k] = r_vec[k];
  return vec;
}

vec2d_d copy_rmatrix(NumericMatrix& matrix)
{
  vec2d_d table(matrix.nrow(), vector<double>(matrix.ncol()));
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
List JoinIndices(IntegerVector r_src_uids, IntegerVector r_trg_uids, List r_uid_count_locs, IntegerVector r_join_gene_signs,
  NumericMatrix r_value_table, int num_cases, int num_controls, int top_k,
  int iterations, IntegerMatrix r_case_mask, std::string method, int path_length, int nthreads,
  SEXP xp_paths0, SEXP xp_paths1, SEXP xp_paths_res){

  if(r_trg_uids.size() != r_src_uids.size())
    stop("uid size mismatch");

  vec2d_u64 case_mask = create_case_mask(r_case_mask, num_cases, num_controls);
  vec2d_d value_table = copy_rmatrix(r_value_table);
  vector<int> join_gene_signs = copy_rvector(r_join_gene_signs);

  join_config conf;
  conf.num_cases = num_cases;
  conf.num_controls = num_controls;
  conf.top_k = top_k;
  conf.path_length = path_length;
  conf.iterations = iterations;
  conf.nthreads = nthreads;


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
  // for(int k = 0; k < 32; k++)
  //   printf(" %d:%d(%d..%d)", uids[k].src, uids[k].trg, uids[k].location, uids[k].count);
  // printf("\n");

  paths_vec* paths1 = (paths_vec*) XPtr<paths_type>(xp_paths1).get();

  // TODO fix leak
  paths_vec* paths0 = NULL;
  if(!Rf_isNull(xp_paths0)){
    paths0 = (paths_vec*) XPtr<paths_type>(xp_paths0).get();
  }else{
    paths0 = new paths_vec;
    paths0->size = r_trg_uids.size();
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
  // printf("srcuid size : %d\n", r_src_uids.size());
  // printf("trguid size : %d\n", r_trg_uids.size());
  // printf("uids cl size: %d\n", r_uid_count_locs.size());
  // printf("join genes  : %d\n", r_join_gene_signs.size());
  // printf("value table : %d x %d\n", ValueTable.rows(), ValueTable.cols());
  // printf("caseorcon   : %d x %d\n", CaseORControl.rows(), CaseORControl.cols());
  // printf("num_cases      : %d\n", num_cases);
  // printf("num_controls   : %d\n", num_controls);
  // printf("K           : %d\n", K);
  // printf("iterations  : %d\n", iterations);
  printf("method      : %s\n", method.c_str());
  printf("pathlen     : %d\n", path_length);
  printf("nthread     : %d\n", nthreads);
  printf("\n");
  printf("################\n\n");

  if(method == "method2") {

    joined_res* res = join_method2_new(conf, uids, join_gene_signs, value_table, case_mask, paths0, paths1, paths_res, total_paths);

    NumericVector permuted_scores(res->permuted_scores.size());
    for(int k = 0; k < res->permuted_scores.size(); k++)
      permuted_scores[k] = res->permuted_scores[k];

    NumericVector scores(res->scores.size());
    IntegerMatrix ids(res->scores.size(), 2);
    for(int k = 0; k < res->scores.size(); k++){
      Score& score = res->scores[k];
      scores[k] = score.score;
      // add 1 because it will be used as an index in R
      ids(k,0) = score.src + 1;
      ids(k,1) = score.trg + 1;
    }

    return List::create(Named("scores") = scores, Named("ids") = ids, Named("TestScores") = permuted_scores);

  } else {
    stop("unknown method: " + method);
  }
}