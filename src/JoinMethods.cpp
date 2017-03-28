#include <fstream>

#include <inttypes.h>
#include <stdio.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <cmath>
#include <limits>
#include <vector>
#include <set>
#include <queue>
#include <cstdio>
#include <iostream>
#include <fstream>

#include <Rcpp.h>

#include "gcre.h"


// [[Rcpp::plugins(cpp11)]]

// using namespace Rcpp;
using namespace std;

/**
 *
 * Create Rcpp::List indexed by the names of gene uids. Each entry in the Rcpp::List
 * is a vector with the first entry being the count of paths starting
 * with the corresponding gene uids and the second entry be the location
 * of the first path in the relations data frame.
 *
 */
// TODO why not do this in R directly? --dmitri
// [[Rcpp::export]]
Rcpp::List getMatchingList(Rcpp::IntegerVector uids, Rcpp::IntegerVector counts, Rcpp::IntegerVector location){
  Rcpp::List count_locs;
  for(int i = 0; i < uids.size(); i++){
    Rcpp::IntegerVector count_loc(2);
    count_loc[0] = counts[i];
    count_loc[1] = location[i];
    count_locs[to_string(uids[i])] = count_loc;
  }
  return count_locs;
}

static vector<int> copy_r(const Rcpp::IntegerVector& r_vec) {
  vector<int> vec(r_vec.size());
  for(int k = 0; k < r_vec.size(); k += 1)
    vec[k] = r_vec[k];
  return vec;
}

static vec2d_d copy_r(const Rcpp::NumericMatrix& matrix) {
  vec2d_d table(matrix.nrow(), vector<double>(matrix.ncol()));
  for(int r = 0; r < matrix.nrow(); r += 1){
    for(int c = 0; c < matrix.ncol(); c += 1)
      table[r][c] = matrix(r,c);
  }
  printf("copied value matrix: %d x %d\n", matrix.nrow(), matrix.ncol());
  return table;
}

static vec2d_i copy_r(const Rcpp::IntegerMatrix& matrix) {
  vec2d_i table(matrix.nrow(), vector<int>(matrix.ncol()));
  for(int r = 0; r < matrix.nrow(); r += 1){
    for(int c = 0; c < matrix.ncol(); c += 1)
      table[r][c] = matrix(r,c);
  }
  printf("copied value matrix: %d x %d\n", matrix.nrow(), matrix.ncol());
  return table;
}


static vector<uid_ref> assemble_uids(int path_length, const Rcpp::IntegerVector& r_src_uids, const Rcpp::IntegerVector& r_trg_uids, const Rcpp::List& r_count_locs, const Rcpp::IntegerVector& r_signs){
  long total_paths = 0;
  vector<uid_ref> uids(r_trg_uids.size(), uid_ref());
  for(int k = 0; k < r_trg_uids.size(); k++){
    auto& uid = uids[k];
    uid.trg = r_trg_uids[k];
    uid.src = r_src_uids[k];
    const Rcpp::IntegerVector& count_loc = r_count_locs[to_string(uid.trg)];
    uid.count = count_loc[0];
    uid.location = count_loc[1];

    // set begin index of result path storage block
    uid.path_idx = total_paths;

    total_paths += uid.count;

    // store computed join signs
    // TODO check correctness (path_length >5?)
    for(int loc = uid.location; loc < uid.location + uid.count; loc++) {
      int sign = 0;
      if(path_length > 3)
        sign = r_signs[k];
      else if(path_length < 3)
        sign = r_signs[loc];
      else if(path_length == 3)
        sign = (r_signs[k] + r_signs[loc] == 0) ? -1 : 1;
      uid.signs.push_back(sign == 1);
    }

  }
  printf("path_length: %ld, uids: %lu (paths: %ld)\n", path_length, uids.size(), total_paths);
  return uids;
}

// [[Rcpp::export]]
Rcpp::List ProcessPaths(Rcpp::IntegerVector r_src_uids1, Rcpp::IntegerVector r_trg_uids1, Rcpp::List r_count_locs1, Rcpp::IntegerVector r_signs1,
  Rcpp::IntegerVector r_src_uids1_2, Rcpp::IntegerVector r_trg_uids1_2, Rcpp::List r_count_locs1_2, Rcpp::IntegerVector r_signs1_2,
  Rcpp::IntegerVector r_src_uids2, Rcpp::IntegerVector r_trg_uids2, Rcpp::List r_count_locs2, Rcpp::IntegerVector r_signs2,
  Rcpp::IntegerVector r_src_uids3, Rcpp::IntegerVector r_trg_uids3, Rcpp::List r_count_locs3, Rcpp::IntegerVector r_signs3,
  Rcpp::IntegerVector r_src_uids4, Rcpp::IntegerVector r_trg_uids4, Rcpp::List r_count_locs4, Rcpp::IntegerVector r_signs4,
  Rcpp::IntegerVector r_src_uids5, Rcpp::IntegerVector r_trg_uids5, Rcpp::List r_count_locs5, Rcpp::IntegerVector r_signs5,
  Rcpp::IntegerVector r_data_inds1, Rcpp::IntegerVector r_data_inds1_2, Rcpp::IntegerVector r_data_inds2, Rcpp::IntegerVector r_data_inds3,
  Rcpp::IntegerMatrix r_data1, Rcpp::IntegerMatrix r_data2, Rcpp::NumericMatrix r_value_table,
  int num_cases, int num_ctrls, int top_k, int iterations, Rcpp::IntegerMatrix r_perm_cases, int method, int path_length, int nthreads){

  printf("meh\n");

  auto uids1a = assemble_uids(1, r_src_uids1, r_trg_uids1, r_count_locs1, r_signs1);
  auto uids1b = assemble_uids(1, r_src_uids1_2, r_trg_uids1_2, r_count_locs1_2, r_signs1_2);
  auto uids2 = assemble_uids(2, r_src_uids2, r_trg_uids2, r_count_locs2, r_signs2);
  auto uids3 = assemble_uids(3, r_src_uids3, r_trg_uids3, r_count_locs3, r_signs3);
  auto uids4 = assemble_uids(4, r_src_uids4, r_trg_uids4, r_count_locs4, r_signs4);
  auto uids5 = assemble_uids(5, r_src_uids5, r_trg_uids5, r_count_locs5, r_signs5);

  auto data_inx1a = copy_r(r_data_inds1);
  auto data_inx1b = copy_r(r_data_inds1_2);
  auto data_inx2 = copy_r(r_data_inds3);
  auto data_inx3 = copy_r(r_data_inds3);

  auto data1 = copy_r(r_data1);
  auto data2 = copy_r(r_data2);

  auto value_table = copy_r(r_value_table);
  auto perm_cases = copy_r(r_perm_cases);

  printf("done\n");
  exit(0);

  return Rcpp::List::create();
}
