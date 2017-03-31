#include <unordered_map>
#include <array>
#include <Rcpp.h>
#include "gcre.h"

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


static UidRelSet assemble_uids(int path_length, const Rcpp::IntegerVector& r_src_uids, const Rcpp::IntegerVector& r_trg_uids, const Rcpp::List& r_count_locs, const Rcpp::IntegerVector& r_signs){

  printf("start uid assemble for length: %d (src: %lu, locs: %lu)\n", path_length, r_src_uids.size(), r_count_locs.size());

  auto start = chrono::system_clock::now();

  // Rcpp named Lists are ridiculously slow, and keeping that format would required a lot of int<->string conversions
  unordered_map<int, array<int,2>> count_locs;
  const Rcpp::CharacterVector& r_names = r_count_locs.names();
  for(const auto& r_str : r_names) {
    string uid = Rcpp::as<string>(r_str);
    const Rcpp::IntegerVector& count_loc = r_count_locs[uid];
    count_locs[stoi(uid)] = {count_loc[0], count_loc[1]};
  }

  long total_paths = 0;

  // TODO not entirely sure this allocates the way I think it does
  vector<uid_ref> uids(r_trg_uids.size(), uid_ref());

  for(int k = 0; k < r_trg_uids.size(); k++){

    auto& uid = uids[k];
    uid.trg = r_trg_uids[k];
    uid.src = r_src_uids[k];
    const auto& cloc = count_locs[uid.trg];
    uid.count = cloc[0];
    uid.location = cloc[1];

    // set begin index of result path storage block
    uid.path_idx = total_paths;
    total_paths += uid.count;

  }

  auto time = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - start);
  printf("path_length: %d, uids: %lu, paths: %ld (%'ld ms)\n", path_length, uids.size(), total_paths, time.count());

  return UidRelSet(path_length, uids, copy_r(r_signs));

}

static Rcpp::List make_score_list(const joined_res& res) {

  Rcpp::NumericVector permuted_scores(res.permuted_scores.size());

  for(int k = 0; k < res.permuted_scores.size(); k++)
    permuted_scores[k] = res.permuted_scores[k];

  Rcpp::NumericVector scores(res.scores.size());
  Rcpp::IntegerMatrix ids(res.scores.size(), 2);
  for(int k = 0; k < res.scores.size(); k++){
    const auto& score = res.scores[k];
    scores[k] = score.score;
    // add 1 because it will be used as an index in R
    ids(k,0) = score.src + 1;
    ids(k,1) = score.trg + 1;
  }

  return Rcpp::List::create(Rcpp::Named("scores") = scores, Rcpp::Named("ids") = ids, Rcpp::Named("TestScores") = permuted_scores);
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
  int num_cases, int num_ctrls, int top_k, int iterations, Rcpp::IntegerMatrix r_perm_cases, std::string method, int path_length, int nthreads){

  JoinExec exec(method, num_cases, num_ctrls, iterations);
  exec.top_k = top_k;
  exec.nthreads = nthreads;

  setlocale(LC_ALL, "");

  printf("\nstarting join:\n\n");
  printf("       method : %d\n", exec.num_ctrls);
  printf("  path_length : %d\n", path_length);
  printf("        cases : %d\n", exec.num_cases);
  printf("     controls : %d\n", exec.num_ctrls);
  printf("      threads : %d\n", nthreads);
  printf("  iterations : %d\n", iterations);

  Timer::print_header();

  auto data_idx1a = copy_r(r_data_inds1);
  auto data_idx1b = copy_r(r_data_inds1_2);
  auto data_idx2 = copy_r(r_data_inds3);
  auto data_idx3 = copy_r(r_data_inds3);

  auto data1 = copy_r(r_data1);
  auto data2 = copy_r(r_data2);

  exec.setValueTable(copy_r(r_value_table));
  exec.setPermutedCases(copy_r(r_perm_cases));

  auto parsed_data1 = exec.createPathSet(data1.size());
  parsed_data1->load(data1);

  auto zero_set = exec.createPathSet(0);

  unique_ptr<PathSet> paths1 = nullptr, paths2 = nullptr, paths3 = nullptr;

  auto results = Rcpp::List::create(Rcpp::Named("lst1"), Rcpp::Named("lst2"), Rcpp::Named("lst3"), Rcpp::Named("lst4"), Rcpp::Named("lst5"));

  if(path_length >= 1) {

    auto uids1a = assemble_uids(1, r_src_uids1, r_trg_uids1, r_count_locs1, r_signs1);
    auto uids1b = assemble_uids(1, r_src_uids1_2, r_trg_uids1_2, r_count_locs1_2, r_signs1_2);

    paths1 = exec.createPathSet(uids1a.count_total_paths());
    auto zero_1 = exec.createPathSet(data_idx1a.size());
    auto input_1 = parsed_data1->select(data_idx1a);

    exec.join(uids1a, *zero_1, *input_1, *paths1);

    auto zero_2 = exec.createPathSet(data_idx1b.size());
    auto parsed_data2 = exec.createPathSet(data2.size());
    parsed_data2->load(data2);
    auto input_2 = parsed_data2->select(data_idx1b);

    auto res = exec.join(uids1b, *zero_2, *input_2, *zero_set);
    results["lst1"] = make_score_list(res);
  }

  if(path_length >= 2) {
    auto uids2 = assemble_uids(2, r_src_uids2, r_trg_uids2, r_count_locs2, r_signs2);
    paths2 = exec.createPathSet(uids2.count_total_paths());
    auto input = parsed_data1->select(data_idx2);
    Timer timer(exec, 2, paths2->size);
    auto res = exec.join(uids2, *paths1, *input, *paths2);
    results["lst2"] = make_score_list(res);
  }

  if(path_length >= 3) {
    auto uids3 = assemble_uids(3, r_src_uids3, r_trg_uids3, r_count_locs3, r_signs3);
    paths3 = exec.createPathSet(uids3.count_total_paths());
    auto input = parsed_data1->select(data_idx3);
    Timer timer(exec, 3, paths3->size);
    auto res = exec.join(uids3, *paths2, *input, *paths3);
    results["lst3"] = make_score_list(res);
  }

  if(path_length >= 4) {
    auto uids4 = assemble_uids(4, r_src_uids4, r_trg_uids4, r_count_locs4, r_signs4);
    Timer timer(exec, 4, uids4.count_total_paths());
    auto res = exec.join(uids4, *paths3, *paths2, *zero_set);
    results["lst4"] = make_score_list(res);
  }

  if(path_length >= 5) {
    auto uids5 = assemble_uids(5, r_src_uids5, r_trg_uids5, r_count_locs5, r_signs5);
    Timer timer(exec, 5, uids5.count_total_paths());
    auto res = exec.join(uids5, *paths3, *paths3, *zero_set);
    results["lst5"] = make_score_list(res);
  }

  std::cout << "[success]" << std::endl;

  return results;
}
