#include <typeinfo>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <set>
#include <time.h>

#include "gcre.h"
#include "test.h"

using namespace std;
using namespace test;

int main(int argc, char* argv[]) {

  char** arge = argc + argv;
  
  string file;
  if(has_opt(argv, arge, "-f")){
    file = string(get_opt(argv, arge, "-f"));
  } else {
    cout << "need data file (-f)" << endl;
    return 1;
  }
  cout << "using data file: " << file << endl;

  ifstream fdata(file);

  int path_length = stoi(read_line(fdata));
  int num_cases = stoi(read_line(fdata));
  int num_ctrls = stoi(read_line(fdata));

  JoinExec exec(num_cases, num_ctrls);

  if(has_opt(argv, arge, "-k"))
    exec.top_k = stoi(get_opt(argv, arge, "-k"));

  if(has_opt(argv, arge, "-p"))
    exec.iterations = max(1, stoi(get_opt(argv, arge, "-p")));

  if(has_opt(argv, arge, "-t"))
    exec.nthreads = stoi(get_opt(argv, arge, "-t"));

  int repeat = 5;
  if(has_opt(argv, arge, "-r"))
    repeat = stoi(get_opt(argv, arge, "-r"));

  int width;

  printf("\nheader:\n\n");
  printf("   path_length : %d\n", path_length);
  printf("         cases : %d\n", exec.num_cases);
  printf("      controls : %d\n", exec.num_ctrls);
  printf("\n");

  vector<uid_ref> uids0 = read_uids(fdata); vector<int> signs0 = read_ints(fdata);
  vector<uid_ref> uids1 = read_uids(fdata); vector<int> signs1 = read_ints(fdata);
  vector<uid_ref> uids2 = read_uids(fdata); vector<int> signs2 = read_ints(fdata);
  vector<uid_ref> uids3 = read_uids(fdata); vector<int> signs3 = read_ints(fdata);
  vector<uid_ref> uids4 = read_uids(fdata); vector<int> signs4 = read_ints(fdata);
  vector<uid_ref> uids5 = read_uids(fdata); vector<int> signs5 = read_ints(fdata);

  vector<int> data_idx0 = read_ints(fdata);
  vector<int> data_idx1 = read_ints(fdata);
  vector<int> data_idx2 = read_ints(fdata);
  vector<int> data_idx3 = read_ints(fdata);

  vec2d_i data1 = read_data(fdata);
  vec2d_i data2 = read_data(fdata);
  vec2d_i perms = read_data(fdata);
  vec2d_d table = read_vals(fdata);

  exec.setValueTable(table);
  exec.setPermutedCases(perms);

  auto parsed_data1 = exec.createPathSet(data1.size());
  parsed_data1->load(data1);

  auto zero_set = exec.createPathSet(0);

  unique_ptr<PathSet> paths1 = nullptr, paths2 = nullptr, paths3 = nullptr;

  if(path_length >= 1) {

    paths1 = exec.createPathSet(JoinExec::count_total_paths(uids0));
    auto zero_1 = exec.createPathSet(data_idx0.size());
    auto input_1 = parsed_data1->select(data_idx1);

    auto res0 = exec.join(uids0, signs0, *zero_1, *input_1, *paths1);

    auto zero_2 = exec.createPathSet(data_idx1.size());
    auto parsed_data2 = exec.createPathSet(data2.size());
    parsed_data2->load(data2);
    auto input_2 = parsed_data2->select(data_idx1);

    auto res1 = exec.join(uids1, signs1, *zero_2, *input_2, *zero_set);

  }

  if(path_length >= 2) {
    paths2 = exec.createPathSet(JoinExec::count_total_paths(uids2));
    auto input = parsed_data1->select(data_idx2);
    auto res2 = exec.join(uids2, signs2, *paths1, *input, *paths2);
  }

  if(path_length >= 3) {
    paths3 = exec.createPathSet(JoinExec::count_total_paths(uids3));
    auto input = parsed_data1->select(data_idx3);
    auto res3 = exec.join(uids3, signs3, *paths2, *input, *paths3);
  }

  if(path_length >= 4) {
    auto res4 = exec.join(uids4, signs4, *paths3, *paths2, *zero_set);
  }

  if(path_length >= 5) {
    auto res5 = exec.join(uids5, signs5, *paths3, *paths3, *zero_set);
  }

/*
  printf("\n");
  for(int r = 0; r < repeat; r++){
    clock_t time = clock();
    res = method.join(conf, uids, signs, value_table, permute_cases, paths0, paths1, NULL, total_paths);
    time = clock() - time;
    printf("[%d] time: %lu (%f ms)\n", r, time, time / (double) CLOCKS_PER_SEC * 1000);
  }
  printf("\n");

  printf("\n################################\n");
  printf("   length : %d  paths: %lu  width: %d  iter: %d  thread: %d\n", conf.path_length, total_paths, 0, conf.iterations, conf.nthreads);
  printf("  results : %lu |", res.scores.size());
  for(int k = 0; k < res.scores.size(); k++){
    Score s = res.scores[k];
    printf(" %f[%d:%d]", s.score, s.src, s.trg);
  }
  printf("\n");
  printf("    perms :");
  for(auto ps : res.permuted_scores)
    printf(" %0.2f", ps);
  printf("\n");

  printf("################################\n\n");
  */
  printf("done\n");
}
