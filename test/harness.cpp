#include <chrono>

#include "gcre.h"
#include "test.h"

using namespace std;
using namespace test;

void add_idx(vector<uid_ref>& refs) {
  int idx = 0;
  for(auto& ref : refs) {
    ref.path_idx = idx;
    idx += ref.count;
  }
}

void add_signs(int path_length, vector<uid_ref>& uids, vector<int> signs){
  for(int idx = 0; idx < uids.size(); idx++){
    auto& uid = uids[idx];
    for(int loc = uid.location; loc < uid.location + uid.count; loc++) {
      int sign = 0;
      if(path_length > 3)
        sign = signs[idx];
      else if(path_length < 3)
        sign = signs[loc];
      else if(path_length == 3)
        sign = (signs[idx] + signs[loc] == 0) ? -1 : 1;
      uid.signs.push_back(sign == 1);
    }
  }
}

int main(int argc, char* argv[]) {

  setlocale(LC_ALL, "");

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

  if(has_opt(argv, arge, "-l"))
    path_length = stoi(get_opt(argv, arge, "-l"));

  if(has_opt(argv, arge, "-k"))
    exec.top_k = stoi(get_opt(argv, arge, "-k"));

  if(has_opt(argv, arge, "-p"))
    exec.iterations = max(0, stoi(get_opt(argv, arge, "-p")));

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

  auto uids0 = read_uids(fdata);
  auto sign0 = read_ints(fdata);
  add_signs(1, uids0, sign0);
  add_idx(uids0);
  
  auto uids1 = read_uids(fdata);
  auto sign1 = read_ints(fdata);
  add_signs(1, uids1, sign1);
  add_idx(uids1);
  
  auto uids2 = read_uids(fdata);
  auto sign2 = read_ints(fdata);
  add_signs(2, uids2, sign2);
  add_idx(uids2);
  
  auto uids3 = read_uids(fdata);
  auto sign3 = read_ints(fdata);
  add_signs(3, uids3, sign3);
  add_idx(uids3);
  
  auto uids4 = read_uids(fdata);
  auto sign4 = read_ints(fdata);
  add_signs(4, uids4, sign4);
  add_idx(uids4);
  
  auto uids5 = read_uids(fdata);
  auto sign5 = read_ints(fdata);
  add_signs(5, uids5, sign5);
  add_idx(uids5);

  auto data_idx0 = read_ints(fdata);
  auto data_idx1 = read_ints(fdata);
  auto data_idx2 = read_ints(fdata);
  auto data_idx3 = read_ints(fdata);

  auto data1 = read_data(fdata);
  auto data2 = read_data(fdata);
  auto perms = read_data(fdata);
  auto table = read_vals(fdata);

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

    auto res0 = exec.join(uids0, *zero_1, *input_1, *paths1);

    auto zero_2 = exec.createPathSet(data_idx1.size());
    auto parsed_data2 = exec.createPathSet(data2.size());
    parsed_data2->load(data2);
    auto input_2 = parsed_data2->select(data_idx1);

    auto res1 = exec.join(uids1, *zero_2, *input_2, *zero_set);

  }

  if(path_length >= 2) {
    paths2 = exec.createPathSet(JoinExec::count_total_paths(uids2));
    auto input = parsed_data1->select(data_idx2);
    auto res2 = exec.join(uids2, *paths1, *input, *paths2);
  }

  if(path_length >= 3) {
    paths3 = exec.createPathSet(JoinExec::count_total_paths(uids3));
    auto input = parsed_data1->select(data_idx3);
    auto res3 = exec.join(uids3, *paths2, *input, *paths3);
  }

  if(path_length >= 4) {

    auto res4 = exec.join(uids4, *paths3, *paths2, *zero_set);

    printf("\n################################\n");
    printf("   length : %d  width: %d  iters: %d  thread: %d\n", path_length, exec.width_ul, exec.iterations, exec.nthreads);
    printf("  results : %lu |", res4.scores.size());
    for(int k = 0; k < res4.scores.size(); k++){
      Score s = res4.scores[k];
      printf(" %f[%d:%d]", s.score, s.src, s.trg);
    }
    printf("\n");
    printf("    perms :");
    for(auto ps : res4.permuted_scores)
      printf(" %0.2f", ps);
    printf("\n");

    printf("################################\n\n");


  }

  if(path_length >= 5) {

    auto start = chrono::system_clock::now();

    auto res5 = exec.join(uids5, *paths3, *paths3, *zero_set);

    auto time = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - start);
    printf("[%d] time: %'ld ms (%'ld x threads)\n", 0, time.count(), time.count() * max(1, exec.nthreads));

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

  */
  printf("done\n");
}
