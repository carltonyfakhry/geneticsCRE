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

UidRelSet add_signs(int path_length, vector<uid_ref>& uids, vector<int> signs){
  add_idx(uids);
  return UidRelSet(path_length, uids, signs);
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

  int iters = 10;
  if(has_opt(argv, arge, "-p"))
    iters = max(0, stoi(get_opt(argv, arge, "-p")));

  string method = "method2";
  if(has_opt(argv, arge, "-m"))
    method = get_opt(argv, arge, "-m");

  JoinExec exec(method, num_cases, num_ctrls, iters);

  if(has_opt(argv, arge, "-l"))
    path_length = stoi(get_opt(argv, arge, "-l"));

  if(has_opt(argv, arge, "-k"))
    exec.top_k = stoi(get_opt(argv, arge, "-k"));

  if(has_opt(argv, arge, "-t"))
    exec.nthreads = stoi(get_opt(argv, arge, "-t"));

  int repeat = 1;
  if(has_opt(argv, arge, "-r"))
    repeat = stoi(get_opt(argv, arge, "-r"));

  int width;

  exec.print_vector_info();

  printf("\nheader:\n\n");
  printf("   path_length : %d\n", path_length);
  printf("         cases : %d\n", exec.num_cases);
  printf("      controls : %d\n", exec.num_ctrls);
  printf("\n");

  auto uids0v = read_uids(fdata);
  auto sign0 = read_ints(fdata);
  auto uids0 = add_signs(1, uids0v, sign0);
  
  auto uids1v = read_uids(fdata);
  auto sign1 = read_ints(fdata);
  auto uids1 = add_signs(1, uids1v, sign1);
  
  auto uids2v = read_uids(fdata);
  auto sign2 = read_ints(fdata);
  auto uids2 = add_signs(2, uids2v, sign2);
  
  auto uids3v = read_uids(fdata);
  auto sign3 = read_ints(fdata);
  auto uids3 = add_signs(3, uids3v, sign3);
  
  auto uids4v = read_uids(fdata);
  auto sign4 = read_ints(fdata);
  auto uids4 = add_signs(4, uids4v, sign4);
  
  auto uids5v = read_uids(fdata);
  auto sign5 = read_ints(fdata);
  auto uids5 = add_signs(5, uids5v, sign5);

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

  Timer::print_header();

  if(path_length >= 1) {

    paths1 = exec.createPathSet(uids0.count_total_paths());
    auto zero_1 = exec.createPathSet(data_idx0.size());
    auto input_1 = parsed_data1->select(data_idx0);

    auto res0 = exec.join(uids0, *zero_1, *input_1, *paths1);

    auto zero_2 = exec.createPathSet(data_idx1.size());
    auto parsed_data2 = exec.createPathSet(data2.size());
    parsed_data2->load(data2);
    auto input_2 = parsed_data2->select(data_idx1);

    auto res1 = exec.join(uids1, *zero_2, *input_2, *zero_set);

  }

  if(path_length >= 2) {
    paths2 = exec.createPathSet(uids2.count_total_paths());
    auto input = parsed_data1->select(data_idx2);
    Timer timer(exec, 2, paths2->size);
    auto res2 = exec.join(uids2, *paths1, *input, *paths2);
  }

  if(path_length >= 3) {
    paths3 = exec.createPathSet(uids3.count_total_paths());
    auto input = parsed_data1->select(data_idx3);
    Timer timer(exec, 3, paths3->size);
    auto res3 = exec.join(uids3, *paths2, *input, *paths3);
  }

  if(path_length >= 4) {

    joined_res res4;
    for(int r = 0; r < repeat; r++){
      Timer timer(exec, 4, uids4.count_total_paths());
      res4 = exec.join(uids4, *paths3, *paths2, *zero_set);
    }

    printf("\n################################\n");
    printf("   length : %d  width: %d  iters: %d  thread: %d\n", path_length, exec.width_ul, exec.iterations, exec.nthreads);
    printf("  results : %lu |", res4.scores.size());
    for(int k = 0; k < res4.scores.size(); k++){
      Score s = res4.scores[k];
      printf(" %f[%d:%d]", s.score, s.src, s.trg);
    }
    printf("\n");
    printf("    perms :");
    for(int k = 0; k < min(24, (int) res4.permuted_scores.size()); k++)
      printf(" %0.2f", res4.permuted_scores[k]);
    printf("\n");

    printf("################################\n");


  }

  if(path_length >= 5) {
    Timer timer(exec, 5, uids5.count_total_paths());
    auto res5 = exec.join(uids5, *paths3, *paths3, *zero_set);
  }

  printf("done\n");
}
