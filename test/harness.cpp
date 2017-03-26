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


  // printf("generating simplified case masks for permutation (it = %d)\n", conf.iterations);
  // srand((unsigned) time(0));

  // vec2d_u16 permute_cases;
  // for(int k = 0; k < conf.iterations; k++){
  //   set<int> cases;
  //   permute_cases.push_back(vector<uint16_t>());
  //   while(cases.size() < conf.num_cases)
  //     cases.insert(rand() % (conf.num_cases + conf.num_controls));
  //   for(auto c : cases)
  //     permute_cases.back().push_back(c);
  // }
  // printf("permuted_cases: %lu (%lu)\n", permute_cases.size(), permute_cases.back().size());


  // PARSE CASE CONTROL
  // std::vector<std::vector<uint64_t> > CaseORControl2 = parseCaseORControl(CaseORControl, nCases, nControls);
  // uint64_t CaseORControl22[iterations*(vlen+vlen2)] __attribute__ ((aligned (16)));

  // CREATE DATA SETS
/*
  std::vector<std::vector<uint64_t> > parsed_data = parseData(data, nCases, nControls, vlen, vlen2);
  std::vector<std::vector<uint64_t> > parsed_data2 = parseData(data2, nCases, nControls, vlen, vlen2);
*/
  // -- FIST ITERATION -- 

  // CREATE PATH SETS
 /*
  List lst1;
  int total_paths = getTotalPaths(trguids1, uids_CountLoc1);
  std::vector<std::vector<uint64_t> > paths_pos1 = getZeroMatrix(total_paths, vlen + vlen2);
  std::vector<std::vector<uint64_t> > paths_neg1;
  if(method == 2) paths_neg1 = getZeroMatrix(total_paths, vlen + vlen2);
  std::vector<std::vector<uint64_t> > paths_pos0 = getZeroMatrix(data_inds1.size(), vlen + vlen2);
  std::vector<std::vector<uint64_t> > paths_neg0;
  if(method == 2) paths_neg0 = getZeroMatrix(data_inds1.size(), vlen + vlen2);
  std::vector<std::vector<uint64_t> > paths_data_pos = matchData(parsed_data, data_inds1);
  std::vector<std::vector<uint64_t> > paths_data_neg;
  if(method == 2) paths_data_neg = getZeroMatrix(data_inds1.size(), vlen + vlen2);
  int nnodes1 = unique(srcuids1).size();

  std::cout << nnodes1 << " source nodes for paths of length " << 1 << " and their permutations will be processed!" << std::endl;

  lst1 = JoinPaths(paths_pos0, paths_neg0, paths_data_pos, paths_data_neg, paths_pos1, paths_neg1,
   srcuids1, trguids1, uids_CountLoc1, joining_gene_sign1, ValueTable2, CaseORControl22,
   nCases, nControls, K, iterations, method, 1, nthreads);
*/
/*
  List lst1_2;

  total_paths = getTotalPaths(trguids1_2, uids_CountLoc1_2);
  std::vector<std::vector<uint64_t> > paths_pos1_2 = getZeroMatrix(total_paths, vlen + vlen2);
  std::vector<std::vector<uint64_t> > paths_neg1_2;
  if(method == 2) paths_neg1_2 = getZeroMatrix(total_paths, vlen + vlen2);
  std::vector<std::vector<uint64_t> > paths_pos0_2 = getZeroMatrix(data_inds1_2.size(), vlen + vlen2);
  std::vector<std::vector<uint64_t> > paths_neg0_2;
  if(method == 2) paths_neg0_2 = getZeroMatrix(data_inds1_2.size(), vlen + vlen2);
  paths_data_pos = matchData(parsed_data2, data_inds1_2);
  if(method == 2) paths_data_neg = getZeroMatrix(data_inds1_2.size(), vlen + vlen2);
  nnodes1 = unique(srcuids1_2).size();

  std::cout << nnodes1 << " source nodes for paths of length " << 1 << " and their permutations will be processed!" << std::endl;

  lst1_2 = JoinPaths(paths_pos0_2, paths_neg0_2, paths_data_pos, paths_data_neg, paths_pos1_2, paths_neg1_2,
   srcuids1_2, trguids1_2, uids_CountLoc1_2, joining_gene_sign1_2, ValueTable2, CaseORControl22,
   nCases, nControls, K, iterations, method, 1, nthreads);


  // Process Paths of length 2
  List lst2;
  std::vector<std::vector<uint64_t> > paths_pos2;
  std::vector<std::vector<uint64_t> > paths_neg2;

  if(pathLength >= 2){

    total_paths = getTotalPaths(trguids2, uids_CountLoc2);
    paths_pos2 = getZeroMatrix(total_paths, vlen + vlen2);
    if(method == 2) paths_neg2 = getZeroMatrix(total_paths, vlen + vlen2);
    paths_data_pos = matchData(parsed_data, data_inds2);
    if(method == 2) paths_data_neg = getZeroMatrix(data_inds2.size(), vlen + vlen2);
    int nnodes2 = unique(srcuids2).size();

    std::cout << nnodes2 << " source nodes for paths of length " << 2 << " and their permutations will be processed!" << std::endl;

    lst2 = JoinPaths(paths_pos1, paths_neg1, paths_data_pos, paths_data_neg, paths_pos2, paths_neg2,
     srcuids2, trguids2, uids_CountLoc2, joining_gene_sign2, ValueTable2,
     CaseORControl22, nCases, nControls, K, iterations, method, 2, nthreads);
  }

  // Process Paths of length 3
  List lst3;
  std::vector<std::vector<uint64_t> > paths_pos3;
  std::vector<std::vector<uint64_t> > paths_neg3;

  if(pathLength >= 3){

    total_paths = getTotalPaths(trguids3, uids_CountLoc3);
    paths_pos3 = getZeroMatrix(total_paths, vlen + vlen2);
    if(method == 2) paths_neg3 = getZeroMatrix(total_paths, vlen + vlen2);
    paths_data_pos = matchData(parsed_data, data_inds3);
    if(method == 2) paths_data_neg = getZeroMatrix(data_inds3.size(), vlen + vlen2);
    int nnodes3 = unique(srcuids3).size();

    std::cout << nnodes3 << " source nodes for paths of length " << 3 << " and their permutations will be processed!" << std::endl;

    lst3 = JoinPaths(paths_pos2, paths_neg2, paths_data_pos, paths_data_neg, paths_pos3, paths_neg3,
     srcuids3, trguids3, uids_CountLoc3, joining_gene_sign3, ValueTable2,
     CaseORControl22, nCases, nControls, K, iterations, method, 3, nthreads);
  }

  // Process Paths of length 4
  List lst4;

  if(pathLength >= 4){

    std::vector<std::vector<uint64_t> > paths_pos4;
    std::vector<std::vector<uint64_t> > paths_neg4;
    int nnodes4 = unique(srcuids4).size();

    std::cout << nnodes4 << " source nodes for paths of length " << 4 << " and their permutations will be processed!" << std::endl;

    lst4 = JoinPaths(paths_pos3, paths_neg3, paths_pos2, paths_neg2, paths_pos4, paths_neg4,
     srcuids4, trguids4, uids_CountLoc4, joining_gene_sign4, ValueTable2,
     CaseORControl22, nCases, nControls, K, iterations, method, 4, nthreads);
  }

  // Process Paths of length 5
  List lst5;

  if(pathLength >= 5){

    std::vector<std::vector<uint64_t> > paths_pos5;
    std::vector<std::vector<uint64_t> > paths_neg5;
    int nnodes5 = unique(srcuids5).size();

    std::cout << nnodes5 << " source nodes for paths of length " << 5 << " and their permutations will be processed!" << std::endl;

    lst5 = JoinPaths(paths_pos3, paths_neg3, paths_pos3, paths_neg3, paths_pos5, paths_neg5,
     srcuids5, trguids5, uids_CountLoc5, joining_gene_sign5, ValueTable2,
     CaseORControl22, nCases, nControls, K, iterations, method, 5, nthreads);

  }

  return List::create(Named("lst1") = lst1_2,
    Named("lst2") = lst2,
    Named("lst3") = lst3,
    Named("lst4") = lst4,
    Named("lst5") = lst5);

*/

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
