#include <typeinfo>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <set>
#include <time.h>
#include "gcre.h"

using namespace std;

static char* get_opt(char** begin, char** end, const string& option) {
  char** itr = find(begin, end, option);
  if (itr != end && ++itr != end) {
    return *itr;
  }
  return 0;
}

static bool has_opt(char** begin, char** end, const string& option) {
  return find(begin, end, option) != end;
}

static void read_vec(istream& fdata, vec2d_u64& vec){
  string line, val;

  getline(fdata, line, ':');
  cout << line << ": ";
  getline(fdata, line, ' ');
  int width = stoi(line);
  getline(fdata, line);
  stringstream buf(line);
  while (getline(buf, val, ' ')) {
    if(vec.size() == 0 || vec.back().size() == width)
      vec.push_back(vec_u64());
    vec.back().push_back(stoul(val));
  }
  cout << vec.size() << endl;
}

int main(int argc, char* argv[]) {
  printf("hello\n\n");
  PathSet_BlockM1 p(8, 4);
  for(int k = 0; k < p.size; k++) {
    for(int v = 0; v < p.width_ul; v++){
      printf(" %lu", p[k][v]);
    }
    printf("\n");
  }
  printf("\n");

  PathSet_BlockM2 p2(8, 4);
  printf("\ndone\n");
  return 0;
}

/*
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

  join_config conf;
  conf.top_k = 4;
  conf.iterations = 10;
  conf.nthreads = 1;

  if(has_opt(argv, arge, "-k"))
    conf.top_k = stoi(get_opt(argv, arge, "-k"));

  if(has_opt(argv, arge, "-p"))
    conf.iterations = max(1, stoi(get_opt(argv, arge, "-p")));

  if(has_opt(argv, arge, "-t"))
    conf.nthreads = stoi(get_opt(argv, arge, "-t"));

  int repeat = 5;
  if(has_opt(argv, arge, "-r"))
    repeat = stoi(get_opt(argv, arge, "-r"));

  // PATH 4
  // CASE 1537
  // CTRL 1537
  // TOTL 1276368

  int width;
  string line, val, v;
  stringstream buf, ub;
  ifstream fdata(file);

  getline(fdata, line);
  conf.path_length = stoi(line.substr(5));
  getline(fdata, line);
  conf.num_cases = stoi(line.substr(5));
  getline(fdata, line);
  conf.num_controls = stoi(line.substr(5));
  getline(fdata, line);
  uint64_t total_paths = stoul(line.substr(5));

  printf("\nheader:\n\n");
  printf("   path_length : %d\n", conf.path_length);
  printf("     num_cases : %d\n", conf.num_cases);
  printf("  num_controls : %d\n", conf.num_controls);
  printf("   total_paths : %lu\n", total_paths);
  printf("\n");

  vector<uid_ref> uids;
  getline(fdata, line, ' ');
  cout << line << ": ";
  getline(fdata, line);
  buf.str(line);
  buf.clear();
  while (getline(buf, val, ' ')) {
    uid_ref uid;
    ub.str(val);
    ub.clear();
    getline(ub, v, ':');
    uid.src = stoi(v);
    getline(ub, v, ':');
    uid.trg = stoi(v);
    getline(ub, v, ':');
    uid.count = stoi(v);
    getline(ub, v, ':');
    uid.location = stoi(v);
    uids.push_back(uid);
  }
  cout << uids.size() << endl;

  vector<int> signs;
  getline(fdata, line, ' ');
  cout << line << ": ";
  getline(fdata, line);
  buf.str(line);
  buf.clear();
  while (getline(buf, val, ' ')) {
    signs.push_back(stoi(val));
  }
  cout << signs.size() << endl;

  vec2d_d value_table;
  getline(fdata, line, ':');
  cout << line << ": ";
  getline(fdata, line, ' ');
  width = stoi(line);
  getline(fdata, line);
  buf.str(line);
  buf.clear();
  while (getline(buf, val, ' ')) {
    if(value_table.size() == 0 || value_table.back().size() == width)
      value_table.push_back(vector<double>());
    value_table.back().push_back(stod(val));
  }
  cout << value_table.size() << endl;

  getline(fdata, line, ':');
  cout << line << ": ";
  getline(fdata, line, ' ');
  width = stoi(line);
  getline(fdata, line);
  buf.str(line);
  buf.clear();
  // don't need to read the permuted case mask
  while (getline(buf, val, ' ')) {
  }
  cout << endl;

  vec2d_u64 paths0pos;
  read_vec(fdata, paths0pos);

  vec2d_u64 paths0neg;
  read_vec(fdata, paths0neg);

  vec2d_u64 paths1pos;
  read_vec(fdata, paths1pos);

  vec2d_u64 paths1neg;
  read_vec(fdata, paths1neg);

  JoinMethod2Native method;

  paths_type* paths0 = method.createPathSet(paths0pos, paths0neg, conf.num_cases, conf.num_controls);
  paths_type* paths1 = method.createPathSet(paths1pos, paths1neg, conf.num_cases, conf.num_controls);

  joined_res tres;
  joined_res& res = tres;

  printf("generating simplified case masks for permutation (it = %d)\n", conf.iterations);
  srand((unsigned) time(0));

  vec2d_u16 permute_cases;
  for(int k = 0; k < conf.iterations; k++){
    set<int> cases;
    permute_cases.push_back(vector<uint16_t>());
    while(cases.size() < conf.num_cases)
      cases.insert(rand() % (conf.num_cases + conf.num_controls));
    for(auto c : cases)
      permute_cases.back().push_back(c);
  }
  printf("permuted_cases: %lu (%lu)\n", permute_cases.size(), permute_cases.back().size());

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
  
  printf("done\n");
}
*/