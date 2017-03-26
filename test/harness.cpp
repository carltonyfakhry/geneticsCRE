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

static string read_line(ifstream& fdata){
  string line;
  getline(fdata, line, ' ');
  cout << line << " (";
  getline(fdata, line);
  cout << line.size() << ")" << endl;
  return line;
}

static vector<uid_ref> read_uids(ifstream& fdata){
  string line, val, v;
  stringstream buf, ub;
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
  return uids;
}

static vector<int> read_ints(ifstream& fdata){
  string line, val;
  stringstream buf;
  vector<int> data;
  getline(fdata, line, ' ');
  cout << line << ": ";
  getline(fdata, line);
  buf.str(line);
  buf.clear();
  while (getline(buf, val, ' ')) {
    data.push_back(stoi(val));
  }
  cout << data.size() << endl;
  return data;
}

static vec2d_i read_data(ifstream& fdata){
  string line, val, v;
  stringstream buf, ub;
  vec2d_i data;

  getline(fdata, line, ' ');
  cout << line << ": ";
  getline(fdata, line);
  buf.str(line);
  buf.clear();

  while (getline(buf, val, ' ')) {
    data.push_back(vector<int>());
    ub.str(val);
    ub.clear();
    while (getline(ub, v, ','))
      data.back().push_back(stoi(v));
  }
  cout << data.size() << " x " << data.front().size() << endl;
  return data;
}


static vec2d_d read_vals(ifstream& fdata){
  string line, val, v;
  stringstream buf, ub;
  vec2d_d data;

  getline(fdata, line, ' ');
  cout << line << ": ";
  getline(fdata, line);
  buf.str(line);
  buf.clear();

  while (getline(buf, val, ' ')) {
    data.push_back(vector<double>());
    ub.str(val);
    ub.clear();
    while (getline(ub, v, ','))
      data.back().push_back(stod(v));
  }
  cout << data.size() << " x " << data.front().size() << endl;
  return data;
}

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

  int width;
  string line, val, v;
  stringstream buf, ub;

  ifstream fdata(file);

  conf.path_length = stoi(read_line(fdata));
  conf.num_cases = stoi(read_line(fdata));
  conf.num_controls = stoi(read_line(fdata));

  printf("\nheader:\n\n");
  printf("   path_length : %d\n", conf.path_length);
  printf("     num_cases : %d\n", conf.num_cases);
  printf("  num_controls : %d\n", conf.num_controls);
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

/*
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
  */
  printf("done\n");
}
