#ifndef GCRE_TEST_H
#define GCRE_TEST_H

#include <fstream>
#include <vector>

#include "gcre.h"

using namespace std;

namespace test {

  char* get_opt(char** begin, char** end, const string& option);
  bool has_opt(char** begin, char** end, const string& option);

  string read_line(ifstream& fdata);
  vector<uid_ref> read_uids(ifstream& fdata);
  vector<int> read_ints(ifstream& fdata);
  vec2d_i read_data(ifstream& fdata);
  vec2d_d read_vals(ifstream& fdata);

}

#endif
