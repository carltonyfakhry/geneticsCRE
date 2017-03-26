#include <typeinfo>
#include <sstream>
#include <cstdlib>
#include "test.h"

using namespace std;

namespace test {

  char* get_opt(char** begin, char** end, const string& option) {
    char** itr = find(begin, end, option);
    if (itr != end && ++itr != end) {
      return *itr;
    }
    return 0;
  }

  bool has_opt(char** begin, char** end, const string& option) {
    return find(begin, end, option) != end;
  }

  string read_line(ifstream& fdata){
    string line;
    getline(fdata, line, ' ');
    cout << line << " (";
    getline(fdata, line);
    cout << line.size() << ")" << endl;
    return line;
  }

  vector<uid_ref> read_uids(ifstream& fdata){
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

  vector<int> read_ints(ifstream& fdata){
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

  vec2d_i read_data(ifstream& fdata){
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

  vec2d_d read_vals(ifstream& fdata){
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

}