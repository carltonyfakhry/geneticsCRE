#include <iostream>
#include <vector>
#include <typeinfo>
#include <stdio.h>
#include <stdint.h>

using namespace std;

struct path_set {
  int len;
  vector<int> v;
  path_set(int len);
  ~path_set();
};

struct uid_ref {
  int src;
  int trg;
  int count;
  int location;
};

path_set::path_set(int len) : len(len) {
  printf("i am created: %d, %d\n", len, len);
  v.push_back(17);
  v.push_back(17);
}

path_set::~path_set(){
  printf("i am destructed\n");
}

int main()
{

  vector<uid_ref> u(2, uid_ref());
  uid_ref& r = * new uid_ref;
  u.push_back(r);
  for(int k = 0; k < u.size(); k++)
    printf(" %d", u[k].trg);
  printf("\n");

  printf("done\n");
}