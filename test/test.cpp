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


  path_set q = path_set(10);
  path_set * p = &q;
  p->len = -8900;
  printf("a: %d\n", q.len);
  printf("b: %d\n", p->len);
  

  printf("v: %d, %d\n", q.v.size(), q.v.capacity());

  path_set s = q;

  printf("q : %d\n", &(q.v));
  printf("q: %d\n", &(q.v[0]));
  printf("q: %d\n", &(q.v[1]));

  printf("q: %d\n", &(s.v));
  printf("q: %d\n", &(s.v[0]));
  printf("q: %d\n", &(s.v[1]));

  printf("v: %d, %d\n", s.v.size(), s.v.capacity());
  s.v.reserve(20481);
  printf("v: %d, %d\n", s.v.size(), s.v.capacity());


  printf("q: %d\n", &(s.v));
  printf("q: %d\n", &(s.v[0]));
  printf("q: %d\n", &(s.v[1]));
  s.v.push_back(9);
  printf("q: %d\n", &(s.v));
  printf("q: %d\n", &(s.v[0]));
  printf("q: %d\n", &(s.v[1]));

  printf("done\n");
}