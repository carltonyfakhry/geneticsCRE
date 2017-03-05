#include <typeinfo>
#include "gcre.h"

using namespace std;

int main(int argc, char* argv[]) {
  JoinMethod2Native m;
  std::cout << "method: " << typeid(m).name() << endl;
  printf("done\n");
}