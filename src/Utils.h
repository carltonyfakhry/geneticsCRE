#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <Rmath.h>
#include <utility>
#include <iterator>
#include <sstream>
#include <queue>
#include <vector>
#include <bitset>

// for timer
#include <ctime>
#include <chrono>
#include <unistd.h>


using namespace Rcpp;


#define ROUND_DOWN(x, s) ((x) & ~((s)-1))



const uint64_t one_64bit = 1;


typedef std::pair<double,std::pair<int, int> > ScoreIndices;
typedef std::vector<std::priority_queue<double, std::vector<double>, std::greater<double> > > VectorPriorityQueues;



class CompareScoreIndices
{
public:
  bool operator()(ScoreIndices n1, ScoreIndices n2) {
    return n1.first>n2.first;
  }
};


class CompareScore
{
public:
  bool operator()(double score1, double score2) {
    return score1>score2;
  }
};



typedef std::priority_queue<ScoreIndices, std::vector<ScoreIndices>, CompareScoreIndices> IndicesScoresQueue;



/**
 *
 * Convert integer to string
 *
 */
std::string IntToString(int a);



/**
 *
 * Convert string to uint64_t
 *
 */
uint64_t StringToInt(std::string str);



/**
 *
 * Create list indexed by the names of gene uids. Each entry in the list
 * is a vector with the first entry being the count of paths starting
 * with the corresponding gene uids and the second entry be the location
 * of the first path in the relations data frame.
 *
 */
Rcpp::List getMatchingList(Rcpp::IntegerVector uids, Rcpp::IntegerVector counts, Rcpp::IntegerVector location);


// timer class for benchmarking --dmitri

using namespace std;

class Timer {

public:

  static void print_header() {
    printf("\nTIME:PID IMPL METHOD WIDTH LENGTH PATHS PERMS THREADS CYCLES MS\n\n");
  }

  Timer(int method, int width_ul, int path_length, int iters, int nthreads) : 
  method(method),
  width_ul(width_ul),
  path_length(path_length),
  iters(iters),
  nthreads(nthreads)
  {}

  ~Timer(){
    auto time = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - start);
    printf("\nTIME:%d CPU-OMP m%d %d %d %lu %d %d %lu %lu\n\n", getpid(), method, width_ul * 64,
      path_length, total_paths, iters, nthreads, std::clock() - clock_start, (unsigned long) time.count());
  }

  uint64_t total_paths = 0;

private:

  const clock_t clock_start = std::clock();
  const chrono::system_clock::time_point start = chrono::system_clock::now();
  const int method;
  const int width_ul;
  const int path_length;
  const int iters;
  const int nthreads;

};