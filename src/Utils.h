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

using namespace Rcpp;

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
