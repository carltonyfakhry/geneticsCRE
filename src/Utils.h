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
