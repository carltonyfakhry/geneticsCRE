#include "Utils.h"

using namespace Rcpp;



/**
 *
 * Convert integer to string
 *
 */

std::string IntToString(int a){
  std::ostringstream temp;
  temp<<a;
  return temp.str();
}



/**
 *
 * Convert string to uint64_t
 *
 */
uint64_t StringToInt(std::string str){
  uint64_t value;
  std::istringstream iss(str);
  iss >> value;
  return value;
}



/**
 *
 * Create list indexed by the names of gene uids. Each entry in the list
 * is a vector with the first entry being the count of paths starting
 * with the corresponding gene uids and the second entry be the location
 * of the first path in the relations data frame.
 *
 */
// [[Rcpp::export]]
Rcpp::List getMatchingList(Rcpp::IntegerVector uids, Rcpp::IntegerVector counts, Rcpp::IntegerVector location){

  Rcpp::List uids_2countsloc = Rcpp::List();

  for(unsigned int i = 0; i < uids.size(); i++){

    std::string uids_str = IntToString(uids[i]);
    Rcpp::IntegerVector temp = Rcpp::IntegerVector(2);
    int count = counts[i];
    temp[0] = count;
    temp[1] = location[i];
    uids_2countsloc[uids_str] = temp;

  }

  return uids_2countsloc;

}

