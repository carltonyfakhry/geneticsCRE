#include <fstream>

#include <inttypes.h>
#include <stdio.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <cmath>
#include <limits>
#include <vector>
#include <set>
#include <queue>
#include <cstdio>
#include <iostream>
#include <fstream>

#include <Rcpp.h>

#include "gcre.h"


// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;


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

int getTotalPaths(IntegerVector trguids, List uids_CountLoc){

  int total_paths = 0;

  for(int i = 0; i < trguids.size(); i++){

    int uid = trguids[i];
    std::string geneuid = IntToString(uid);
    IntegerVector temp = as<IntegerVector>(uids_CountLoc[geneuid]);
    total_paths += temp[0];

  }

  return total_paths;

}

int getTotalCountsCountLoc(List uids_CountLoc){

  int total_paths = 0;

  for(List::iterator it = uids_CountLoc.begin(); it != uids_CountLoc.end(); ++it){

    IntegerVector temp = as<IntegerVector>(*it);
    total_paths += temp[0];

  }

  return total_paths;

}

std::vector<std::vector<double> > copyValueTable(NumericMatrix ValueTable){

  std::vector<std::vector<double> > ValueTable2(ValueTable.nrow(), std::vector<double>(ValueTable.ncol()));

  for(int i = 0; i < ValueTable.nrow(); i++)
    for(int j = 0; j < ValueTable.ncol(); j++)
      ValueTable2[i][j] = ValueTable(i,j);

    return ValueTable2;

  }

/**
*
* Convert CaseORControl from 0/1 values to 64 bit representations.
*

std::vector<std::vector<uint64_t> > parseCaseORControl(IntegerMatrix CaseORControl, int nCases, int nControls){


  int vlen = (int) ceil(nCases/64.0);
  int vlen2 = (int) ceil(nControls/64.0);

  std::vector<std::vector<uint64_t> > CaseORControl2(CaseORControl.nrow(), std::vector<uint64_t>(vlen + vlen2, 0));
  for(int i = 0; i < CaseORControl.nrow(); i++){

    for(int j = 0; j < nCases; j++){
      int index = j/64;
      if(CaseORControl(i,j) == 1)
        CaseORControl2[i][index] |= one_64bit << j % 64;
    }

    for(int j = nCases; j < nCases + nControls; j++){
      int index = vlen + (j-nCases)/64;
      if(CaseORControl(i,j) == 1)
        CaseORControl2[i][index] |= one_64bit << (j-nCases) % 64;
    }

  }

  return CaseORControl2;

}
*/

/*

std::vector<std::vector<uint64_t> > parseData(IntegerMatrix data, int nCases, int nControls, int vlen, int vlen2){

  std::vector<std::vector<uint64_t> > parseddata(data.nrow(), std::vector<uint64_t>(vlen + vlen2));

  for(int i = 0; i < data.nrow(); i++){

    for(int j = 0; j < data.ncol(); j++){

      if(data(i,j) != 0){

        if(j < nCases)
          parseddata[i][j/64] |= one_64bit << j % 64;
        else
          parseddata[i][vlen + (j-nCases)/64] |= one_64bit << (j-nCases) % 64;

      }

    }

  }

  return parseddata;

}


std::vector<std::vector<uint64_t> > matchData(const std::vector<std::vector<uint64_t> > &parseddata, IntegerVector data_inds){

  std::vector<std::vector<uint64_t> > matcheddata(data_inds.size(), std::vector<uint64_t>(parseddata[0].size(),0));

  for(int i = 0; i < data_inds.size(); i++){

    int data_index = data_inds[i];
    matcheddata[i] = parseddata[data_index];

  }

  return matcheddata;

}
*/

  vector<uid_ref> make_uids(IntegerVector& r_src_uids, IntegerVector& r_trg_uids, List& r_uid_count_locs){
    long total_paths = 0;
    vector<uid_ref> uids(r_trg_uids.size(), uid_ref());
    for(int k = 0; k < r_trg_uids.size(); k++){
      uid_ref& uid = uids[k];
      uid.trg = r_trg_uids[k];
      uid.src = r_src_uids[k];
      IntegerVector count_loc = r_uid_count_locs[to_string(uid.trg)];
      uid.count = count_loc[0];
      uid.location = count_loc[1];
      total_paths += uid.count;
    }
    printf("uids: %lu (paths: %ld)\n", uids.size(), total_paths);
    return uids;
  }

  void print_uids(FILE* fp, int num, vector<uid_ref>&& uids, IntegerVector& signs){
    std::fprintf(fp, "UIDS%d", num);
    for(int k = 0; k < uids.size(); k++){
      uid_ref& u = uids[k];
      std::fprintf(fp, " %d:%d:%d:%d", u.src, u.trg, u.count, u.location);
    }
    std::fprintf(fp, "\n");
    std::fprintf(fp, "SIGN%d", num);
    for(int k = 0; k < signs.size(); k++){

      std::fprintf(fp, " %d", signs[k]);
    }
    std::fprintf(fp, "\n");
  }

  void print_dataind(FILE* fp, int num, IntegerVector& data){
    std::fprintf(fp, "DIDX%d", num);
    for(int k = 0; k < data.size(); k++){
      std::fprintf(fp, " %d", data[k]);
    }
    std::fprintf(fp, "\n");
  }

  void print_matrix(FILE* fp, string name, int num, IntegerMatrix& data){
    std::fprintf(fp, "%s%d", name.c_str(), num);
    for(int r = 0; r < data.nrow(); r++){
      std::fprintf(fp, " ");
      for(int k = 0; k < data.ncol(); k++){
        if(k != 0)
          std::fprintf(fp, ",");
        std::fprintf(fp, "%d", data(r, k));
      }
    }
    std::fprintf(fp, "\n");
  }

  void print_table(FILE* fp, NumericMatrix& data){
    std::fprintf(fp, "TABLE");
    for(int r = 0; r < data.nrow(); r++){
      std::fprintf(fp, " ");
      for(int k = 0; k < data.ncol(); k++){
        if(k != 0)
          std::fprintf(fp, ",");
        if(data(r,k) < 0.000001)
          std::fprintf(fp, "0");
        else
          std::fprintf(fp, "%lf", data(r,k));
      }
    }
    std::fprintf(fp, "\n");
  }

// [[Rcpp::export]]
  List ProcessPaths(IntegerVector srcuids1, IntegerVector trguids1, List uids_CountLoc1, IntegerVector joining_gene_sign1,
    IntegerVector srcuids1_2, IntegerVector trguids1_2, List uids_CountLoc1_2, IntegerVector joining_gene_sign1_2,
    IntegerVector srcuids2, IntegerVector trguids2, List uids_CountLoc2, IntegerVector joining_gene_sign2,
    IntegerVector srcuids3, IntegerVector trguids3, List uids_CountLoc3, IntegerVector joining_gene_sign3,
    IntegerVector srcuids4, IntegerVector trguids4, List uids_CountLoc4, IntegerVector joining_gene_sign4,
    IntegerVector srcuids5, IntegerVector trguids5, List uids_CountLoc5, IntegerVector joining_gene_sign5,
    IntegerVector data_inds1, IntegerVector data_inds1_2, IntegerVector data_inds2, IntegerVector data_inds3,
    IntegerMatrix data, IntegerMatrix data2, NumericMatrix ValueTable, int nCases, int nControls, int K,
    int iterations, IntegerMatrix CaseORControl, int method, int pathLength, int nthreads){

    printf("meh\n");

    FILE* fp = std::fopen("/data/gcre/data1.ser", "w");
    if(!fp) {
      std::perror("File opening failed");
      exit(EXIT_FAILURE);
    }

    std::fprintf(fp, "PATHS %d\n", pathLength);
    std::fprintf(fp, "CASES %d\n", nCases);
    std::fprintf(fp, "CTRLS %d\n", nControls);

    print_uids(fp, 0, make_uids(srcuids1, trguids1, uids_CountLoc1), joining_gene_sign1);
    print_uids(fp, 1, make_uids(srcuids1_2, trguids1_2, uids_CountLoc1_2), joining_gene_sign1_2);
    print_uids(fp, 2, make_uids(srcuids2, trguids2, uids_CountLoc2), joining_gene_sign2);
    print_uids(fp, 3, make_uids(srcuids3, trguids3, uids_CountLoc3), joining_gene_sign3);
    print_uids(fp, 4, make_uids(srcuids4, trguids4, uids_CountLoc4), joining_gene_sign4);
    print_uids(fp, 5, make_uids(srcuids5, trguids5, uids_CountLoc5), joining_gene_sign5);

    print_dataind(fp, 0, data_inds1);
    print_dataind(fp, 1, data_inds1_2);
    print_dataind(fp, 2, data_inds2);
    print_dataind(fp, 3, data_inds3);

    print_matrix(fp, "DATA", 1, data);
    print_matrix(fp, "DATA", 2, data2);
    print_matrix(fp, "PERM", 0, CaseORControl);

    print_table(fp, ValueTable);

    std::fclose(fp);
    printf("done\n");
    exit(0);

    return List::create();
  }
