#include <fstream>
#include "Utils.h"
#include "gcre.h"

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

/**
 *
 * Get the total number of paths to be processed.
 *
 */

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



/**
 *
 * Get the total counts in uids_CountLoc.
 *
 */

int getTotalCountsCountLoc(List uids_CountLoc){

  int total_paths = 0;

  for(List::iterator it = uids_CountLoc.begin(); it != uids_CountLoc.end(); ++it){

    IntegerVector temp = as<IntegerVector>(*it);
    total_paths += temp[0];

  }

  return total_paths;

}



/**
 *
 * This function returns a 2D vector filled with zeros.
 *
 */

std::vector<std::vector<uint64_t> > getZeroMatrix(int dim1, int dim2){
  std::vector<std::vector<uint64_t> > temp(dim1, std::vector<uint64_t>(dim2, 0));
  return temp;
}



/**
 *
 * This functions writes paths to a file.
 *
 */
void StorePaths(std::vector<std::vector<uint64_t> > &paths, std::string file_path){

  std::ofstream data_file;
  data_file.open(file_path);

  for(unsigned int i = 0; i < paths.size(); i++){

    for(unsigned int j = 0; j < paths[0].size(); j++){

      if(j < paths[0].size() - 1)
        data_file << paths[i][j] << ",";
      else
        data_file << paths[i][j] << std::endl;

    }

  }

  data_file.close();

}



/**
 *
 * Reads paths from a file.
 *
 */

std::vector<std::vector<uint64_t> > readPaths(std::string file_path){

  std::ifstream data_file(file_path);
  std::vector<std::vector<uint64_t> > paths;

  std::string str;
  int i = 0;
  while (std::getline(data_file, str)){
    std::vector<uint64_t> path;
    std::string str2("");
    for(unsigned int j = 0; j < str.size(); j++){
      if(str[j] != ','){
        str2 += str[j];
      }
      else{
        uint64_t val = StringToInt(str2);
        path.push_back(val);
        str2 = "";
      }
    }
    uint64_t val = StringToInt(str2);
    path.push_back(val);
    str2 = "";
    paths.push_back(path);
    i++;
  }

  data_file.close();

  return paths;

}



/**
 *
 * This function parses paths into their 64 bit representations.
 *
 */
// [[Rcpp::export]]
void parsePaths(IntegerMatrix data, int nCases, int nControls, std::string file_path){

  int vlen = (int) ceil(nCases/64.0);
  int vlen2 = (int) ceil(nControls/64.0);

  std::vector<std::vector<uint64_t> > paths(data.nrow(), std::vector<uint64_t>(vlen + vlen2, 0));

  for(int i = 0; i < data.nrow(); i++){

    for(int j = 0; j < data.ncol(); j++){

      if(data(i,j) != 0){

        if(j < nCases)
          paths[i][j/64] |= one_64bit << j % 64;
        else
          paths[i][vlen + (j-nCases)/64] |= one_64bit << (j-nCases) % 64;

      }

    }

  }

  StorePaths(paths, file_path);

}



/**
 *
 * Convert CaseORControl from 0/1 values to 64 bit representations.
 *
 */

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

// [[Rcpp::export]]
List JoinIndices(IntegerVector srcuid, IntegerVector trguids2, List uids_CountLoc, IntegerVector joining_gene_sign,
  NumericMatrix ValueTable, int nCases, int nControls, int K,
  int iterations, IntegerMatrix CaseORControl, std::string method, int pathLength, int nthreads, std::string pos_path1,
  std::string neg_path1, std::string conflict_path1, std::string pos_path2, std::string neg_path2, std::string conflict_path2,
  std::string dest_path_pos, std::string dest_path_neg, std::string dest_path_conflict){

/*
  printf("################\n");
  printf("srcuid size : %d\n", srcuid.size());
  printf("trguid size : %d\n", trguids2.size());
  printf("uids cl size: %d\n", uids_CountLoc.size());
  printf("join gene s : %d\n", joining_gene_sign.size());
  printf("value table : %d x %d\n", ValueTable.rows(), ValueTable.cols());
  printf("caseorcon   : %d x %d\n", CaseORControl.rows(), CaseORControl.cols());
  printf("ncases      : %d\n", nCases);
  printf("ncontrols   : %d\n", nControls);
  printf("K           : %d\n", K);
  printf("iterations  : %d\n", iterations);
  printf("method      : %s\n", method.c_str());
  printf("pathlen     : %d\n", pathLength);
  printf("nthread     : %d\n", nthreads);
  printf("\n");

  printf("pos1 : %s\n", pos_path1.c_str());
  printf("neg1 : %s\n", neg_path1.c_str());
  printf("con1 : %s\n", conflict_path1.c_str());
  printf("pos2 : %s\n", pos_path2.c_str());
  printf("neg2 : %s\n", neg_path2.c_str());
  printf("con2 : %s\n", conflict_path2.c_str());
  printf("pos3 : %s\n", dest_path_pos.c_str());
  printf("neg3 : %s\n", dest_path_neg.c_str());
  printf("con3 : %s\n", dest_path_conflict.c_str());

  printf("################\n\n");
*/

  if(method == "method2") {
    return JoinIndicesMethod2(srcuid, trguids2, uids_CountLoc, joining_gene_sign,
      ValueTable, nCases, nControls, K,
      iterations, CaseORControl, pathLength, nthreads, pos_path1,
      neg_path1, conflict_path1, pos_path2, neg_path2, conflict_path2,
      dest_path_pos, dest_path_neg, dest_path_conflict);
  } else {
    stop("unknown method: " + method);
  }
}
