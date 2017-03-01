#include <fstream>
#include <Rcpp.h>
#include "gcre.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

const uint64_t one_64bit = 1;

/**
 *
 * Get the total number of paths to be processed.
 *
 */
int getTotalPaths(IntegerVector trguids, List uids_CountLoc){
  int total_paths = 0;
  for(int i = 0; i < trguids.size(); i++){
    int uid = trguids[i];
    string geneuid = to_string(uid);
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

vector<vector<uint64_t>> getZeroMatrix(int dim1, int dim2){
  vector<vector<uint64_t>> temp(dim1, vector<uint64_t>(dim2, 0));
  return temp;
}



/**
 *
 * This functions writes paths to a file.
 *
 */
void StorePaths(vector<vector<uint64_t>> &paths, string file_path){

  ofstream data_file;
  data_file.open(file_path);

  for(unsigned int i = 0; i < paths.size(); i++){

    for(unsigned int j = 0; j < paths[0].size(); j++){

      if(j < paths[0].size() - 1)
        data_file << paths[i][j] << ",";
      else
        data_file << paths[i][j] << endl;

    }

  }

  data_file.close();

}



/**
 *
 * Reads paths from a file.
 *
 */

vector<vector<uint64_t>> readPaths(string file_path){

  ifstream data_file(file_path);
  vector<vector<uint64_t>> paths;

  string str;
  int i = 0;
  while (getline(data_file, str)){
    vector<uint64_t> path;
    string str2("");
    for(unsigned int j = 0; j < str.size(); j++){
      if(str[j] != ','){
        str2 += str[j];
      }
      else{
        path.push_back((uint64_t) stoul(str2));
        str2 = "";
      }
    }
    path.push_back((uint64_t) stoul(str2));
    str2 = "";
    paths.push_back(path);
    i++;
  }

  data_file.close();

  return paths;

}


// conversions from R

vector<int> copy_rvector(IntegerVector r_vec)
{
  vector<int> vec(r_vec.size());
  for(int k = 0; k < r_vec.size(); k += 1)
    vec[k] = r_vec[k];
  return vec;
}

vec2dd copy_rmatrix(NumericMatrix matrix)
{
  vec2dd table(matrix.nrow(), vector<double>(matrix.ncol()));
  for(int r = 0; r < matrix.nrow(); r += 1){
    for(int c = 0; c < matrix.ncol(); c += 1)
      table[r][c] = matrix(r,c);
  }
  printf("copied value matrix: %d x %d\n", matrix.nrow(), matrix.ncol());
  return table;
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

  vector<vector<uint64_t>> paths(data.nrow(), vector<uint64_t>(vlen + vlen2, 0));

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

vector<vector<uint64_t>> parseCaseORControl(IntegerMatrix CaseORControl, int nCases, int nControls){


  int vlen = (int) ceil(nCases/64.0);
  int vlen2 = (int) ceil(nControls/64.0);

  vector<vector<uint64_t>> CaseORControl2(CaseORControl.nrow(), vector<uint64_t>(vlen + vlen2, 0));
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

/**
 *
 * Create list indexed by the names of gene uids. Each entry in the list
 * is a vector with the first entry being the count of paths starting
 * with the corresponding gene uids and the second entry be the location
 * of the first path in the relations data frame.
 *
 */
// [[Rcpp::export]]
List getMatchingList(IntegerVector uids, IntegerVector counts, IntegerVector location){
  List uids_2countsloc = List();
  for(unsigned int i = 0; i < uids.size(); i++){
    string uids_str = to_string(uids[i]);
    Rcpp::IntegerVector temp = IntegerVector(2);
    int count = counts[i];
    temp[0] = count;
    temp[1] = location[i];
    uids_2countsloc[uids_str] = temp;
  }
  return uids_2countsloc;
}

// [[Rcpp::export]]
List JoinIndices(IntegerVector r_src_uids, IntegerVector r_trg_uids, List uids_CountLoc, IntegerVector r_join_gene_signs,
  NumericMatrix r_value_table, int nCases, int nControls, int K,
  int iterations, IntegerMatrix CaseORControl, std::string method, int pathLength, int nthreads, std::string pos_path1,
  std::string neg_path1, std::string conflict_path1, std::string pos_path2, std::string neg_path2, std::string conflict_path2,
  std::string dest_path_pos, std::string dest_path_neg, std::string dest_path_conflict){

  vec2dd value_table = copy_rmatrix(r_value_table);
  vector<int> src_uids = copy_rvector(r_src_uids);
  vector<int> trg_uids = copy_rvector(r_trg_uids);
  vector<int> join_gene_signs = copy_rvector(r_join_gene_signs);

  

  printf("################\n");
  // printf("srcuid size : %d\n", srcuid.size());
  // printf("trguid size : %d\n", trguids2.size());
  // printf("uids cl size: %d\n", uids_CountLoc.size());
  // printf("join genes  : %d\n", joining_gene_sign.size());
  // printf("value table : %d x %d\n", ValueTable.rows(), ValueTable.cols());
  // printf("caseorcon   : %d x %d\n", CaseORControl.rows(), CaseORControl.cols());
  // printf("ncases      : %d\n", nCases);
  // printf("ncontrols   : %d\n", nControls);
  // printf("K           : %d\n", K);
  // printf("iterations  : %d\n", iterations);
  printf("method      : %s\n", method.c_str());
  printf("pathlen     : %d\n", pathLength);
  printf("nthread     : %d\n", nthreads);
  printf("\n");

  // printf("pos1 : %s\n", pos_path1.c_str());
  // printf("neg1 : %s\n", neg_path1.c_str());
  // printf("con1 : %s\n", conflict_path1.c_str());
  // printf("pos2 : %s\n", pos_path2.c_str());
  // printf("neg2 : %s\n", neg_path2.c_str());
  // printf("con2 : %s\n", conflict_path2.c_str());
  // printf("pos3 : %s\n", dest_path_pos.c_str());
  // printf("neg3 : %s\n", dest_path_neg.c_str());
  // printf("con3 : %s\n", dest_path_conflict.c_str());

  printf("################\n\n");


  if(method == "method2") {
    return join_method2(src_uids, trg_uids, uids_CountLoc, join_gene_signs,
      value_table, nCases, nControls, K,
      iterations, CaseORControl, pathLength, nthreads, pos_path1,
      neg_path1, conflict_path1, pos_path2, neg_path2, conflict_path2,
      dest_path_pos, dest_path_neg, dest_path_conflict, getTotalPaths(r_trg_uids, uids_CountLoc));
  } else if(method == "method2-old") {
    return JoinIndicesMethod2(r_src_uids, r_trg_uids, uids_CountLoc, r_join_gene_signs,
      r_value_table, nCases, nControls, K,
      iterations, CaseORControl, pathLength, nthreads, pos_path1,
      neg_path1, conflict_path1, pos_path2, neg_path2, conflict_path2,
      dest_path_pos, dest_path_neg, dest_path_conflict);
  } else {
    stop("unknown method: " + method);
  }
}
