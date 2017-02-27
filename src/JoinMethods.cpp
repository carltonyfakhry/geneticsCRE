#include "Utils.h"
#include <omp.h>
#include <fstream>
#include <emmintrin.h>
#include <nmmintrin.h>


// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]


using namespace Rcpp;



/**
 *
 * Brian Kernighan's algorithm for counting the number of bits
 * in a int64_t integer.
 *
 */
inline int bitCount(uint64_t n){

  int count = 0;

  while(n){

    n &= (n-one_64bit);
    count++;

  }

  return count;

}



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



/**
 *
 * Join the paths for Method 1.
 *
 */
// [[Rcpp::export]]
List JoinIndicesMethod1(IntegerVector srcuid, IntegerVector trguids2, List uids_CountLoc, IntegerVector joining_gene_sign,
                        NumericMatrix ValueTable, int nCases, int nControls, int K,
                        int iterations, IntegerMatrix CaseORControl, int method, int pathLength, int nthreads, std::string pos_path1,
                        std::string pos_path2, std::string dest_path_pos){

  int vlen = (int) ceil(nCases/64.0);
  int vlen2 = (int) ceil(nControls/64.0);

  // Copy ValueTable into a C++ 2D vector to be used in an openmp for loop
  std::vector<std::vector<double> > ValueTable2(ValueTable.nrow(), std::vector<double>(ValueTable.ncol()));
  for(int i = 0; i < ValueTable.nrow(); i++)
    for(int j = 0; j < ValueTable.ncol(); j++)
      ValueTable2[i][j] = ValueTable(i,j);

  // Store CaseORControl information into 64bit integers
  std::vector<std::vector<uint64_t> > CaseORControl2 = parseCaseORControl(CaseORControl, nCases, nControls);
  uint64_t CaseORControl22[iterations][vlen + vlen2] __attribute__ ((aligned (16)));
  for(int i = 0; i < iterations; i++)
    for(int j = 0; j < vlen + vlen2; j++)
      CaseORControl22[i][j] = CaseORControl2[i][j];


  // Copy the paths into 2D vectors
  std::vector<std::vector<uint64_t> > paths_pos1 = (pos_path1 == "") ? getZeroMatrix(trguids2.size(), vlen + vlen2) : (readPaths(pos_path1));
  std::vector<std::vector<uint64_t> > temp_paths_pos2;
  std::vector<std::vector<uint64_t> > &paths_pos2 = (pathLength == 5) ? paths_pos1 : (temp_paths_pos2 = readPaths(pos_path2));

  // A priority queue for the indices of the top K paths in the data
  IndicesScoresQueue indicesQ;
  for(int j = 0; j < K; j++) indicesQ.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));

  // If pathLength <=3 then we will store the joining of the paths
  int total_paths = (pathLength > 3) ? 0 : getTotalPaths(trguids2, uids_CountLoc);
  std::vector<std::vector<uint64_t> > joined_paths(total_paths, std::vector<uint64_t>(vlen + vlen2, 0));

  // nthreads can't be set to be more than the maximum number of possible threads
  int max_nthreads = omp_get_max_threads();
  nthreads = (nthreads > max_nthreads || nthreads == -1) ? max_nthreads : nthreads;
  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);

  // Create a priority queue which holds the top K scores per thread
  std::vector<IndicesScoresQueue> LocalIndicesQ;
  for(int k = 0; k < nthreads; k++){
    IndicesScoresQueue localindicesq;
    for(int j = 0; j < K; j++) localindicesq.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));
    LocalIndicesQ.push_back(localindicesq);
  }

  // Global 2D vector which holds the maximum scores per permutation for each corresponding thread
  std::vector<std::vector<double> > max_scores(nthreads, std::vector<double>(iterations, -INFINITY));

  int flipped_pivot_length = paths_pos2.size();
  int prev_src = -1;
  int total_srcs = 0;
  total_paths = 0;

  for(int i = 0; i < trguids2.size(); i++){

    // Check if source node has changed
    int src = srcuid[i];
    if(src != prev_src){
      prev_src = src;
      total_srcs++;
      if((total_srcs % 100) == 0){
        Rcout << "    " <<  total_srcs << " source nodes for paths of length " << pathLength << " and their permutations have been processed!" << std::endl;
        checkUserInterrupt();
      }
    }

    // Find the locations and the number of the paths matching with uid
    int uid = trguids2[i];
    std::string geneuid = IntToString(uid); // convert integer to string
    IntegerVector uid_count_loc = uids_CountLoc[geneuid];
    int count = uid_count_loc[0];
    if(count == 0) continue;
    int location = uid_count_loc[1];

    // Get the data of the first path
    std::vector<uint64_t> &path_pos1 = paths_pos1[i];

    // Copy the values from max_scores into a local container local_max_scores
    std::vector<std::vector<double> > local_max_scores(nthreads, std::vector<double>(iterations));
    for(int j = 0; j < nthreads; j++)
      for(int k = 0; k < iterations; k++)
        local_max_scores[j][k] = max_scores[j][k];

    #pragma omp parallel for shared(path_pos1,paths_pos2,flipped_pivot_length,location,count,ValueTable2,LocalIndicesQ,local_max_scores,CaseORControl2,vlen,vlen2,nCases,nControls,iterations) if(pathLength > 3)
    for(int j = location; j < (location + count); j++){

      int tid = omp_get_thread_num();
      IndicesScoresQueue &tid_localindicesq = LocalIndicesQ[tid];
      std::vector<uint64_t> joined_pos(path_pos1.size());

      // Get the data for path2 which is to be joined with path1
      std::vector<uint64_t> &path_pos2 = paths_pos2[j];

      // Join the paths
      for(int k = 0; k < vlen + vlen2; k++) joined_pos[k] = path_pos1[k] | path_pos2[k];

      // Store the joined paths if pathLength <= 3
      if(pathLength <= 3){
        joined_paths[total_paths] = joined_pos;
        total_paths++;
      }

      int cases = 0;
      int controls = 0;
      for(int k = 0; k < vlen; k++) cases += __builtin_popcountll(joined_pos[k]);
      for(int k = vlen; k < vlen + vlen2; k++) controls += __builtin_popcountll(joined_pos[k]);

      double score = ValueTable2[cases][controls];
      double flipped_score = ValueTable2[controls][cases];

      if(score > tid_localindicesq.top().first){
        tid_localindicesq.push(std::pair<double, std::pair<int, int> >(score, std::pair<int,int>(i, j)));
        tid_localindicesq.pop();
      }

      if(flipped_score > tid_localindicesq.top().first){
        tid_localindicesq.push(std::pair<double, std::pair<int, int> >(flipped_score, std::pair<int,int>(i, j + flipped_pivot_length)));
        tid_localindicesq.pop();
      }

      std::vector<double> &tid_max_scores = local_max_scores[tid];
      uint64_t joined_pos2[vlen + vlen2] __attribute__ ((aligned (16)));
      for(int k = 0; k < vlen + vlen2; k++) joined_pos2[k] = joined_pos[k];

      for(int m = 0; m < iterations; m++){

        double &max_score = tid_max_scores[m];
        double max_score2 = tid_max_scores[m];
        double perm_score = 0;
        double perm_flipped_score = 0;

        std::vector<uint64_t> &caseorcontrol = CaseORControl2[m];
        uint64_t *p_joined_pos = joined_pos2;
        uint64_t *p_caseorcontrol2 = CaseORControl22[m];

        int cases_m = 0;

        int k = 0;
        for(; k < ROUND_DOWN(vlen,2); k+=2){
          __m128i* ptr1 = (__m128i*) (p_caseorcontrol2);
          __m128i* ptr2 = (__m128i*) (p_joined_pos);
          __m128i val1_4 = _mm_loadu_si128(ptr1);
          __m128i val2_4 = _mm_loadu_si128(ptr2);
          uint64_t perm_joins[] __attribute__ ((aligned (16))) = {0,0};
          _mm_store_si128((__m128i*) perm_joins, _mm_and_si128(val1_4, val2_4));
          cases_m += __builtin_popcountll(perm_joins[0]);
          cases_m += __builtin_popcountll(perm_joins[1]);
          p_joined_pos+=2;
          p_caseorcontrol2+=2;
        }
        if(k < vlen - 1){
          p_joined_pos++;
          p_caseorcontrol2++;
          __m128i* ptr1 = (__m128i*) (p_caseorcontrol2);
          __m128i* ptr2 = (__m128i*) (p_joined_pos);
          __m128i val1_4 = _mm_loadu_si128(ptr1);
          __m128i val2_4 = _mm_loadu_si128(ptr2);
          uint64_t perm_joins[] __attribute__ ((aligned (16))) = {0,0};
          _mm_store_si128((__m128i*) perm_joins, _mm_and_si128(val1_4, val2_4));
          cases_m += __builtin_popcountll(perm_joins[0]);
          cases_m += __builtin_popcountll(perm_joins[1]);
          p_caseorcontrol2++;
          p_joined_pos++;
        }

        // for(int k = 0; k < vlen; k++){
        //   uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
        //   cases_m += __builtin_popcountll(permuted_path_k);
        // }

        int controls_m = 0;

        k = vlen;
        for(; k < ROUND_DOWN(vlen + vlen2,2); k+=2){
          __m128i* ptr1 = (__m128i*) (p_caseorcontrol2);
          __m128i* ptr2 = (__m128i*) (p_joined_pos);
          __m128i val1_4 = _mm_loadu_si128(ptr1);
          __m128i val2_4 = _mm_loadu_si128(ptr2);
          uint64_t perm_joins[] __attribute__ ((aligned (16))) = {0,0};
          _mm_storeu_si128((__m128i*) perm_joins, _mm_and_si128(val1_4, val2_4));
          controls_m += __builtin_popcountll(perm_joins[0]);
          controls_m += __builtin_popcountll(perm_joins[1]);
          p_joined_pos+=2;
          p_caseorcontrol2+=2;
        }
        if(k < vlen + vlen2 - 1){
          p_joined_pos++;
          p_caseorcontrol2++;
          __m128i* ptr1 = (__m128i*) (p_caseorcontrol2);
          __m128i* ptr2 = (__m128i*) (p_joined_pos);
          __m128i val1_4 = _mm_loadu_si128(ptr1);
          __m128i val2_4 = _mm_loadu_si128(ptr2);
          uint64_t perm_joins[] __attribute__ ((aligned (16))) = {0,0};
          _mm_storeu_si128((__m128i*) perm_joins, _mm_and_si128(val1_4, val2_4));
          controls_m += __builtin_popcountll(perm_joins[0]);
          controls_m += __builtin_popcountll(perm_joins[1]);
        }

        // for(int k = vlen; k < vlen + vlen2; k++){
        //   uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
        //   controls_m += __builtin_popcountll(permuted_path_k);
        // }

        int new_cases_m = cases_m + (controls - controls_m);
        int new_controls_m = controls_m + (cases - cases_m);

        perm_score = ValueTable2[new_cases_m][new_controls_m];
        if(perm_score > max_score2){
          max_score = perm_score;
          max_score2 = perm_score;
        }
        perm_flipped_score = ValueTable2[new_controls_m][new_cases_m];
        if(perm_flipped_score > max_score2){
          max_score = perm_flipped_score;
        }

      }

    }

    // Update the values in max_scores given the values of local_max_scores
    for(int j = 0; j < nthreads; j++){
      for(int k = 0; k < iterations; k++)
        if(max_scores[j][k] < local_max_scores[j][k]) max_scores[j][k] = local_max_scores[j][k];
    }

  }

  NumericVector permutedscores(iterations);
  for(int i = 0; i < nthreads; i++){
    for(int j = 0; j < iterations; j++){
      if(i == 0){
        permutedscores[j] = max_scores[i][j];
      }else if(max_scores[i][j] > permutedscores[j]){
        permutedscores[j] = max_scores[i][j];
      }
    }
  }

  for(unsigned int j = 0; j < LocalIndicesQ.size(); j++){
    IndicesScoresQueue &localqueue = LocalIndicesQ[j];
    while(localqueue.size() != 0) {
      indicesQ.push(localqueue.top());
      localqueue.pop();
      indicesQ.pop();
    }
  }

  // Create a vector for the indices of the k paths with the highest scores
  IntegerMatrix ids(indicesQ.size(),2);
  NumericVector scores(indicesQ.size());
  int j = 0;
  while(indicesQ.size() != 0){
    ScoreIndices temp = indicesQ.top();
    scores[j] = temp.first;
    ids(j,0) = temp.second.first + 1; // Add 1 because it will be used as an index in the R code
    ids(j,1) = temp.second.second + 1;
    indicesQ.pop();
    j++;
  }

  if(dest_path_pos != "")
    StorePaths(joined_paths, dest_path_pos);

  return List::create(Named("scores") = scores,
                      Named("ids") = ids,
                      Named("TestScores") = permutedscores);

}



/**
 *
 * Join the paths for Method 2.
 *
 */
// [[Rcpp::export]]
List JoinIndicesMethod2(IntegerVector srcuid, IntegerVector trguids2, List uids_CountLoc, IntegerVector joining_gene_sign,
                        NumericMatrix ValueTable, int nCases, int nControls, int K,
                        int iterations, IntegerMatrix CaseORControl, int method, int pathLength, int nthreads, std::string pos_path1,
                        std::string neg_path1, std::string conflict_path1, std::string pos_path2, std::string neg_path2, std::string conflict_path2,
                        std::string dest_path_pos, std::string dest_path_neg, std::string dest_path_conflict){

  int vlen = (int) ceil(nCases/64.0);
  int vlen2 = (int) ceil(nControls/64.0);

  // Copy ValueTable into a C++ 2D vector to be used in an openmp for loop
  std::vector<std::vector<double> > ValueTable2(ValueTable.nrow(), std::vector<double>(ValueTable.ncol()));
  for(int i = 0; i < ValueTable.nrow(); i++)
    for(int j = 0; j < ValueTable.ncol(); j++)
      ValueTable2[i][j] = ValueTable(i,j);

  // If pathLength <=3 then we will store the joining of the different paths
  int total_paths = (pathLength > 3) ? 0 : getTotalPaths(trguids2, uids_CountLoc);
  std::vector<std::vector<uint64_t> > joined_paths_pos(total_paths, std::vector<uint64_t>(vlen + vlen2, 0));
  std::vector<std::vector<uint64_t> > joined_paths_neg(total_paths, std::vector<uint64_t>(vlen + vlen2, 0));
  std::vector<std::vector<uint64_t> > joined_paths_conflict(total_paths, std::vector<uint64_t>(vlen + vlen2, 0));

  // Store CaseORControl information into 64bit integers
  std::vector<std::vector<uint64_t> > CaseORControl2 = parseCaseORControl(CaseORControl, nCases, nControls);

  // Copy the paths into 2D vectors
  std::vector<std::vector<uint64_t> > temp_paths_pos2;
  std::vector<std::vector<uint64_t> > temp_paths_neg2;
  std::vector<std::vector<uint64_t> > temp_paths_neg22;
  std::vector<std::vector<uint64_t> > temp_paths_conflict2;
  std::vector<std::vector<uint64_t> > temp_paths_conflict22;
  int total_srcuidsRels2 = getTotalCountsCountLoc(uids_CountLoc);
  std::vector<std::vector<uint64_t> > paths_pos1 = (pos_path1 == "") ? getZeroMatrix(trguids2.size(), vlen + vlen2) : readPaths(pos_path1); // Only empty paths for paths of length 1
  std::vector<std::vector<uint64_t> > &paths_pos2 = (pathLength == 5) ? paths_pos1 : temp_paths_pos2 = readPaths(pos_path2);
  std::vector<std::vector<uint64_t> > paths_neg1 = (neg_path1 == "") ? getZeroMatrix(trguids2.size(), vlen + vlen2) : readPaths(neg_path1); // Only empty paths for paths of length 1
  std::vector<std::vector<uint64_t> > &paths_neg2 = (neg_path2 == "") ? temp_paths_neg2 = getZeroMatrix(total_srcuidsRels2, vlen + vlen2) : temp_paths_neg2 = (pathLength == 5) ? paths_neg1 : temp_paths_neg22 = readPaths(neg_path2); // Only empty paths for paths of length 1
  std::vector<std::vector<uint64_t> > paths_conflict1 = (conflict_path1 == "") ? getZeroMatrix(trguids2.size(), vlen + vlen2) : readPaths(conflict_path1); // Only empty paths for paths of length 1
  std::vector<std::vector<uint64_t> > &paths_conflict2 = (conflict_path2 == "") ? temp_paths_conflict2 = getZeroMatrix(total_srcuidsRels2, vlen + vlen2) : temp_paths_conflict2 = (pathLength == 5) ? paths_conflict1 : temp_paths_conflict22 = readPaths(conflict_path2); // Only empty paths for paths of length 1

  // A priority queue for the indices of the top K paths in the data
  IndicesScoresQueue indicesQ;
  for(int j = 0; j < K; j++) indicesQ.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));

  // nthreads can't be set to be more than the maximum number of possible threads
  int max_nthreads = omp_get_max_threads();
  nthreads = (nthreads > max_nthreads || nthreads == -1) ? max_nthreads : nthreads;
  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);

  // Create a priority queue which holds the top K scores per thread
  std::vector<IndicesScoresQueue> LocalIndicesQ;
  for(int k = 0; k < nthreads; k++){
    IndicesScoresQueue localindicesq;
    for(int j = 0; j < K; j++) localindicesq.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));
    LocalIndicesQ.push_back(localindicesq);
  }

  // Global 2D vector which holds the maximum scores per permutation for each corresponding thread
  std::vector<std::vector<double> > max_scores(nthreads, std::vector<double>(iterations, -INFINITY));

  int flipped_pivot_length = paths_pos2.size();
  int prev_src = -1;
  int total_srcs = 0;
  total_paths = 0;

  for(int i = 0; i < trguids2.size(); i++){

    // Check if source node has changed
    int src = srcuid[i];
    if(src != prev_src){
      prev_src = src;
      total_srcs++;
      if((total_srcs % 100) == 0){
        Rcout << "    " <<  total_srcs << " source nodes for paths of length " << pathLength << " and their permutations have been processed!" << std::endl;
        checkUserInterrupt();
      }
    }

    // Find the locations and the number of the paths matching with uid
    int uid = trguids2[i];
    std::string geneuid = IntToString(uid); // convert integer to string
    IntegerVector uid_count_loc = uids_CountLoc[geneuid];
    int count = uid_count_loc[0];
    if(count == 0) continue;
    int location = uid_count_loc[1];
    int sign;
    if(pathLength > 3) sign = joining_gene_sign[i];

    // Get the data of the first path
    std::vector<uint64_t> &path_pos1 = paths_pos1[i];
    std::vector<uint64_t> &path_neg1 = paths_neg1[i];
    std::vector<uint64_t> &path_conflict1 = paths_conflict1[i];

    // Copy the values from max_scores into a local container local_max_scores to increase performance
    std::vector<std::vector<double> > local_max_scores(nthreads, std::vector<double>(iterations));
    for(int j = 0; j < nthreads; j++)
      for(int k = 0; k < iterations; k++)
        local_max_scores[j][k] = max_scores[j][k];

    #pragma omp parallel for shared(path_pos1,path_neg1,path_conflict1,paths_pos2,paths_neg2,paths_conflict2,flipped_pivot_length,location,count,ValueTable2,LocalIndicesQ,local_max_scores,CaseORControl2,vlen,vlen2,nCases,nControls,iterations) if(pathLength > 3)
    for(int j = location; j < (location + count); j++){

      int tid = omp_get_thread_num();
      IndicesScoresQueue &tid_localindicesq = LocalIndicesQ[tid];

      if(pathLength < 3) sign = joining_gene_sign[j];

      if(pathLength == 3){
        int sign2 = joining_gene_sign[i];
        int sign3 = joining_gene_sign[j];
        sign = ((sign2 == -1 && sign3 == 1) || (sign2 == 1 && sign3 == -1)) ? -1 : 1;
      }

      std::vector<uint64_t> joined_pos(path_pos1.size());
      std::vector<uint64_t> joined_neg(path_pos1.size());
      std::vector<uint64_t> joined_conflict(path_pos1.size());

      std::vector<uint64_t> &path_pos2 = (sign == 1) ? paths_pos2[j] : paths_neg2[j];
      std::vector<uint64_t> &path_neg2 = (sign == 1) ? paths_neg2[j] : paths_pos2[j];
      std::vector<uint64_t> &path_conflict2 = paths_conflict2[j];

      for(int k = 0; k < vlen + vlen2; k++){
        uint64_t temp_pos = path_pos1[k] | path_pos2[k];
        uint64_t temp_neg = path_neg1[k] | path_neg2[k];
        uint64_t temp_conflict = (path_conflict1[k] | path_conflict2[k]) | (temp_pos & temp_neg);
        joined_conflict[k] = temp_conflict;
        joined_pos[k] = temp_pos ^ (temp_conflict & temp_pos);
        joined_neg[k] = temp_neg ^ (temp_conflict & temp_neg);
      }

      if(pathLength <= 3){
        joined_paths_pos[total_paths] = joined_pos;
        joined_paths_neg[total_paths] = joined_neg;
        joined_paths_conflict[total_paths] = joined_conflict;
        total_paths++;
      }

      int cases_pos = 0;
      int cases_neg = 0;
      int controls_pos = 0;
      int controls_neg = 0;
      for(int k = 0; k < vlen; k++){
        cases_pos += __builtin_popcountll(joined_pos[k]);
        controls_neg += __builtin_popcountll(joined_neg[k]);
      }
      for(int k = vlen; k < vlen + vlen2; k++){
        controls_pos += __builtin_popcountll(joined_pos[k]);
        cases_neg += __builtin_popcountll(joined_neg[k]);
      }

      int cases = cases_pos + cases_neg;
      int controls = controls_pos + controls_neg;

      double score = ValueTable2[cases][controls];
      double flipped_score = ValueTable2[controls][cases];

      if(score > tid_localindicesq.top().first){
        tid_localindicesq.push(std::pair<double, std::pair<int, int> >(score, std::pair<int,int>(i, j)));
        tid_localindicesq.pop();
      }

      if(flipped_score > tid_localindicesq.top().first){
        tid_localindicesq.push(std::pair<double, std::pair<int, int> >(flipped_score, std::pair<int,int>(i, j + flipped_pivot_length)));
        tid_localindicesq.pop();
      }

      std::vector<double> &tid_max_scores = local_max_scores[tid];

      for(int m = 0; m < iterations; m++){

        double &max_score = tid_max_scores[m];
        double max_score2 = tid_max_scores[m];
        double perm_score = 0;
        double perm_flipped_score = 0;

        std::vector<uint64_t> &caseorcontrol = CaseORControl2[m];
        int cases_pos_m = 0;
        int controls_pos_m = 0;
        int cases_neg_m = 0;
        int controls_neg_m = 0;
        for(int k = 0; k < vlen; k++){
          uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
          cases_pos_m += __builtin_popcountll(permuted_path_k);
          permuted_path_k = joined_neg[k] & caseorcontrol[k];
          controls_neg_m += __builtin_popcountll(permuted_path_k);
        }
        for(int k = vlen; k < vlen + vlen2; k++){
          uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
          controls_pos_m += __builtin_popcountll(permuted_path_k);
          permuted_path_k = joined_neg[k] & caseorcontrol[k];
          cases_neg_m += __builtin_popcountll(permuted_path_k);
        }

        int new_cases = cases_pos_m + cases_neg_m + (controls_pos - controls_pos_m) + (controls_neg - controls_neg_m);
        int new_controls = controls_pos_m + controls_neg_m + (cases_pos - cases_pos_m) + (cases_neg - cases_neg_m);

        perm_score = ValueTable2[new_cases][new_controls];
        if(perm_score > max_score2){
          max_score = perm_score;
          max_score2 = perm_score;
        }
        perm_flipped_score = ValueTable2[new_controls][new_cases];
        if(perm_flipped_score > max_score2){
          max_score = perm_flipped_score;
        }

      }

    }

    // Update the values in max_scores given the values of local_max_scores
    for(int j = 0; j < nthreads; j++){
      for(int k = 0; k < iterations; k++)
        if(max_scores[j][k] < local_max_scores[j][k]) max_scores[j][k] = local_max_scores[j][k];
    }

  }

  NumericVector permutedscores(iterations);
  for(int i = 0; i < nthreads; i++){
    for(int j = 0; j < iterations; j++){
      if(i == 0){
        permutedscores[j] = max_scores[i][j];
      }else if(max_scores[i][j] > permutedscores[j]){
        permutedscores[j] = max_scores[i][j];
      }
    }
  }

  for(unsigned int j = 0; j < LocalIndicesQ.size(); j++){
    IndicesScoresQueue &localqueue = LocalIndicesQ[j];
    while(localqueue.size() != 0) {
      indicesQ.push(localqueue.top());
      localqueue.pop();
      indicesQ.pop();
    }
  }

  // Create a vector for the indices of the k paths with the highest scores
  IntegerMatrix ids(indicesQ.size(),2);
  NumericVector scores(indicesQ.size());
  int j = 0;
  while(indicesQ.size() != 0){
    ScoreIndices temp = indicesQ.top();
    scores[j] = temp.first;
    ids(j,0) = temp.second.first + 1; // Add 1 because it will be used as an index in the R code
    ids(j,1) = temp.second.second + 1;
    indicesQ.pop();
    j++;
  }

  if(dest_path_pos != ""){
    StorePaths(joined_paths_pos, dest_path_pos);
    StorePaths(joined_paths_neg, dest_path_neg);
    StorePaths(joined_paths_conflict, dest_path_conflict);
  }

  return List::create(Named("scores") = scores,
                      Named("ids") = ids,
                      Named("TestScores") = permutedscores);


}

