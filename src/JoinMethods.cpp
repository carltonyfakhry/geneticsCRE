#include "Utils.h"
#include <omp.h>
#include <fstream>

// [[Rcpp::plugins(openmp)]]
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
 * This function returns a 2D vector filled with zeros.
 *
 */

std::vector<std::vector<uint64_t> > getZeroMatrix(int dim1, int dim2){
  std::vector<std::vector<uint64_t> > temp(dim1, std::vector<uint64_t>(dim2, 0));
  return temp;
}



/**
 *
 * Write paths into a file.
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
 * Parse paths to their 64 bit representations.
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
                        NumericMatrix ValueTable, NumericMatrix queues4init, int nCases, int nControls, int K,
                        int iterations, IntegerMatrix CaseORControl, int method, int pathLength, int nthreads, std::string pos_path1,
                        std::string pos_path2, std::string dest_path_pos){

  int vlen = (int) ceil(nCases/64.0);
  int vlen2 = (int) ceil(nControls/64.0);

  std::vector<std::vector<double> > ValueTable2(ValueTable.nrow(), std::vector<double>(ValueTable.ncol()));
  for(int i = 0; i < ValueTable.nrow(); i++)
    for(int j = 0; j < ValueTable.ncol(); j++)
      ValueTable2[i][j] = ValueTable(i,j);

  // Store CaseORControl information into 64bit integers
  std::vector<std::vector<uint64_t> > CaseORControl2 = parseCaseORControl(CaseORControl, nCases, nControls);

  // Copy the paths into 2D vectors
  std::vector<std::vector<uint64_t> > paths_pos1 = (pos_path1 == "") ? getZeroMatrix(trguids2.size(), vlen + vlen2) : (readPaths(pos_path1));
  std::vector<std::vector<uint64_t> > paths_pos2 = (pathLength == 5) ? paths_pos1 : (readPaths(pos_path2));

  // std::cout << "parsedpaths" << std::endl;

  // A priority queue for the indices of the top K paths
  IndicesScoresQueue indicesQ;
  for(int j = 0; j < K; j++) indicesQ.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));

  // If pathLength <=3 then we will store the joining of the different paths
  int total_paths = (pathLength > 3) ? 0 : getTotalPaths(trguids2, uids_CountLoc);
  std::vector<std::vector<uint64_t> > joined_paths(total_paths, std::vector<uint64_t>(vlen + vlen2, 0));

  int flipped_pivot_length = pos_path2.size();
  int prev_src = -1;
  int total_srcs = 0;
  total_paths = 0;

  // int max_available_threads = omp_get_max_threads();
  // if(nthreads == -1 || nthreads > max_available_threads){
  //   // int max_available_cores = omp_get_num_procs();
  //   omp_set_num_threads(max_available_threads);
  //   nthreads = max_available_cores;
  // }
  nthreads = omp_get_max_threads();
  // omp_set_dynamic(0);
  omp_set_num_threads(nthreads);
  // std::cout << nthreads << std::endl;

  // Create a priority queue which holds the top K scores per thread
  std::vector<IndicesScoresQueue> LocalIndicesQ;
  for(int k = 0; k < nthreads; k++){
    IndicesScoresQueue localindicesq;
    for(int j = 0; j < K; j++) localindicesq.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));
    LocalIndicesQ.push_back(localindicesq);
  }

  // Global container which holds the maximum scores per permutation for each corresponding thread
  std::vector<std::vector<double> > max_scores(nthreads, std::vector<double>(iterations, -INFINITY));

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

    // Copy the values from max_scores into a local container local_max_scores to increase performance
    std::vector<std::vector<double> > local_max_scores(nthreads, std::vector<double>(iterations));
    for(int j = 0; j < nthreads; j++)
      for(int k = 0; k < iterations; k++)
        local_max_scores[j][k] = max_scores[j][k];

    #pragma omp parallel for schedule(dynamic,1) shared(path_pos1,paths_pos2,flipped_pivot_length,location,count,ValueTable2,LocalIndicesQ,local_max_scores,CaseORControl2,vlen,vlen2,nCases,nControls,iterations) if(pathLength > 3)
    for(int j = location; j < (location + count); j++){

      int tid = omp_get_thread_num();

      IndicesScoresQueue &tid_localindicesq = LocalIndicesQ[tid];

      std::vector<uint64_t> joined_pos(path_pos1.size());

      // Get the data for path2 which is to be joined with path1
      std::vector<uint64_t> &path_pos2 = paths_pos2[j];
      for(int k = 0; k < vlen + vlen2; k++) joined_pos[k] = path_pos1[k] | path_pos2[k];

      if(pathLength <= 3){
        joined_paths[total_paths] = joined_pos;
        total_paths++;
      }

      int cases = 0;
      int controls = 0;
      for(int k = 0; k < vlen; k++) cases += __builtin_popcount(joined_pos[k]);
      for(int k = vlen; k < vlen + vlen2; k++) controls += __builtin_popcount(joined_pos[k]);

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
        int cases_m = 0;
        int controls_m = 0;

        for(int k = 0; k < vlen; k++){
          uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
          cases_m += __builtin_popcount(permuted_path_k);
        }

        for(int k = vlen; k < vlen + vlen2; k++){
          uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
          controls_m += __builtin_popcount(permuted_path_k);
        }

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

  return List::create();

}



/**
 *
 * Join the paths for Method 2.
 *
 */
// [[Rcpp::export]]
List JoinIndicesMethod2(IntegerVector srcuid, IntegerVector trguids2, List uids_CountLoc, IntegerVector joining_gene_sign,
                        NumericMatrix ValueTable, NumericMatrix queues4init, int nCases, int nControls, int K,
                        int iterations, IntegerMatrix CaseORControl, int method, int pathLength, int nthreads, std::string pos_path1,
                        std::string neg_path1, std::string conflict_path1, std::string pos_path2, std::string neg_path2, std::string conflict_path2,
                        std::string dest_path_pos, std::string dest_path_neg, std::string dest_path_conflict){

  int vlen = (int) ceil(nCases/64.0);
  int vlen2 = (int) ceil(nControls/64.0);

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
  std::vector<std::vector<uint64_t> > paths_pos1 = (pos_path1 == "") ? getZeroMatrix(trguids2.size(), vlen + vlen2) : (readPaths(pos_path1)); // Only empty path for paths of length 1
  std::vector<std::vector<uint64_t> > paths_pos2 = (pos_path2 == "") ? getZeroMatrix(total_paths, vlen + vlen2) : (readPaths(pos_path2));
  std::vector<std::vector<uint64_t> > paths_neg1 = (neg_path1 == "") ? getZeroMatrix(trguids2.size(), vlen + vlen2) : (readPaths(neg_path1)); // Only empty path for paths of length 1
  std::vector<std::vector<uint64_t> > paths_neg2 = (neg_path2 == "") ? getZeroMatrix(total_paths, vlen + vlen2) : (readPaths(neg_path2)); // Only empty path for paths of length 1
  std::vector<std::vector<uint64_t> > paths_conflict1 = (conflict_path1 == "") ? getZeroMatrix(trguids2.size(), vlen + vlen2) : (readPaths(conflict_path1)); // Only empty path for paths of length 1
  std::vector<std::vector<uint64_t> > paths_conflict2 = (conflict_path2 == "") ? getZeroMatrix(total_paths, vlen + vlen2) : (readPaths(conflict_path2)); // Only empty path for paths of length 1

  // A priority queue for the indices of the top K paths
  IndicesScoresQueue indicesQ;
  for(int j = 0; j < K; j++) indicesQ.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));

  int flipped_pivot_length = pos_path2.size();
  int prev_src = -1;
  int total_srcs = 0;
  total_paths = 0;

  // if(nthreads == -1){
  //   int max_available_cores = omp_get_num_procs();
  //   omp_set_num_threads(max_available_cores);
  //   nthreads = max_available_cores;
  // }

  nthreads = omp_get_max_threads();
  // omp_set_dynamic(0);
  omp_set_num_threads(nthreads);

  // Create a priority queue which holds the top K scores per thread
  std::vector<IndicesScoresQueue> LocalIndicesQ;
  for(int k = 0; k < nthreads; k++){
    IndicesScoresQueue localindicesq;
    for(int j = 0; j < K; j++) localindicesq.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));
    LocalIndicesQ.push_back(localindicesq);
  }

  // Global container which holds the maximum scores per permutation for each corresponding thread
  std::vector<std::vector<double> > max_scores(nthreads, std::vector<double>(iterations, -INFINITY));

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

    #pragma omp parallel for schedule(dynamic,1) shared(path_pos1,path_neg1,path_conflict1,paths_pos2,paths_neg2,paths_conflict2,flipped_pivot_length,location,count,ValueTable2,LocalIndicesQ,local_max_scores,CaseORControl2,vlen,vlen2,nCases,nControls,iterations) if(pathLength > 3)
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
        cases_pos += __builtin_popcount(joined_pos[k]);
        controls_neg += __builtin_popcount(joined_neg[k]);
      }
      for(int k = vlen; k < vlen + vlen2; k++){
        controls_pos += __builtin_popcount(joined_pos[k]);
        cases_neg += __builtin_popcount(joined_neg[k]);
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
          cases_pos_m += __builtin_popcount(permuted_path_k);
          permuted_path_k = joined_neg[k] & caseorcontrol[k];
          controls_neg_m += __builtin_popcount(permuted_path_k);
        }
        for(int k = vlen; k < vlen + vlen2; k++){
          uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
          controls_pos_m += __builtin_popcount(permuted_path_k);
          permuted_path_k = joined_neg[k] & caseorcontrol[k];
          cases_neg_m += __builtin_popcount(permuted_path_k);
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



// /**
//  *
//  * Initialize the priority queues.
//  *
//  */
//
// std::vector<std::priority_queue<double, std::vector<double>, std::greater<double> > > InitQueues(int iterations, int K, NumericMatrix queues4init){
// // std::vector<std::priority_queue<double, std::vector<double>, std::greater<double> > > InitQueues(int iterations, int K){
//
//   std::vector<std::priority_queue<double, std::vector<double>, std::greater<double> > > topKScoresQs;
//
//   for(int i = 0; i < iterations; i++){
//
//     std::priority_queue<double, std::vector<double>, std::greater<double> > tempQ;
//     if(queues4init.nrow() != 0){
//       for(int j = 0; j < K; j++) tempQ.push(-INFINITY);
//     }else{
//       for(int j = 0; j < K; j++) tempQ.push(queues4init(i+1,j));
//     }
//     topKScoresQs.push_back(tempQ);
//
//   }
//
//   return topKScoresQs;
//
// }


// std::vector<std::vector<uint64_t> > parsePaths(IntegerVector cols_pos, IntegerVector counts_pos,
//                                                int nCases, int nControls){
//
//   int vlen = (int) ceil(nCases/64.0);
//   int vlen2 = (int) ceil(nControls/64.0);
//   int total_count_pos = 0;
//
//   std::vector<std::vector<uint64_t> > paths(counts_pos.size(), std::vector<uint64_t>(vlen + vlen2, 0));
//
//   for(unsigned i = 0; i < paths.size(); i++){
//
//     std::vector<uint64_t> path = paths[i];
//     int count_pos = counts_pos[i];
//     total_count_pos += count_pos;
//
//     for(int j = total_count_pos - count_pos; j < total_count_pos; j++){
//
//       int col = cols_pos[j];
//       int which_double_index = (col < nCases) ? col/64 : vlen + (col-nCases)/64;
//       path[which_double_index] |= one_64bit << col % 64;
//
//     }
//
//     paths[i] = path;
//
//   }
//
//   return paths;
//
// }



// // [[Rcpp::export]]
// NumericMatrix AprioriTopKScores(IntegerVector counts_neg1, IntegerVector cols_neg1,
//                          IntegerVector counts_pos1, IntegerVector cols_pos1,
//                          IntegerVector counts_conflict1, IntegerVector cols_conflict1,
//                          IntegerVector counts_neg2, IntegerVector cols_neg2, IntegerVector locations_neg2,
//                          IntegerVector counts_pos2, IntegerVector cols_pos2, IntegerVector locations_pos2,
//                          IntegerVector counts_conflict2, IntegerVector cols_conflict2, IntegerVector locations_conflict2,
//                          IntegerVector inds1, IntegerVector trguids2, List uids_CountLoc, IntegerVector joining_gene_sign,
//                          NumericMatrix ValueTable, int K, int nCases, int nControls, IntegerMatrix CaseORControl,
//                          int method, int iterations){
//
//
//   int vlen = (int) ceil(nCases/64.0);
//   int vlen2 = (int) ceil(nControls/64.0);
//
//   NumericMatrix Scores4QueueInitiation(iterations + 1, K);
//   std::vector<double> min_scores(iterations + 1);
//   for(int i = 0; i < iterations + 1; i++) min_scores[i] = -INFINITY;
//
//   // Store CaseORControl information into 64bit integers
//   std::vector<std::vector<uint64_t> > CaseORControl2 = parseCaseORControl(CaseORControl, nCases, nControls);
//
//   // Init iterations + 1 queues for the top K scores of the paths of the initial and permuted data
//   NumericMatrix init4Queues(0,0);
//   VectorPriorityQueues topKScoresQs = InitQueues(iterations + 1, K, init4Queues);
//   // VectorPriorityQueues topKScoresQs = InitQueues(iterations + 1, K);
//
//   // Copy data into bitsets for increased efficiency
//   std::vector<std::vector<uint64_t> > paths_pos1 = parsePaths(cols_pos1, counts_pos1, nCases, nControls);
//   std::vector<std::vector<uint64_t> > paths_pos2 = parsePaths(cols_pos2, counts_pos2, nCases, nControls);
//   std::vector<std::vector<uint64_t> > paths_neg1 = parsePaths(cols_neg1, counts_neg1, nCases, nControls);
//   std::vector<std::vector<uint64_t> > paths_neg2 = parsePaths(cols_neg2, counts_neg2, nCases, nControls);
//   std::vector<std::vector<uint64_t> > paths_conflict1 = parsePaths(cols_conflict1, counts_conflict1, nCases, nControls);
//   std::vector<std::vector<uint64_t> > paths_conflict2 = parsePaths(cols_conflict2, counts_conflict2, nCases, nControls);
//
//   for(int i = 0; i < inds1.size(); i++){
//
//     std::vector<uint64_t> path_pos1 = paths_pos1[inds1[i]];
//     std::vector<uint64_t> path_neg1 = paths_neg1[inds1[i]];
//     std::vector<uint64_t> path_conflict1 = paths_conflict1[inds1[i]];
//
//     // Find the locations and the number of the paths matching with uid
//     int uid = trguids2[i];
//     std::string geneuid = IntToString(uid); // convert integer to string
//     IntegerVector uid_count_loc = uids_CountLoc[geneuid];
//     int count = uid_count_loc[0];
//     if(count == 0) continue;
//     int location = uid_count_loc[1];
//
//     for(int j = location; j < location + count; j++){
//
//       std::vector<uint64_t> path_pos2 = paths_pos2[j];
//       std::vector<uint64_t> path_neg2 = paths_neg2[j];
//       std::vector<uint64_t> path_conflict2 = paths_conflict2[j];
//
//       std::vector<uint64_t> joined_pos(vlen + vlen2);
//       std::vector<uint64_t> joined_neg(vlen + vlen2);
//       std::vector<uint64_t> joined_conflict(vlen + vlen2);
//
//       for(int k = 0; k < vlen + vlen2; k++) joined_pos[k] = path_pos1[k] | path_pos2[k];
//
//       int cases = 0;
//       int controls = 0;
//       for(int k = 0; k < vlen; k++) cases += bitCount(joined_pos[k]);
//       for(int k = vlen; k < vlen + vlen2; k++) controls += bitCount(joined_pos[k]);
//
//       double score = ValueTable(cases,controls);
//       double flipped_score = ValueTable(controls,cases);
//
//       VectorPriorityQueues::iterator it = topKScoresQs.begin();
//
//       if(score > min_scores[0]){
//         (*it).push(score);
//         (*it).pop();
//         min_scores[0] = (*it).top();
//       }
//       if(flipped_score > min_scores[0]){
//         (*it).push(score);
//         (*it).pop();
//         min_scores[0] = (*it).top();
//       }
//
//       for(int m = 0; m < iterations; m++){
//         std::vector<uint64_t> &caseorcontrol = CaseORControl2[m];
//         int cases_m = 0;
//         int controls_m = 0;
//         for(int k = 0; k < vlen; k++){
//           uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
//           cases_m += bitCount(permuted_path_k);
//         }
//         for(int k = vlen; k < vlen + vlen2; k++){
//           uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
//           controls_m += bitCount(permuted_path_k);
//         }
//         cases_m += controls - controls_m;
//         controls_m += cases - cases_m;
//         VectorPriorityQueues::iterator it2 = topKScoresQs.begin() + m + 1;
//         double score = ValueTable(cases_m,controls_m);
//         if(score > min_scores[m+1]){
//           (*it2).push(score);
//           (*it2).pop();
//           min_scores[m+1] = (*it2).top();
//         }
//         score = ValueTable(controls_m,cases_m);
//         if(score > min_scores[m+1]){
//           (*it2).push(score);
//           (*it2).pop();
//           min_scores[m+1] = (*it2).top();
//         }
//       }
//
//     }
//
//   }
//
//   int i = 0;
//   int j = 0;
//   for(VectorPriorityQueues::iterator it = topKScoresQs.begin(); it != topKScoresQs.end(); ++it){
//     std::priority_queue<double, std::vector<double>, std::greater<double> > Queue = *it;
//     while(Queue.size() != 0){
//       Scores4QueueInitiation(i,j) = Queue.top();
//       Queue.pop();
//       j++;
//     }
//     j = 0;
//     i++;
//   }
//
//   return Scores4QueueInitiation;
//
// }



// /**
//  *
//  * Convert paths from the 64 bit representation to a path data structure.
//  *
//  */
// void StoreJoinedPath(std::vector<uint64_t> &joined_path, std::vector<short unsigned int> &cols,
//                      std::vector<int> &counts, std::vector<int> &locations,
//                      int &ncols, int &path_number, int vlen, int vlen2, int nCases, int nControls){
//
//   int count = 0;
//
//   for(int l = 0; l < vlen; l++){
//
//     uint64_t temp = joined_path[l];
//
//     for(int m = 0; m < 64; m++){
//
//       uint64_t mask = temp & (one_64bit << m);
//
//       if(mask != 0){
//         cols.push_back(l*64 + m);
//         count++;
//       }
//
//     }
//
//   }
//
//   for(int l = 0; l < vlen2; l++){
//
//     uint64_t temp = joined_path[vlen + l];
//
//     for(int m = 0; m < 64; m++){
//
//       uint64_t mask = temp & (one_64bit << m);
//
//       if(mask != 0){
//         cols.push_back(l*64 + m + nCases);
//         count++;
//       }
//
//     }
//
//   }
//
//   counts[path_number] = count;
//   locations[path_number] = ncols;
//   ncols += count;
//
// }




// /**
//  *
//  * Join the path data structures for Method 1.
//  *
//  */
// // [[Rcpp::export]]
// List JoinIndicesMethod1(IntegerVector counts_neg1, IntegerVector cols_neg1,
//                         IntegerVector counts_pos1, IntegerVector cols_pos1,
//                         IntegerVector counts_conflict1, IntegerVector cols_conflict1,
//                         IntegerVector counts_neg2, IntegerVector cols_neg2, IntegerVector locations_neg2,
//                         IntegerVector counts_pos2, IntegerVector cols_pos2, IntegerVector locations_pos2,
//                         IntegerVector counts_conflict2, IntegerVector cols_conflict2, IntegerVector locations_conflict2,
//                         IntegerVector srcuid, IntegerVector trguids2, List uids_CountLoc, IntegerVector joining_gene_sign,
//                         NumericMatrix ValueTable, NumericMatrix queues4init, int nCases, int nControls, int K,
//                         int iterations, IntegerMatrix CaseORControl, int method, int pathLength, int nthreads){
//
//
//   int vlen = (int) ceil(nCases/64.0);
//   int vlen2 = (int) ceil(nControls/64.0);
//
//
//   std::vector<std::vector<double> > ValueTable2(ValueTable.nrow(), std::vector<double>(ValueTable.ncol()));
//   for(int i = 0; i < ValueTable.nrow(); i++)
//     for(int j = 0; j < ValueTable.ncol(); j++)
//       ValueTable2[i][j] = ValueTable(i,j);
//
//   // Store CaseORControl information into 64bit integers
//   std::vector<std::vector<uint64_t> > CaseORControl2 = parseCaseORControl(CaseORControl, nCases, nControls);
//
//   // Copy data into bitsets for increased efficiency
//   // std::string pos_path1 = path + "paths_pos1_" + IntToString(pathLength);
//   std::vector<std::vector<uint64_t> > paths_pos1 = parsePaths(cols_pos1, counts_pos1, nCases, nControls);
//   // std::string pos_path1 = path + "paths_pos1_" + IntToString(pathLength);
//   std::vector<std::vector<uint64_t> > paths_pos2 = parsePaths(cols_pos2, counts_pos2, nCases, nControls);
//
//   // Vector of iterations number of  priority queues for top K scores of the permutated paths
//   VectorPriorityQueues topKScoresQs = InitQueues(iterations, K, queues4init);
//   // VectorPriorityQueues topKScoresQs = InitQueues(iterations, K);
//
//   // Create a vector to maintain the the minimum scores in each of the priority queues in topKScoresQ
//   std::vector<double> min_scores(iterations + 1);
//   if(queues4init.nrow() == 0){
//     for(int i = 0; i <= iterations; i++) min_scores[i] = -INFINITY;
//   }else{
//     for(int i = 0; i <= iterations; i++) min_scores[i] = queues4init(i,0);
//   }
//
//   // A priority queue for the indices of the top K paths
//   IndicesScoresQueue indicesQ;
//   for(int j = 0; j < K; j++) indicesQ.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));
//
//   // If pathLength <=3 then we will store the joining of the different paths
//   int temp1 = 0;
//   int temp2 = 0;
//   int temp3 = 0;
//   int &ncols_pos = temp1;
//   int &ncols_neg = temp2;
//   int &ncols_conflict = temp3;
//   int total_paths = (pathLength > 3) ? 0 : getTotalPaths(trguids2, uids_CountLoc);
//   std::vector<int> new_counts_pos(total_paths);
//   std::vector<int> new_counts_neg(total_paths);
//   std::vector<int> new_counts_conflict(total_paths);
//   std::vector<int> new_locations_pos(total_paths);
//   std::vector<int> new_locations_neg(total_paths);
//   std::vector<int> new_locations_conflict(total_paths);
//   std::vector<short unsigned int> new_cols_pos;
//   std::vector<short unsigned int> new_cols_neg;
//   std::vector<short unsigned int> new_cols_conflict;
//
//   int flipped_pivot_length = counts_pos2.size();
//   int prev_src = -1;
//   int total_srcs = 0;
//   int temp4 = 0;
//   int &path_number = temp4;
//   int total_hits = 0;
//
//   if(nthreads == -1){
//     int max_available_cores = omp_get_num_procs();
//     omp_set_num_threads(max_available_cores);
//     nthreads = max_available_cores;
//   }
//
//   // std::vector<VectorPriorityQueues> LocalQueues;
//   std::vector<IndicesScoresQueue> LocalIndicesQ;
//   for(int k = 0; k < nthreads; k++){
//     // VectorPriorityQueues queue = InitQueues(iterations, K, queues4init);
//     // LocalQueues.push_back(queue);
//     IndicesScoresQueue localindicesq;
//     for(int j = 0; j < K; j++) localindicesq.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));
//     LocalIndicesQ.push_back(localindicesq);
//   }
//
//   std::vector<std::vector<double> > max_scores(nthreads, std::vector<double>(iterations, -INFINITY));
//   // if(queues4init.nrow() != 0){
//   //   for(int i = 0; i < nthreads; i++){
//   //     int k = 0;
//   //     for(int j = 0; j < iterations; j++, k++){
//   //       int ind1 = k+1;
//   //       int ind2 = K-1;
//   //       max_scores[i][j] = queues4init(ind1,ind2);
//   //     }
//   //   }
//   // }
//
//   for(int i = 0; i < trguids2.size(); i++){
//
//     // std::cout <<  i << std::endl;
//
//     // Rcpp::checkUserInterrupt();
//     // Check if source node has changed
//     int src = srcuid[i];
//     if(src != prev_src){
//       prev_src = src;
//       total_srcs++;
//       if((total_srcs % 100) == 0){
//         Rcout << "    " <<  total_srcs << " source nodes for paths of length " << pathLength << " and their permutations have been processed!" << std::endl;
//         // std::cout << total_hits << std::endl;
//         checkUserInterrupt();
//       }
//     }
//
//     // Find the locations and the number of the paths matching with uid
//     int uid = trguids2[i];
//     std::string geneuid = IntToString(uid); // convert integer to string
//     IntegerVector uid_count_loc = uids_CountLoc[geneuid];
//     int count = uid_count_loc[0];
//     if(count == 0) continue;
//     int location = uid_count_loc[1];
//
//     // Get the data of the first path
//     std::vector<uint64_t> &path_pos1 = paths_pos1[i];
//
//     std::vector<std::vector<double> > local_max_scores(nthreads, std::vector<double>(iterations));
//     for(int j = 0; j < nthreads; j++){
//       for(int k = 0; k < iterations; k++)
//         local_max_scores[j][k] = max_scores[j][k];
//     }
//
//     // std::vector<std::vector<double> > &local_max_scores = max_scores;
//
//     #pragma omp parallel for schedule(dynamic,1) shared(path_pos1,paths_pos2,flipped_pivot_length,location,count,ValueTable2,LocalIndicesQ,local_max_scores,CaseORControl2,new_cols_pos,new_counts_pos,new_locations_pos,ncols_pos,path_number,vlen,vlen2,nCases,nControls,iterations) if(pathLength > 3)
//     for(int j = location; j < (location + count); j++){
//
//       int tid = omp_get_thread_num();
//       // VectorPriorityQueues &tid_Queues = LocalQueues[tid];
//       IndicesScoresQueue &tid_localindicesq = LocalIndicesQ[tid];
//
//       std::vector<uint64_t> joined_pos(path_pos1.size());
//       std::vector<double> scores(iterations);
//       std::vector<double> flipped_scores(iterations);
//
//       std::vector<uint64_t> &path_pos2 = paths_pos2[j];
//       for(int k = 0; k < vlen + vlen2; k++) joined_pos[k] = path_pos1[k] | path_pos2[k];
//
//       if(pathLength <= 3){
//         StoreJoinedPath(joined_pos, new_cols_pos, new_counts_pos, new_locations_pos, ncols_pos, path_number, vlen, vlen2, nCases, nControls);
//         path_number++;
//       }
//
//       int cases = 0;
//       int controls = 0;
//       for(int k = 0; k < vlen; k++) cases += bitCount(joined_pos[k]);
//       for(int k = vlen; k < vlen + vlen2; k++) controls += bitCount(joined_pos[k]);
//
//       double score = ValueTable2[cases][controls];
//       double flipped_score = ValueTable2[controls][cases];
//
//       if(score > tid_localindicesq.top().first){
//         tid_localindicesq.push(std::pair<double, std::pair<int, int> >(score, std::pair<int,int>(i, j)));
//         tid_localindicesq.pop();
//       }
//
//       if(flipped_score > tid_localindicesq.top().first){
//         tid_localindicesq.push(std::pair<double, std::pair<int, int> >(flipped_score, std::pair<int,int>(i, j + flipped_pivot_length)));
//         tid_localindicesq.pop();
//       }
//
//       std::vector<double> &tid_max_scores = local_max_scores[tid];
//
//       for(int m = 0; m < iterations; m++){
//
//         double &max_score = tid_max_scores[m];
//         double max_score2 = tid_max_scores[m];
//         double perm_score = 0;
//         double perm_flipped_score = 0;
//
//         std::vector<uint64_t> &caseorcontrol = CaseORControl2[m];
//         int cases_m = 0;
//         int controls_m = 0;
//
//         for(int k = 0; k < vlen; k++){
//           uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
//           cases_m += bitCount(permuted_path_k);
//         }
//
//         for(int k = vlen; k < vlen + vlen2; k++){
//           uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
//           controls_m += bitCount(permuted_path_k);
//         }
//
//         int new_cases_m = cases_m + (controls - controls_m);
//         int new_controls_m = controls_m + (cases - cases_m);
//
//         perm_score = ValueTable2[new_cases_m][new_controls_m];
//         if(perm_score > max_score2){
//           max_score = perm_score;
//           max_score2 = perm_score;
//           // total_hits++;
//         }
//         perm_flipped_score = ValueTable2[new_controls_m][new_cases_m];
//         if(perm_flipped_score > max_score2){
//           max_score = perm_flipped_score;
//         }
//
//         // if(score > (*it).top()){
//         //   (*it).push(score);
//         //   (*it).pop();
//         //   total_hits++;
//         // }
//         // double &tid_max_score = tid_max_scores[m];
//         // if(score > tid_max_scores[m]){
//         //   tid_max_scores[m] = score;
//         //   // total_hits++;
//         // }
//         // max_score = perm_score;
//         // scores[m] = perm_score;
//         // flipped_scores[m] = perm_flipped_score;
//         // if(score > tid_max_scores[m]){
//         //   tid_max_scores[m] = score;
//         // }
//         // if(score > (*it).top()){
//         //   (*it).push(score);
//         //   (*it).pop();
//         //   total_hits++;
//         // }
//
//       }
//
//     }
//
//     for(int j = 0; j < nthreads; j++){
//       for(int k = 0; k < iterations; k++)
//         if(max_scores[j][k] < local_max_scores[j][k]) max_scores[j][k] = local_max_scores[j][k];
//     }
//
//   }
//
//   // for(int j = 0; j < LocalQueues.size(); j++){
//   //   VectorPriorityQueues &localqueue = LocalQueues[j];
//   //   VectorPriorityQueues::iterator it = localqueue.begin();
//   //   VectorPriorityQueues::iterator it2 = topKScoresQs.begin();
//   //   for(int k = 0; k < topKScoresQs.size(); k++,it++,it2++){
//   //     while((*it).size() != 0) {
//   //       (*it2).push((*it).top());
//   //       (*it).pop();
//   //       (*it2).pop();
//   //     }
//   //   }
//   // }
//
//   NumericVector permutedscores(iterations);
//   for(int i = 0; i < nthreads; i++){
//     for(int j = 0; j < iterations; j++){
//       if(i == 0){
//         permutedscores[j] = max_scores[i][j];
//       }else if(max_scores[i][j] > permutedscores[j]){
//         permutedscores[j] = max_scores[i][j];
//       }
//     }
//   }
//
//   for(int j = 0; j < LocalIndicesQ.size(); j++){
//     IndicesScoresQueue &localqueue = LocalIndicesQ[j];
//     while(localqueue.size() != 0) {
//       indicesQ.push(localqueue.top());
//       localqueue.pop();
//       indicesQ.pop();
//     }
//   }
//
//   // std::cout << total_hits << std::endl;
//   // // Get the vector of all the top k scores including scores computed from permuted data
//   // std::vector<std::priority_queue<double, std::vector<double>, std::greater<double> > >::iterator it = topKScoresQs.begin();
//   // NumericVector permutedscores(K*(iterations + 1));
//   // int i = 0;
//   // for(it = topKScoresQs.begin(); it != topKScoresQs.end(); ++it){
//   //   std::priority_queue<double, std::vector<double>, std::greater<double> > Queue = *it;
//   //   while(Queue.size() != 0){
//   //     permutedscores[i] = Queue.top();
//   //     Queue.pop();
//   //     i++;
//   //   }
//   // }
//
//   // Create a vector for the indices of the k paths with the highest scores
//   IntegerMatrix ids(indicesQ.size(),2);
//   NumericVector scores(indicesQ.size());
//   int j = 0;
//   while(indicesQ.size() != 0){
//     ScoreIndices temp = indicesQ.top();
//     // permutedscores[i] = temp.first;
//     scores[j] = temp.first;
//     ids(j,0) = temp.second.first + 1; // Add 1 because it will be used as an index in the R code
//     ids(j,1) = temp.second.second + 1;
//     indicesQ.pop();
//     j++;
//     // i++;
//   }
//
//   // std::cout << "total hits are " << total_hits << std::endl;
//   if(pathLength > 3){
//     return List::create(Named("scores") = scores,
//                         Named("ids") = ids,
//                         Named("TestScores") = permutedscores);
//
//   }else{
//
//     // Include the computed joined paths
//     IntegerVector cols_pos = IntegerVector(new_cols_pos.begin(), new_cols_pos.end());
//     IntegerVector cols_neg = IntegerVector(new_cols_neg.begin(), new_cols_neg.end());
//     IntegerVector cols_conflict = IntegerVector(new_cols_conflict.begin(), new_cols_conflict.end());
//     IntegerVector counts_pos = IntegerVector(new_counts_pos.begin(), new_counts_pos.end());
//     IntegerVector counts_neg = IntegerVector(new_counts_neg.begin(), new_counts_neg.end());
//     IntegerVector counts_conflict = IntegerVector(new_counts_conflict.begin(), new_counts_conflict.end());
//     IntegerVector locations_pos = IntegerVector(new_locations_pos.begin(), new_locations_pos.end());
//     IntegerVector locations_neg = IntegerVector(new_locations_neg.begin(), new_locations_neg.end());
//     IntegerVector locations_conflict = IntegerVector(new_locations_conflict.begin(), new_locations_conflict.end());
//
//     List cummulation_vectors = List::create(Named("counts_pos") = counts_pos,
//                                             Named("counts_neg") = counts_neg,
//                                             Named("counts_conflict") = counts_conflict,
//                                             Named("locations_pos") = locations_pos,
//                                             Named("locations_neg") = locations_neg,
//                                             Named("locations_conflict") = locations_conflict,
//                                             Named("cols_pos") = cols_pos,
//                                             Named("cols_neg") = cols_neg,
//                                             Named("cols_conflict") = cols_conflict);
//     return List::create(Named("scores") = scores,
//                         Named("ids") = ids,
//                         Named("TestScores") = permutedscores,
//                         Named("cummulation_vectors") = cummulation_vectors);
//   }
//
//
// }



// /**
//  *
//  * Join the path data structures for Method 2.
//  *
//  */
// // [[Rcpp::export]]
// List JoinIndicesMethod2(IntegerVector counts_neg1, IntegerVector cols_neg1,
//                  IntegerVector counts_pos1, IntegerVector cols_pos1,
//                  IntegerVector counts_conflict1, IntegerVector cols_conflict1,
//                  IntegerVector counts_neg2, IntegerVector cols_neg2, IntegerVector locations_neg2,
//                  IntegerVector counts_pos2, IntegerVector cols_pos2, IntegerVector locations_pos2,
//                  IntegerVector counts_conflict2, IntegerVector cols_conflict2, IntegerVector locations_conflict2,
//                  IntegerVector srcuid, IntegerVector trguids2, List uids_CountLoc, IntegerVector joining_gene_sign,
//                  NumericMatrix ValueTable, NumericMatrix queues4init, int nCases, int nControls, int K,
//                  int iterations, IntegerMatrix CaseORControl, int method, int pathLength, int nthreads){
//
//
//   int vlen = (int) ceil(nCases/64.0);
//   int vlen2 = (int) ceil(nControls/64.0);
//
//   std::vector<std::vector<double> > ValueTable2(ValueTable.nrow(), std::vector<double>(ValueTable.ncol()));
//   for(int i = 0; i < ValueTable.nrow(); i++)
//     for(int j = 0; j < ValueTable.ncol(); j++)
//       ValueTable2[i][j] = ValueTable(i,j);
//
//   // Store CaseORControl information into 64bit integers
//   std::vector<std::vector<uint64_t> > CaseORControl2 = parseCaseORControl(CaseORControl, nCases, nControls);
//
//   // Copy data into bitsets for increased efficiency  std::string path("path");
//   std::string path("path");
//   std::vector<std::vector<uint64_t> > paths_pos1 = parsePaths(cols_pos1, counts_pos1, nCases, nControls);
//   std::vector<std::vector<uint64_t> > paths_pos2 = parsePaths(cols_pos2, counts_pos2, nCases, nControls);
//   std::vector<std::vector<uint64_t> > paths_neg1 = parsePaths(cols_neg1, counts_neg1, nCases, nControls);
//   std::vector<std::vector<uint64_t> > paths_neg2 = parsePaths(cols_neg2, counts_neg2, nCases, nControls);
//   std::vector<std::vector<uint64_t> > paths_conflict1 = parsePaths(cols_conflict1, counts_conflict1, nCases, nControls);
//   std::vector<std::vector<uint64_t> > paths_conflict2 = parsePaths(cols_conflict2, counts_conflict2, nCases, nControls);
//
//   // // Vector of iterations number of  priority queues for top K scores of the permutated paths
//   // VectorPriorityQueues topKScoresQs = InitQueues(iterations, K, queues4init);
//   // // VectorPriorityQueues topKScoresQs = InitQueues(iterations, K);
//
//   // // Create a vector to maintain the the minimum scores in each of the priority queues in topKScoresQ
//   // std::vector<double> min_scores(iterations + 1);
//   // if(queues4init.nrow() == 0){
//   //   for(int i = 0; i <= iterations; i++) min_scores[i] = -INFINITY;
//   // }else{
//   // for(int i = 0; i <= iterations; i++) min_scores[i] = queues4init(i,0);
//   // }
//
//   // A priority queue for the indices of the top K paths
//   IndicesScoresQueue indicesQ;
//   for(int j = 0; j < K; j++) indicesQ.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));
//
//   // If pathLength <=3 then we will store the joining of the different paths
//   int ncols_pos = 0;
//   int ncols_neg = 0;
//   int ncols_conflict = 0;
//   int total_paths = (pathLength > 3) ? 0 : getTotalPaths(trguids2, uids_CountLoc);
//   std::vector<int> new_counts_pos(total_paths);
//   std::vector<int> new_counts_neg(total_paths);
//   std::vector<int> new_counts_conflict(total_paths);
//   std::vector<int> new_locations_pos(total_paths);
//   std::vector<int> new_locations_neg(total_paths);
//   std::vector<int> new_locations_conflict(total_paths);
//   std::vector<short unsigned int> new_cols_pos;
//   std::vector<short unsigned int> new_cols_neg;
//   std::vector<short unsigned int> new_cols_conflict;
//
//   int flipped_pivot_length = counts_pos2.size();
//   int prev_src = -1;
//   int total_srcs = 0;
//   int temp4 = 0;
//   int &path_number = temp4;
//
//   // std::vector<VectorPriorityQueues> LocalQueues;
//   std::vector<IndicesScoresQueue> LocalIndicesQ;
//   for(int k = 0; k < nthreads; k++){
//     // VectorPriorityQueues queue = InitQueues(iterations, K, queues4init);
//     // LocalQueues.push_back(queue);
//     IndicesScoresQueue localindicesq;
//     for(int j = 0; j < K; j++) localindicesq.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));
//     LocalIndicesQ.push_back(localindicesq);
//   }
//
//   std::vector<std::vector<double> > max_scores(nthreads, std::vector<double>(iterations, -INFINITY));
//
//   for(int i = 0; i < trguids2.size(); i++){
//
//     // Check if source node has changed
//     int src = srcuid[i];
//     if(src != prev_src){
//       prev_src = src;
//       total_srcs++;
//       if((total_srcs % 100) == 0){
//         Rcout << "    " <<  total_srcs << " source nodes for paths of length " << pathLength << " and their permutations have been processed!" << std::endl;
//         checkUserInterrupt();
//       }
//     }
//
//     // Find the locations and the number of the paths matching with uid
//     int uid = trguids2[i];
//     std::string geneuid = IntToString(uid); // convert integer to string
//     IntegerVector uid_count_loc = uids_CountLoc[geneuid];
//     int count = uid_count_loc[0];
//     if(count == 0) continue;
//     int location = uid_count_loc[1];
//     int sign;
//     if(pathLength > 3) sign = joining_gene_sign[i];
//
//     // Get the data of the first path
//     std::vector<uint64_t> &path_pos1 = paths_pos1[i];
//     std::vector<uint64_t> &path_neg1 = paths_neg1[i];
//     std::vector<uint64_t> &path_conflict1 = paths_conflict1[i];
//
//     // std::vector<std::vector<double> > local_max_scores(nthreads, std::vector<double>(iterations, -INFINITY));
//
//     std::vector<std::vector<double> > local_max_scores(nthreads, std::vector<double>(iterations));
//     // std::vector<std::vector<double> > local_max_scores = max_scores;
//     for(int j = 0; j < nthreads; j++)
//       for(int k = 0; k < iterations; k++)
//         local_max_scores[j][k] = max_scores[j][k];
//
//     #pragma omp parallel for schedule(dynamic,1) shared(path_pos1,path_neg1,path_conflict1,paths_pos2,paths_neg2,paths_conflict2,flipped_pivot_length,location,count,local_max_scores,ValueTable2,CaseORControl2,new_cols_pos,new_counts_pos,new_locations_pos,ncols_pos,path_number,vlen,vlen2,nCases,nControls,iterations,sign) if(pathLength > 3) num_threads(10)
//     for(int j = location; j < (location + count); j++){
//
//
//       int tid = omp_get_thread_num();
//       // VectorPriorityQueues &tid_Queues = LocalQueues[tid];
//       IndicesScoresQueue &tid_localindicesq = LocalIndicesQ[tid];
//       std::vector<double> &tid_max_scores = local_max_scores[tid];
//
//       if(pathLength < 3) sign = joining_gene_sign[j];
//
//       if(pathLength == 3){
//         int sign2 = joining_gene_sign[i];
//         int sign3 = joining_gene_sign[j];
//         sign = ((sign2 == -1 && sign3 == 1) || (sign2 == 1 && sign3 == -1)) ? -1 : 1;
//       }
//
//       std::vector<uint64_t> joined_pos(path_pos1.size());
//       std::vector<uint64_t> joined_neg(path_pos1.size());
//       std::vector<uint64_t> joined_conflict(path_pos1.size());
//
//       std::vector<uint64_t> &path_pos2 = (sign == 1) ? paths_pos2[j] : paths_neg2[j];
//       std::vector<uint64_t> &path_neg2 = (sign == 1) ? paths_neg2[j] : paths_pos2[j];
//       std::vector<uint64_t> &path_conflict2 = paths_conflict2[j];
//
//       for(int k = 0; k < vlen + vlen2; k++){
//         uint64_t temp_pos = path_pos1[k] | path_pos2[k];
//         uint64_t temp_neg = path_neg1[k] | path_neg2[k];
//         uint64_t temp_conflict = (path_conflict1[k] | path_conflict2[k]) | (temp_pos & temp_neg);
//         joined_conflict[k] = temp_conflict;
//         joined_pos[k] = temp_pos ^ (temp_conflict & temp_pos);
//         joined_neg[k] = temp_neg ^ (temp_conflict & temp_neg);
//       }
//
//       if(pathLength <= 3){
//         StoreJoinedPath(joined_pos, new_cols_pos, new_counts_pos, new_locations_pos, ncols_pos, path_number, vlen, vlen2, nCases, nControls);
//         StoreJoinedPath(joined_neg, new_cols_neg, new_counts_neg, new_locations_neg, ncols_neg, path_number, vlen, vlen2, nCases, nControls);
//         StoreJoinedPath(joined_conflict, new_cols_conflict, new_counts_conflict, new_locations_conflict, ncols_conflict, path_number, vlen, vlen2, nCases, nControls);
//         path_number++;
//       }
//
//       int cases_pos = 0;
//       int cases_neg = 0;
//       int controls_pos = 0;
//       int controls_neg = 0;
//       for(int k = 0; k < vlen; k++){
//         cases_pos += bitCount(joined_pos[k]);
//         controls_neg += bitCount(joined_neg[k]);
//       }
//       for(int k = vlen; k < vlen + vlen2; k++){
//         controls_pos += bitCount(joined_pos[k]);
//         cases_neg += bitCount(joined_neg[k]);
//       }
//
//       int cases = cases_pos + cases_neg;
//       int controls = controls_pos + controls_neg;
//
//       double score = ValueTable2[cases][controls];
//       double flipped_score = ValueTable2[controls][cases];
//
//       if(score > tid_localindicesq.top().first){
//         tid_localindicesq.push(std::pair<double, std::pair<int, int> >(score, std::pair<int,int>(i, j)));
//         tid_localindicesq.pop();
//       }
//
//       if(flipped_score > tid_localindicesq.top().first){
//         tid_localindicesq.push(std::pair<double, std::pair<int, int> >(flipped_score, std::pair<int,int>(i, j + flipped_pivot_length)));
//         tid_localindicesq.pop();
//       }
//
//       // if(score > min_scores[0]){
//       //   #pragma omp critical(dataupdate)
//       //   {
//       //     indicesQ.push(std::pair<double, std::pair<int, int> >(score, std::pair<int,int>(i, j)));
//       //     indicesQ.pop();
//       //     min_scores[0] = indicesQ.top().first;
//       //   }
//       // }
//       // if(flipped_score > min_scores[0]){
//       //   #pragma omp critical(dataupdate)
//       //   {
//       //     indicesQ.push(std::pair<double, std::pair<int, int> >(flipped_score, std::pair<int,int>(i, j + flipped_pivot_length)));
//       //     indicesQ.pop();
//       //     min_scores[0] = indicesQ.top().first;
//       //   }
//       // }
//
//       for(int m = 0; m < iterations; m++){
//
//         std::vector<uint64_t> &caseorcontrol = CaseORControl2[m];
//         int cases_pos_m = 0;
//         int controls_pos_m = 0;
//         int cases_neg_m = 0;
//         int controls_neg_m = 0;
//         for(int k = 0; k < vlen; k++){
//           uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
//           cases_pos_m += bitCount(permuted_path_k);
//           permuted_path_k = joined_neg[k] & caseorcontrol[k];
//           controls_neg_m += bitCount(permuted_path_k);
//         }
//         for(int k = vlen; k < vlen + vlen2; k++){
//           uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
//           controls_pos_m += bitCount(permuted_path_k);
//           permuted_path_k = joined_neg[k] & caseorcontrol[k];
//           cases_neg_m += bitCount(permuted_path_k);
//         }
//
//         int new_cases = cases_pos_m + cases_neg_m + (controls_pos - controls_pos_m) + (controls_neg - controls_neg_m);
//         int new_controls = controls_pos_m + controls_neg_m + (cases_pos - cases_pos_m) + (cases_neg - cases_neg_m);
//
//         // VectorPriorityQueues::iterator it = topKScoresQs.begin() + m;
//         double score = ValueTable2[new_cases][new_controls];
//         // if(score > min_scores[m+1]){
//         //   #pragma omp critical(dataupdate)
//         //   {
//         //     (*it).push(score);
//         //     (*it).pop();
//         //     min_scores[m+1] = (*it).top();
//         //   }
//         // }
//
//         double &tid_max_score = tid_max_scores[m];
//         if(score > tid_max_score){
//           tid_max_score = score;
//         }
//         score = ValueTable2[new_controls][new_cases];
//         if(score > tid_max_score){
//           tid_max_score = score;
//         }
//         // if(score > min_scores[m+1]){
//         //   #pragma omp critical(dataupdate)
//         //   {
//         //     (*it).push(score);
//         //     (*it).pop();
//         //     min_scores[m+1] = (*it).top();
//         //   }
//         // }
//
//       }
//
//     }
//
//     for(int j = 0; j < nthreads; j++){
//       for(int k = 0; k < iterations; k++){
//         if(max_scores[j][k] < local_max_scores[j][k]) max_scores[j][k] = local_max_scores[j][k];
//       }
//     }
//
//   }
//
//
//   NumericVector permutedscores(iterations);
//   for(int i = 0; i < nthreads; i++){
//     for(int j = 0; j < iterations; j++){
//       if(i == 0){
//         permutedscores[j] = max_scores[i][j];
//       }else if(max_scores[i][j] > permutedscores[j]){
//         permutedscores[j] = max_scores[i][j];
//       }
//     }
//   }
//
//   for(int j = 0; j < LocalIndicesQ.size(); j++){
//     IndicesScoresQueue &localqueue = LocalIndicesQ[j];
//     while(localqueue.size() != 0) {
//       indicesQ.push(localqueue.top());
//       localqueue.pop();
//       indicesQ.pop();
//     }
//   }
//
//   // Get the vector of all the top k scores including scores computed from permuted data
//   // std::vector<std::priority_queue<double, std::vector<double>, std::greater<double> > >::iterator it = topKScoresQs.begin();
//   // NumericVector permutedscores(K*(iterations + 1));
//   // int i = 0;
//   // for(it = topKScoresQs.begin(); it != topKScoresQs.end(); ++it){
//   //   std::priority_queue<double, std::vector<double>, std::greater<double> > Queue = *it;
//   //   while(Queue.size() != 0){
//   //     permutedscores[i] = Queue.top();
//   //     Queue.pop();
//   //     i++;
//   //   }
//   // }
//
//   // Create a vector for the indices of the k paths with the highest scores
//   IntegerMatrix ids(indicesQ.size(),2);
//   NumericVector scores(indicesQ.size());
//   int j = 0;
//   while(indicesQ.size() != 0){
//     ScoreIndices temp = indicesQ.top();
//     // permutedscores[i] = temp.first;
//     scores[j] = temp.first;
//     ids(j,0) = temp.second.first + 1; // Add 1 because it will be used as an index in the R code
//     ids(j,1) = temp.second.second + 1;
//     indicesQ.pop();
//     j++;
//     // i++;
//   }
//
//   // std::cout << "total hits are " << total_hits << std::endl;
//   if(pathLength > 3){
//     return List::create(Named("scores") = scores,
//                         Named("ids") = ids,
//                         Named("TestScores") = permutedscores);
//
//   }else{
//
//     // Include the computed joined paths
//     IntegerVector cols_pos = IntegerVector(new_cols_pos.begin(), new_cols_pos.end());
//     IntegerVector cols_neg = IntegerVector(new_cols_neg.begin(), new_cols_neg.end());
//     IntegerVector cols_conflict = IntegerVector(new_cols_conflict.begin(), new_cols_conflict.end());
//     IntegerVector counts_pos = IntegerVector(new_counts_pos.begin(), new_counts_pos.end());
//     IntegerVector counts_neg = IntegerVector(new_counts_neg.begin(), new_counts_neg.end());
//     IntegerVector counts_conflict = IntegerVector(new_counts_conflict.begin(), new_counts_conflict.end());
//     IntegerVector locations_pos = IntegerVector(new_locations_pos.begin(), new_locations_pos.end());
//     IntegerVector locations_neg = IntegerVector(new_locations_neg.begin(), new_locations_neg.end());
//     IntegerVector locations_conflict = IntegerVector(new_locations_conflict.begin(), new_locations_conflict.end());
//
//     List cummulation_vectors = List::create(Named("counts_pos") = counts_pos,
//                                             Named("counts_neg") = counts_neg,
//                                             Named("counts_conflict") = counts_conflict,
//                                             Named("locations_pos") = locations_pos,
//                                             Named("locations_neg") = locations_neg,
//                                             Named("locations_conflict") = locations_conflict,
//                                             Named("cols_pos") = cols_pos,
//                                             Named("cols_neg") = cols_neg,
//                                             Named("cols_conflict") = cols_conflict);
//     return List::create(Named("scores") = scores,
//                         Named("ids") = ids,
//                         Named("TestScores") = permutedscores,
//                         Named("cummulation_vectors") = cummulation_vectors);
//   }
//
//
// }




// /**
//  *
//  * Join the path data structures.
//  *
//  */
// // [[Rcpp::export]]
// List JoinIndices3(IntegerVector counts_neg1, IntegerVector cols_neg1,
//                  IntegerVector counts_pos1, IntegerVector cols_pos1,
//                  IntegerVector counts_conflict1, IntegerVector cols_conflict1,
//                  IntegerVector counts_neg2, IntegerVector cols_neg2, IntegerVector locations_neg2,
//                  IntegerVector counts_pos2, IntegerVector cols_pos2, IntegerVector locations_pos2,
//                  IntegerVector counts_conflict2, IntegerVector cols_conflict2, IntegerVector locations_conflict2,
//                  IntegerVector srcuid, IntegerVector trguids2, List uids_CountLoc, IntegerVector joining_gene_sign,
//                  NumericMatrix ValueTable, NumericMatrix queues4init, int nCases, int nControls, int K,
//                  int iterations, IntegerMatrix CaseORControl, const int method, const int pathLength){
//
//
//   int vlen = (int) ceil(nCases/64.0);
//   int vlen2 = (int) ceil(nControls/64.0);
//
//   std::vector<std::vector<double> > ValueTable2(ValueTable.nrow(), std::vector<double>(ValueTable.ncol()));
//   for(int i = 0; i < ValueTable.nrow(); i++)
//     for(int j = 0; j < ValueTable.ncol(); j++)
//       ValueTable2[i][j] = ValueTable(i,j);
//
//   // Store CaseORControl information into 64bit integers
//   std::vector<std::vector<uint64_t> > CaseORControl2 = parseCaseORControl(CaseORControl, nCases, nControls);
//
//   // Copy data into bitsets for increased efficiency
//   std::vector<std::vector<uint64_t> > paths1 = parsePaths(cols_pos1, counts_pos1, nCases, nControls);
//   std::vector<std::vector<uint64_t> > paths2 = parsePaths(cols_pos2, counts_pos2, nCases, nControls);
//   // std::vector<std::vector<uint64_t> > paths_neg1 = parsePaths(cols_neg1, counts_neg1, nCases, nControls);
//   // std::vector<std::vector<uint64_t> > paths_neg2 = parsePaths(cols_neg2, counts_neg2, nCases, nControls);
//   // std::vector<std::vector<uint64_t> > paths_conflict1 = parsePaths(cols_conflict1, counts_conflict1, nCases, nControls);
//   // std::vector<std::vector<uint64_t> > paths_conflict2 = parsePaths(cols_conflict2, counts_conflict2, nCases, nControls);
//
//   // Vector of iterations number of  priority queues for top K scores of the permutated paths
//   VectorPriorityQueues topKScoresQs = InitQueues(iterations, K, queues4init);
//   IntegerMatrix test2(2,2);
//
//   std::map<int, std::pair<int, int> > uids_CountLoc2;
//   for(int i = 0; i < trguids2.size(); i++){
//
//     int uid = trguids2[i];
//     std::string geneuid = IntToString(uid); // convert integer to string
//     IntegerVector uid_count_loc = uids_CountLoc[geneuid];
//     int count = uid_count_loc[0];
//     int location = uid_count_loc[1];
//     uids_CountLoc2[uid] = std::pair<int, int>(count, location);
//
//   }
//
//   // Create a vector to maintain the the minimum scores in each of the priority queues in topKScoresQ
//   std::vector<double> min_scores(iterations + 1);
//   if(queues4init.nrow() == 0){
//     for(int i = 0; i <= iterations; i++) min_scores[i] = -INFINITY;
//   }else{
//     for(int i = 0; i <= iterations; i++) min_scores[i] = queues4init(i,0);
//   }
//
//   // A priority queue for the indices of the top K paths
//   IndicesScoresQueue indicesQ;
//   for(int j = 0; j < K; j++) indicesQ.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));
//
//   // If pathLength <=3 then we will store the joining of the different paths
//   int total_pos_count = 0;
//   int total_neg_count = 0;
//   int total_conflict_count = 0;
//   int temp1 = 0;
//   int &ncols_pos = temp1;
//   int ncols_neg = 0;
//   int ncols_conflict = 0;
//   int total_paths = (pathLength > 3) ? 0 : getTotalPaths(trguids2, uids_CountLoc);
//   std::vector<int> new_counts_pos(total_paths);
//   std::vector<int> new_counts_neg(total_paths);
//   std::vector<int> new_counts_conflict(total_paths);
//   std::vector<int> new_locations_pos(total_paths);
//   std::vector<int> new_locations_neg(total_paths);
//   std::vector<int> new_locations_conflict(total_paths);
//   std::vector<short unsigned int> new_cols_pos;
//   std::vector<short unsigned int> new_cols_neg;
//   std::vector<short unsigned int> new_cols_conflict;
//
//   int flipped_pivot_length = counts_pos2.size();
//   int prev_src = -1;
//   int total_srcs = 0;
//   int temp2 = 0;
//   int &path_number = temp2;
//   int total_iters = 0;
//
//   #pragma omp parallel for shared(vlen,vlen2,ncols_pos,path_number,iterations,flipped_pivot_length,new_counts_pos,new_locations_pos,new_cols_pos,trguids2,uids_CountLoc2,indicesQ,min_scores,paths1,paths2,CaseORControl2,topKScoresQs,ValueTable2) num_threads(5)
//   for(int iter = 0; iter < iterations + 1; iter++){
//
//     IntegerVector test22 = test2.row(0);
//
//     for(int i = 0; i < trguids2.size(); i++){
//
//       // Find the locations and the number of the paths matching with uid
//       int uid = trguids2[i];
//       std::pair<int, int> uid_count_loc = uids_CountLoc2[uid];
//       int count = uid_count_loc.first;
//       if(count == 0) continue;
//       int location = uid_count_loc.second;
//
//       // Get the data of the first path
//       std::vector<uint64_t> &path1 = paths1[i];
//
//       for(int j = location; j < (location + count); j++){
//
//         std::vector<uint64_t> joined_path(path1.size());
//         std::vector<uint64_t> &path2 = paths2[j];
//         for(int k = 0; k < vlen + vlen2; k++) joined_path[k] = path1[k] | path2[k];
//
//         if(pathLength <= 3 && iter == 0){
//           StoreJoinedPath(joined_path, new_cols_pos, new_counts_pos, new_locations_pos, ncols_pos, path_number, vlen, vlen2, nCases, nControls);
//           path_number++;
//         }
//
//         int cases = 0;
//         int controls = 0;
//         for(int k = 0; k < vlen; k++) cases += bitCount(joined_path[k]);
//         for(int k = vlen; k < vlen + vlen2; k++) controls += bitCount(joined_path[k]);
//
//         double score = ValueTable2[cases][controls];
//         double flipped_score = ValueTable2[controls][cases];
//
//         if(iter == 0){
//           if(score > min_scores[0]){
//             indicesQ.push(std::pair<double, std::pair<int, int> >(score, std::pair<int,int>(i, j)));
//             indicesQ.pop();
//             min_scores[0] = indicesQ.top().first;
//           }
//           if(flipped_score > min_scores[0]){
//             indicesQ.push(std::pair<double, std::pair<int, int> >(flipped_score, std::pair<int,int>(i, j + flipped_pivot_length)));
//             indicesQ.pop();
//             min_scores[0] = indicesQ.top().first;
//           }
//         }
//         else{
//
//           std::vector<uint64_t> &caseorcontrol = CaseORControl2[iter-1];
//           int cases_m = 0;
//           int controls_m = 0;
//
//           for(int k = 0; k < vlen; k++){
//             uint64_t permuted_path_k = joined_path[k] & caseorcontrol[k];
//             cases_m += bitCount(permuted_path_k);
//           }
//
//           for(int k = vlen; k < vlen + vlen2; k++){
//             uint64_t permuted_path_k = joined_path[k] & caseorcontrol[k];
//             controls_m += bitCount(permuted_path_k);
//           }
//
//           int new_cases_m = cases_m + (controls - controls_m);
//           int new_controls_m = controls_m + (cases - cases_m);
//
//           VectorPriorityQueues::iterator it = topKScoresQs.begin() + (iter - 1);
//
//           double score = ValueTable2[new_cases_m][new_controls_m];
//           if(score > min_scores[iter]){
//             (*it).push(score);
//             (*it).pop();
//             min_scores[iter] = (*it).top();
//           }
//           score = ValueTable2[new_controls_m][new_cases_m];
//           if(score > min_scores[iter]){
//             (*it).push(score);
//             (*it).pop();
//             min_scores[iter] = (*it).top();
//           }
//
//         }
//
//       }
//
//     }
//
//     #pragma omp critical(dataupdate)
//     {
//       total_iters++;
//       if((total_iters + 1) % 10 == 0){
//         std::cout << "    " <<  (total_iters + 1) << " number of permutations of paths of length " << pathLength << " have been processed!" << std::endl;
//       }
//     }
//
//   }
//
//   // Get the vector of all the top k scores including scores computed from permuted data
//   std::vector<std::priority_queue<double, std::vector<double>, std::greater<double> > >::iterator it = topKScoresQs.begin();
//   NumericVector permutedscores(K*(iterations + 1));
//   int i = 0;
//   for(it = topKScoresQs.begin(); it != topKScoresQs.end(); ++it){
//     std::priority_queue<double, std::vector<double>, std::greater<double> > Queue = *it;
//     while(Queue.size() != 0){
//       permutedscores[i] = Queue.top();
//       Queue.pop();
//       i++;
//     }
//   }
//
//   // Create a vector for the indices of the k paths with the highest scores
//   IntegerMatrix ids(indicesQ.size(),2);
//   NumericVector scores(indicesQ.size());
//   int j = 0;
//   while(indicesQ.size() != 0){
//     ScoreIndices temp = indicesQ.top();
//     permutedscores[i] = temp.first;
//     scores[j] = temp.first;
//     ids(j,0) = temp.second.first + 1; // Add 1 because it will be used as an index in the R code
//     ids(j,1) = temp.second.second + 1;
//     indicesQ.pop();
//     j++;
//     i++;
//   }
//
//   if(pathLength > 3){
//     return List::create(Named("scores") = scores,
//                         Named("ids") = ids,
//                         Named("TestScores") = permutedscores);
//
//   }else{
//
//     // Include the computed joined paths
//     IntegerVector cols_pos = IntegerVector(new_cols_pos.begin(), new_cols_pos.end());
//     IntegerVector cols_neg = IntegerVector(new_cols_neg.begin(), new_cols_neg.end());
//     IntegerVector cols_conflict = IntegerVector(new_cols_conflict.begin(), new_cols_conflict.end());
//     IntegerVector counts_pos = IntegerVector(new_counts_pos.begin(), new_counts_pos.end());
//     IntegerVector counts_neg = IntegerVector(new_counts_neg.begin(), new_counts_neg.end());
//     IntegerVector counts_conflict = IntegerVector(new_counts_conflict.begin(), new_counts_conflict.end());
//     IntegerVector locations_pos = IntegerVector(new_locations_pos.begin(), new_locations_pos.end());
//     IntegerVector locations_neg = IntegerVector(new_locations_neg.begin(), new_locations_neg.end());
//     IntegerVector locations_conflict = IntegerVector(new_locations_conflict.begin(), new_locations_conflict.end());
//
//     List cummulation_vectors = List::create(Named("counts_pos") = counts_pos,
//                                             Named("counts_neg") = counts_neg,
//                                             Named("counts_conflict") = counts_conflict,
//                                             Named("locations_pos") = locations_pos,
//                                             Named("locations_neg") = locations_neg,
//                                             Named("locations_conflict") = locations_conflict,
//                                             Named("cols_pos") = cols_pos,
//                                             Named("cols_neg") = cols_neg,
//                                             Named("cols_conflict") = cols_conflict);
//     return List::create(Named("scores") = scores,
//                         Named("ids") = ids,
//                         Named("TestScores") = permutedscores,
//                         Named("cummulation_vectors") = cummulation_vectors);
//   }
//
//
// }
