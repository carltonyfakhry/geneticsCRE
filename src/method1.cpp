#include <fstream>
#include "Utils.h"
#include "gcre.h"

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List JoinIndicesMethod1(IntegerVector srcuid, IntegerVector trguids2, List uids_CountLoc, IntegerVector joining_gene_sign,
  NumericMatrix ValueTable, int nCases, int nControls, int K,
  int iterations, IntegerMatrix CaseORControl, int pathLength, int nthreads, std::string pos_path1,
  std::string pos_path2, std::string dest_path_pos){

  int vlen = (int) ceil(nCases/64.0);
  int vlen2 = (int) ceil(nControls/64.0);

  // Copy ValueTable into a C++ 2D vector to be used in an openmp for loop
  std::vector<std::vector<double> > ValueTable2(ValueTable.nrow(), std::vector<double>(ValueTable.ncol()));
  for(int i = 0; i < ValueTable.nrow(); i++) {
    for(int j = 0; j < ValueTable.ncol(); j++) {
      ValueTable2[i][j] = ValueTable(i,j);
    }
  }

  // Store CaseORControl information into 64bit integers
  std::vector<std::vector<uint64_t> > CaseORControl2 = parseCaseORControl(CaseORControl, nCases, nControls);

  // Copy the paths into 2D vectors
  std::vector<std::vector<uint64_t> > paths_pos1 = (pos_path1 == "") ? getZeroMatrix(trguids2.size(), vlen + vlen2) : (readPaths(pos_path1));
  std::vector<std::vector<uint64_t> > temp_paths_pos2;
  std::vector<std::vector<uint64_t> > &paths_pos2 = (pathLength == 5) ? paths_pos1 : (temp_paths_pos2 = readPaths(pos_path2));

  // A priority queue for the indices of the top K paths in the data
  IndicesScoresQueue indicesQ;
  for(int j = 0; j < K; j++)
    indicesQ.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));

  // If pathLength <=3 then we will store the joining of the paths
  int total_paths = (pathLength > 3) ? 0 : getTotalPaths(trguids2, uids_CountLoc);
  std::vector<std::vector<uint64_t> > joined_paths(total_paths, std::vector<uint64_t>(vlen + vlen2, 0));

  nthreads = 1;

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

      for(int j = location; j < (location + count); j++){

        // int tid = omp_get_thread_num();
        int tid = 0;
        IndicesScoresQueue &tid_localindicesq = LocalIndicesQ[tid];
        std::vector<uint64_t> joined_pos(path_pos1.size());

        // Get the data for path2 which is to be joined with path1
        std::vector<uint64_t> &path_pos2 = paths_pos2[j];

        // Join the paths
        for(int k = 0; k < vlen + vlen2; k++)
          joined_pos[k] = path_pos1[k] | path_pos2[k];

        // Store the joined paths if pathLength <= 3
        if(pathLength <= 3){
          joined_paths[total_paths] = joined_pos;
          total_paths++;
        }

        int cases = 0;
        int controls = 0;
        for(int k = 0; k < vlen; k++)
          cases += __builtin_popcountl(joined_pos[k]);
        for(int k = vlen; k < vlen + vlen2; k++)
          controls += __builtin_popcountl(joined_pos[k]);

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

          for(int k = 0; k < vlen; k++){
            uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
            cases_m += __builtin_popcountl(permuted_path_k);
          }

          int controls_m = 0;

          for(int k = vlen; k < vlen + vlen2; k++){
            uint64_t permuted_path_k = joined_pos[k] & caseorcontrol[k];
            controls_m += __builtin_popcountl(permuted_path_k);
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

  return List::create(Named("scores") = scores,Named("ids") = ids, Named("TestScores") = permutedscores);

}
