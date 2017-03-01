#include <fstream>
#include "Utils.h"
#include "gcre.h"

List join_method2(IntegerVector srcuid, IntegerVector trguids2, List uids_CountLoc, IntegerVector joining_gene_sign,
  vecd2d value_table, int nCases, int nControls, int K,
  int iterations, IntegerMatrix CaseORControl, int pathLength, int nthreads, std::string pos_path1,
  std::string neg_path1, std::string conflict_path1, std::string pos_path2, std::string neg_path2, std::string conflict_path2,
  std::string dest_path_pos, std::string dest_path_neg, std::string dest_path_conflict){

  int vlen = (int) ceil(nCases/64.0);
  int vlen2 = (int) ceil(nControls/64.0);

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
  for(int j = 0; j < K; j++)
    indicesQ.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));

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
    std::string geneuid = std::to_string(uid);
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
    for(int j = 0; j < nthreads; j++){
      for(int k = 0; k < iterations; k++)
        local_max_scores[j][k] = max_scores[j][k];
    }
    for(int j = location; j < (location + count); j++){
      int tid = 0;
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

      double score = value_table[cases][controls];
      double flipped_score = value_table[controls][cases];

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

        perm_score = value_table[new_cases][new_controls];
        if(perm_score > max_score2){
          max_score = perm_score;
          max_score2 = perm_score;
        }
        perm_flipped_score = value_table[new_controls][new_cases];
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

  return List::create(Named("scores") = scores, Named("ids") = ids, Named("TestScores") = permutedscores);
}
