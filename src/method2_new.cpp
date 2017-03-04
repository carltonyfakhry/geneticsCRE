#include <fstream>
#include "Utils.h"
#include "gcre.h"

joined_res* join_method2(vector<int>& src_uids, vector<int>& trg_uids, Rcpp::List& uids_CountLoc, vector<int>& join_gene_signs,
  vec2d_d& value_table, int nCases, int nControls, int K,
  int iterations, vec2d_u64& case_mask, int pathLength, int nthreads,
  paths_vec* paths0, paths_vec* paths1, paths_vec* paths_res,
  int total_paths){

  std::cout.imbue(std::locale("en_US.UTF8"));

  bool keep_paths = paths_res != NULL;

  if(keep_paths){
    paths_res->size = total_paths;
    paths_res->width_ul = paths1->width_ul;
    paths_res->num_cases = nCases;
    paths_res->pos.resize(paths_res->size, vec_u64(paths_res->width_ul, 0));
    paths_res->neg.resize(paths_res->size, vec_u64(paths_res->width_ul, 0));
    paths_res->con.resize(paths_res->size, vec_u64(paths_res->width_ul, 0));
    printf("  ** resized stored paths: %d x %d\n", paths_res->size, paths_res->width_ul);
  }

  int vlen = (int) ceil(nCases/64.0);
  int vlen2 = (int) ceil(nControls/64.0);

  vec2d_u64 &paths_pos1 = paths0->pos;
  vec2d_u64 &paths_pos2 = paths1->pos;
  vec2d_u64 &paths_neg1 = paths0->neg;
  vec2d_u64 &paths_neg2 = paths1->neg;
  vec2d_u64 &paths_conflict1 = paths0->con;
  vec2d_u64 &paths_conflict2 = paths1->con;

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
  int path_idx = 0;

  long total_comps = 0;
  long zero_comps = 0;

  for(int i = 0; i < trg_uids.size(); i++){
    // Check if source node has changed
    int src = src_uids[i];
    if(src != prev_src){
      prev_src = src;
      total_srcs++;
      // if((total_srcs % 100) == 0){
      //   Rcout << "    " <<  total_srcs << " source nodes for paths of length " << pathLength << " and their permutations have been processed!" << std::endl;
      //   checkUserInterrupt();
      // }
    }

    // Find the locations and the number of the paths matching with uid
    int uid = trg_uids[i];
    std::string geneuid = std::to_string(uid);
    IntegerVector uid_count_loc = uids_CountLoc[geneuid];
    int count = uid_count_loc[0];
    if(count == 0) continue;
    int location = uid_count_loc[1];
    int sign;
    if(pathLength > 3) sign = join_gene_signs[i];

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

      if(pathLength < 3) sign = join_gene_signs[j];

      if(pathLength == 3){
        int sign2 = join_gene_signs[i];
        int sign3 = join_gene_signs[j];
        sign = ((sign2 == -1 && sign3 == 1) || (sign2 == 1 && sign3 == -1)) ? -1 : 1;
      }

      vec_u64 joined_pos(path_pos1.size());
      vec_u64 joined_neg(path_pos1.size());
      vec_u64 joined_conflict(path_pos1.size());

      // TODO ??
      std::vector<uint64_t> &path_pos2 = (sign == 1) ? paths_pos2[j] : paths_neg2[j];
      std::vector<uint64_t> &path_neg2 = (sign == 1) ? paths_neg2[j] : paths_pos2[j];
      std::vector<uint64_t> &path_conflict2 = paths_conflict2[j];

      for(int k = 0; k < vlen + vlen2; k++){
        uint64_t temp_pos = path_pos1[k] | path_pos2[k];
        uint64_t temp_neg = path_neg1[k] | path_neg2[k];

        total_comps += 1;
        if(temp_pos == 0 || temp_neg == 0)
          zero_comps += 1;

        uint64_t temp_conflict = (path_conflict1[k] | path_conflict2[k]) | (temp_pos & temp_neg);
        joined_conflict[k] = temp_conflict;
        joined_pos[k] = temp_pos ^ (temp_conflict & temp_pos);
        joined_neg[k] = temp_neg ^ (temp_conflict & temp_neg);
      }

      if(keep_paths){
        paths_res->pos[path_idx] = joined_pos;
        paths_res->neg[path_idx] = joined_neg;
        paths_res->con[path_idx] = joined_conflict;
      }

      path_idx += 1;

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

        vec_u64& caseorcontrol = case_mask[m];
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

  joined_res* res = new joined_res;

  res->permuted_scores.resize(iterations);
  for(int i = 0; i < nthreads; i++){
    for(int j = 0; j < iterations; j++){
      if(i == 0){
        res->permuted_scores[j] = max_scores[i][j];
      }else if(max_scores[i][j] > res->permuted_scores[j]){
        res->permuted_scores[j] = max_scores[i][j];
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
  // IntegerMatrix ids(indicesQ.size(),2);
  // TODO move to R wrapper
  res->ids.resize(indicesQ.size(), vec_u64(2));
  res->scores.resize(indicesQ.size());
  int j = 0;
  while(indicesQ.size() != 0){
    ScoreIndices temp = indicesQ.top();
    res->scores[j] = temp.first;
    res->ids[j][0] = temp.second.first + 1; // Add 1 because it will be used as an index in the R code
    res->ids[j][1] = temp.second.second + 1;
    indicesQ.pop();
    j++;
  }

  return res;
}
