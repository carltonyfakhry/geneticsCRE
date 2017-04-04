// #include "Utils.h"
// #include <omp.h>
// #include <fstream>
// #include <emmintrin.h>
// #include <nmmintrin.h>
// #include <immintrin.h>
// #include <smmintrin.h>
//
// // [[Rcpp::plugins(openmp)]]
// // [[Rcpp::plugins(cpp11)]]
//
//
// using namespace Rcpp;
//
//
// #define ALIGN __attribute__ (( aligned(16) ))
//
//
//
//
//
//
// /**
// *
// * Get the total number of paths to be processed.
// *
// */
//
// int getTotalPaths2(IntegerVector trguids, List uids_CountLoc){
//
//   int total_paths = 0;
//
//   for(int i = 0; i < trguids.size(); i++){
//
//     int uid = trguids[i];
//     std::string geneuid = IntToString(uid);
//     IntegerVector temp = as<IntegerVector>(uids_CountLoc[geneuid]);
//     total_paths += temp[0];
//
//   }
//
//   return total_paths;
//
// }
//
//
//
// /**
// *
// * Copy ValueTable into an aligned array.
// *
// */
//
// void copyValueTable2(double *ValueTable2, NumericMatrix ValueTable){
//
//   for(int i = 0; i < ValueTable.nrow(); i++)
//     for(int j = 0; j < ValueTable.ncol(); j++){}
//       // *ValueTable2 = ValueTable(i,j);
//
// }
//
//
//
// /**
// *
// * Convert CaseORControl from 0/1 values to their 64 bit representations and
// * store them into an aligned array.
// *
// */
//
// void parseCaseORControl2(uint64_t *CaseORControl2, IntegerMatrix CaseORControl, int nCases, int nControls, int vlen1, int vlen2){
//
//   int vlen = vlen1 + vlen2;
//
//   for(int i = 0; i < CaseORControl.nrow(); i++){
//
//     for(int j = 0; j < nCases; j++)
//       if(CaseORControl(i,j) == 1)
//         CaseORControl2[i * vlen + j/64] |= one_64bit << j % 64;
//
//     for(int j = nCases; j < nCases + nControls; j++)
//       if(CaseORControl(i,j) == 1)
//         CaseORControl2[i * vlen + vlen1 + (j-nCases)/64] |= one_64bit << (j-nCases) % 64;
//
//   }
//
// }
//
//
//
// /**
//  *
//  * Parse data into its 64 bit representation, and store it into an aligned vector.
//  *
//  */
//
// void parseData2(uint64_t *parsed_data, IntegerMatrix data, int nCases, int nControls, int vlen1, int vlen2){
//
//   int vlen = vlen1 + vlen2;
//
//   for(int i = 0; i < data.nrow(); i++)
//     for(int j = 0; j < data.ncol(); j++)
//       if(data(i,j) != 0){
//
//         if(j < nCases)
//           parsed_data[i * vlen + j/64] |= one_64bit << j % 64;
//         else
//           parsed_data[i * vlen + vlen + (j-nCases)/64] |= one_64bit << (j-nCases) % 64;
//
//       }
//
// }


// std::vector<std::vector<uint64_t> > matchData2(const std::vector<std::vector<uint64_t> > &parseddata, IntegerVector data_inds){
//
//   std::vector<std::vector<uint64_t> > matcheddata(data_inds.size(), std::vector<uint64_t>(parseddata[0].size(),0));
//
//   for(int i = 0; i < data_inds.size(); i++){
//
//     int data_index = data_inds[i];
//     matcheddata[i] = parseddata[data_index];
//
//   }
//
//   return matcheddata;
//
// }



// inline void Method1(int i, int j, int &total_paths, std::vector<uint64_t> &path_pos1, std::vector<uint64_t> &path_pos2,
//                     std::vector<std::vector<uint64_t> > &joined_paths, const int flipped_pivot_length, int location, int count,
//                     const std::vector<std::vector<double> > &ValueTable,
//                     IndicesScoresQueue &tid_localindicesq, std::vector<double> &tid_max_scores, uint64_t* p_caseorcontrol,
//                     const int vlen, const int vlen2, const int nCases, const int nControls, const int iterations, const int pathLength){
//
//   uint64_t joined_pos[path_pos1.size()] __attribute__ ((aligned (16)));
//
//   // Join the paths
//   int cases = 0;
//   int controls = 0;
//   for(int k = 0; k < vlen + vlen2; k++){
//     uint64_t temp_pos = path_pos1[k] | path_pos2[k];
//     joined_pos[k] = temp_pos;
//     if(k < vlen) cases += __builtin_popcountll(temp_pos);
//     else controls += __builtin_popcountll(temp_pos);
//     if(pathLength <= 3)
//       joined_paths[total_paths][k] = temp_pos;
//   }
//
//   double score = ValueTable[cases][controls];
//   double flipped_score = ValueTable[controls][cases];
//
//   if(score > tid_localindicesq.top().first){
//     tid_localindicesq.push(std::pair<double, std::pair<int, int> >(score, std::pair<int,int>(i, j)));
//     tid_localindicesq.pop();
//   }
//
//   if(flipped_score > tid_localindicesq.top().first){
//     tid_localindicesq.push(std::pair<double, std::pair<int, int> >(flipped_score, std::pair<int,int>(i, j + flipped_pivot_length)));
//     tid_localindicesq.pop();
//   }
//
//   __m128i *p_caseorcontrol2 = (__m128i*) p_caseorcontrol;
//
//   for(int m = 0; m < iterations; m++){
//
//     double &max_score = tid_max_scores[m];
//     double max_score2 = tid_max_scores[m];
//     double perm_score = 0;
//     double perm_flipped_score = 0;
//
//     __m128i *p_joined_pos = (__m128i*) joined_pos;
//
//     int cases_m = 0;
//     int controls_m = 0;
//
//     int k = 0;
//     for(; k < ROUND_DOWN(vlen,2); k+=2,p_joined_pos++,p_caseorcontrol2++){
//       __m128i val1_4 = _mm_load_si128(p_joined_pos);
//       __m128i val2_4 = _mm_load_si128(p_caseorcontrol2);
//       __m128i val = _mm_and_si128(val1_4, val2_4);
//       cases_m += __builtin_popcountll(_mm_extract_epi64(val, 0));
//       cases_m += __builtin_popcountll(_mm_extract_epi64(val, 1));
//     }
//     if(k < vlen){
//       cases_m += __builtin_popcountll(joined_pos[k] & *(p_caseorcontrol + m*(vlen + vlen2) + k));
//       controls_m += __builtin_popcountll(joined_pos[k + 1] & *(p_caseorcontrol + m*(vlen + vlen2) + k + 1));
//       k+=2;
//       p_joined_pos++;
//       p_caseorcontrol2++;
//     }
//
//     for(; k < ROUND_DOWN(vlen + vlen2,2); k+=2,p_joined_pos++,p_caseorcontrol2++){
//       __m128i val1_4 = _mm_load_si128(p_joined_pos);
//       __m128i val2_4 = _mm_load_si128(p_caseorcontrol2);
//       __m128i val = _mm_and_si128(val1_4, val2_4);
//       controls_m += __builtin_popcountll(_mm_extract_epi64(val, 0));
//       controls_m += __builtin_popcountll(_mm_extract_epi64(val, 1));
//     }
//     if(k < vlen + vlen2){
//       controls_m += __builtin_popcountll(joined_pos[k] & *(p_caseorcontrol + m*(vlen + vlen2) + k));
//     }
//
//     int new_cases_m = cases_m + (controls - controls_m);
//     int new_controls_m = controls_m + (cases - cases_m);
//
//     perm_score = ValueTable[new_cases_m][new_controls_m];
//     if(perm_score > max_score2){
//       max_score = perm_score;
//       max_score2 = perm_score;
//     }
//     perm_flipped_score = ValueTable[new_controls_m][new_cases_m];
//     if(perm_flipped_score > max_score2){
//       max_score = perm_flipped_score;
//     }
//
//   }
//
// }
//
//
// inline void Method2(int i, int j, int &total_paths, std::vector<uint64_t> &path_pos1, std::vector<uint64_t> &path_neg1,
//                     std::vector<uint64_t> &path_pos2, std::vector<uint64_t> &path_neg2,
//                     std::vector<std::vector<uint64_t> > &joined_paths_pos, std::vector<std::vector<uint64_t> > &joined_paths_neg,
//                     const int flipped_pivot_length, int location, int count,
//                     const std::vector<std::vector<double> > &ValueTable,
//                     IndicesScoresQueue &tid_localindicesq, std::vector<double> &tid_max_scores, uint64_t* p_caseorcontrol,
//                     const int vlen, const int vlen2, const int nCases, const int nControls, const int iterations, const int pathLength){
//
//   uint64_t joined_pos[path_pos1.size()] __attribute__ ((aligned (16)));
//   uint64_t joined_neg[path_pos1.size()] __attribute__ ((aligned (16)));
//
//   int cases_pos = 0;
//   int cases_neg = 0;
//   int controls_pos = 0;
//   int controls_neg = 0;
//
//   for(int k = 0; k < vlen + vlen2; k++){
//     uint64_t temp_pos = path_pos1[k] | path_pos2[k];
//     uint64_t temp_neg = path_neg1[k] | path_neg2[k];
//     uint64_t temp_conflict = temp_neg & temp_pos;
//     uint64_t temp_pos2 = temp_pos ^ temp_conflict;
//     uint64_t temp_neg2 = temp_neg ^ temp_conflict;
//     if(k < vlen){
//       cases_pos += __builtin_popcountll(temp_pos2);
//       controls_neg += __builtin_popcountll(temp_neg2);
//     }else{
//       controls_pos += __builtin_popcountll(temp_pos2);
//       cases_neg += __builtin_popcountll(temp_neg2);
//     }
//     joined_pos[k] = temp_pos2;
//     joined_neg[k] = temp_neg2;
//     if(pathLength <= 3){
//       joined_paths_pos[total_paths][k] = temp_pos;
//       joined_paths_neg[total_paths][k] = temp_neg;
//     }
//   }
//
//   if(pathLength <= 3){
//     total_paths++;
//   }
//
//   int cases = cases_pos + cases_neg;
//   int controls = controls_pos + controls_neg;
//
//   double score = ValueTable[cases][controls];
//   double flipped_score = ValueTable[controls][cases];
//
//   if(score > tid_localindicesq.top().first){
//     tid_localindicesq.push(std::pair<double, std::pair<int, int> >(score, std::pair<int,int>(i, j)));
//     tid_localindicesq.pop();
//   }
//
//   if(flipped_score > tid_localindicesq.top().first){
//     tid_localindicesq.push(std::pair<double, std::pair<int, int> >(flipped_score, std::pair<int,int>(i, j + flipped_pivot_length)));
//     tid_localindicesq.pop();
//   }
//
//   __m128i *p_caseorcontrol2 = (__m128i*) p_caseorcontrol;
//
//   for(int m = 0; m < iterations; m++){
//
//     double &max_score = tid_max_scores[m];
//     double max_score2 = tid_max_scores[m];
//     double perm_score = 0;
//     double perm_flipped_score = 0;
//
//     __m128i *p_joined_pos = (__m128i*) joined_pos;
//     __m128i *p_joined_neg = (__m128i*) joined_neg;
//
//     int cases_m = 0;
//     int controls_m = 0;
//
//     int k = 0;
//     for(; k < ROUND_DOWN(vlen,2); k+=2,p_joined_pos++,p_joined_neg++,p_caseorcontrol2++){
//       __m128i val1_4 = _mm_load_si128(p_joined_pos);
//       __m128i val2_4 = _mm_load_si128(p_caseorcontrol2);
//       __m128i val = _mm_and_si128(val1_4, val2_4);
//       cases_m += __builtin_popcountll(_mm_extract_epi64(val, 0));
//       cases_m += __builtin_popcountll(_mm_extract_epi64(val, 1));
//       val1_4 = _mm_load_si128(p_joined_neg);
//       val = _mm_and_si128(val1_4, val2_4);
//       controls_m += __builtin_popcountll(_mm_extract_epi64(val, 0));
//       controls_m += __builtin_popcountll(_mm_extract_epi64(val, 1));
//     }
//     if(k < vlen){
//       cases_m += __builtin_popcountll(joined_pos[k] & *(p_caseorcontrol + m*(vlen + vlen2) + k));
//       controls_m += __builtin_popcountll(joined_pos[k + 1] & *(p_caseorcontrol + m*(vlen + vlen2) + k + 1));
//       controls_m += __builtin_popcountll(joined_neg[k] & *(p_caseorcontrol + m*(vlen + vlen2) + k));
//       cases_m += __builtin_popcountll(joined_neg[k + 1] & *(p_caseorcontrol + m*(vlen + vlen2) + k + 1));
//       k+=2;
//       p_joined_pos++;
//       p_joined_neg++;
//       p_caseorcontrol2++;
//     }
//
//     for(; k < ROUND_DOWN(vlen + vlen2,2); k+=2,p_joined_pos++,p_joined_neg++,p_caseorcontrol2++){
//       __m128i val1_4 = _mm_load_si128(p_joined_pos);
//       __m128i val2_4 = _mm_load_si128(p_caseorcontrol2);
//       __m128i val = _mm_and_si128(val1_4, val2_4);
//       controls_m += __builtin_popcountll(_mm_extract_epi64(val, 0));
//       controls_m += __builtin_popcountll(_mm_extract_epi64(val, 1));
//       val1_4 = _mm_load_si128(p_joined_neg);
//       val = _mm_and_si128(val1_4, val2_4);
//       cases_m += __builtin_popcountll(_mm_extract_epi64(val, 0));
//       cases_m += __builtin_popcountll(_mm_extract_epi64(val, 1));
//     }
//     if(k < vlen + vlen2){
//       controls_m += __builtin_popcountll(joined_pos[k] & *(p_caseorcontrol + m*(vlen + vlen2) + k));
//       cases_m += __builtin_popcountll(joined_neg[k] & *(p_caseorcontrol + m*(vlen + vlen2) + k));
//     }
//
//     int new_cases_m = cases_m + (controls - controls_m);
//     int new_controls_m = controls_m + (cases - cases_m);
//
//     perm_score = ValueTable[new_cases_m][new_controls_m];
//     if(perm_score > max_score2){
//       max_score = perm_score;
//       max_score2 = perm_score;
//     }
//     perm_flipped_score = ValueTable[new_controls_m][new_cases_m];
//     if(perm_flipped_score > max_score2){
//       max_score = perm_flipped_score;
//     }
//
//   }
//
// }
//
//
// List JoinPaths(std::vector<std::vector<uint64_t> > &paths_pos1, std::vector<std::vector<uint64_t> > &paths_neg1,
//                std::vector<std::vector<uint64_t> > &paths_pos2, std::vector<std::vector<uint64_t> > &paths_neg2,
//                std::vector<std::vector<uint64_t> > &joined_paths_pos, std::vector<std::vector<uint64_t> > &joined_paths_neg,
//                IntegerVector srcuid, IntegerVector trguids, List uids_CountLoc, IntegerVector joining_gene_sign,
//                const std::vector<std::vector<double> > &ValueTable, uint64_t *p_caseorcontrol,
//                const int nCases, const int nControls, int K, const int iterations, const int method, const int pathLength, int nthreads){
//
//
//   const int vlen = (int) ceil(nCases/64.0);
//   const int vlen2 = (int) ceil(nControls/64.0);
//
//   // A priority queue for the indices of the top K paths in the data
//   IndicesScoresQueue indicesQ;
//   for(int j = 0; j < K; j++) indicesQ.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));
//
//   // nthreads can't be set to be more than the maximum number of possible threads
//   int max_nthreads = omp_get_max_threads();
//   nthreads = (nthreads > max_nthreads || nthreads == -1) ? max_nthreads : nthreads;
//   omp_set_dynamic(0);
//   omp_set_num_threads(nthreads);
//
//   // Create a priority queue which holds the top K scores per thread
//   std::vector<IndicesScoresQueue> LocalIndicesQ;
//   for(int k = 0; k < nthreads; k++){
//     IndicesScoresQueue localindicesq;
//     for(int j = 0; j < K; j++) localindicesq.push(std::make_pair<double,std::pair<int,int> >(-INFINITY, std::pair<int, int>(-1,-1)));
//     LocalIndicesQ.push_back(localindicesq);
//   }
//
//   // Global 2D vector which holds the maximum scores per permutation for each corresponding thread
//   std::vector<std::vector<double> > max_scores(nthreads, std::vector<double>(iterations, -INFINITY));
//
//   const int flipped_pivot_length = paths_pos2.size();
//   int prev_src = -1;
//   int total_srcs = 0;
//   int temp  = 0;
//   int &total_paths = temp;
//
//   for(int i = 0; i < trguids.size(); i++){
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
//     int uid = trguids[i];
//     std::string geneuid = IntToString(uid); // convert integer to string
//     IntegerVector uid_count_loc = uids_CountLoc[geneuid];
//     int count = uid_count_loc[0];
//     if(count == 0) continue;
//     int location = uid_count_loc[1];
//     int sign;
//     if(pathLength > 3 && method == 2) sign = joining_gene_sign[i];
//
//     // Get the data of the first path
//     std::vector<uint64_t> path_pos1 = paths_pos1[i];
//     std::vector<uint64_t> path_neg1;
//     if(method == 2){
//       path_neg1 = paths_neg1[i];
//     }
//
//     // Copy the values from max_scores into a local container local_max_scores to increase performance
//     std::vector<std::vector<double> > local_max_scores(nthreads, std::vector<double>(iterations));
//     for(int j = 0; j < nthreads; j++)
//       for(int k = 0; k < iterations; k++)
//         local_max_scores[j][k] = max_scores[j][k];
//
// #pragma omp parallel for shared(i,joining_gene_sign,paths_pos2,paths_neg2,location,count,ValueTable,LocalIndicesQ,local_max_scores,p_caseorcontrol) if(pathLength > 3)
//     for(int j = location; j < (location + count); j++){
//
//       int tid = omp_get_thread_num();
//       IndicesScoresQueue &tid_localindicesq = LocalIndicesQ[tid];
//       std::vector<double> &tid_max_scores = local_max_scores[tid];
//
//
//       if(method == 1){
//
//         std::vector<uint64_t> path_pos2 = paths_pos2[j];
//
//         Method1(i,j,total_paths,path_pos1,path_pos2,joined_paths_pos,flipped_pivot_length,
//                 location,count,ValueTable,tid_localindicesq,tid_max_scores,p_caseorcontrol,
//                 vlen,vlen2,nCases,nControls,iterations,pathLength);
//
//
//       }else{
//
//         if(pathLength < 3) sign = joining_gene_sign[j];
//
//         if(pathLength == 3){
//           int sign2 = joining_gene_sign[i];
//           int sign3 = joining_gene_sign[j];
//           sign = ((sign2 == -1 && sign3 == 1) || (sign2 == 1 && sign3 == -1)) ? -1 : 1;
//         }
//
//         std::vector<uint64_t> path_pos2 = (sign == 1) ? paths_pos2[j] : paths_neg2[j];
//         std::vector<uint64_t> path_neg2 = (sign == 1) ? paths_neg2[j] : paths_pos2[j];
//
//         Method2(i,j,total_paths,path_pos1,path_neg1,path_pos2,path_neg2,joined_paths_pos,joined_paths_neg,flipped_pivot_length,
//                 location,count,ValueTable,tid_localindicesq,tid_max_scores,p_caseorcontrol,vlen,vlen2,nCases,nControls,iterations,pathLength);
//
//       }
//
//     }
//
//     // Update the values in max_scores given the values of local_max_scores
//     for(int j = 0; j < nthreads; j++){
//       for(int k = 0; k < iterations; k++)
//         if(max_scores[j][k] < local_max_scores[j][k]) max_scores[j][k] = local_max_scores[j][k];
//     }
//
//   }
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
//   for(unsigned int j = 0; j < LocalIndicesQ.size(); j++){
//     IndicesScoresQueue &localqueue = LocalIndicesQ[j];
//     while(localqueue.size() != 0) {
//       indicesQ.push(localqueue.top());
//       localqueue.pop();
//       indicesQ.pop();
//     }
//   }
//
//   // Create a vector for the indices of the k paths with the highest scores
//   IntegerMatrix ids(indicesQ.size(),2);
//   NumericVector scores(indicesQ.size());
//   int j = 0;
//   while(indicesQ.size() != 0){
//     ScoreIndices temp = indicesQ.top();
//     scores[j] = temp.first;
//     ids(j,0) = temp.second.first + 1; // Add 1 because it will be used as an index in the R code
//     ids(j,1) = temp.second.second + 1;
//     indicesQ.pop();
//     j++;
//   }
//
//   return List::create(Named("scores") = scores,
//                       Named("ids") = ids,
//                       Named("TestScores") = permutedscores);
//
//
// }



// // [[Rcpp::export]]
// List ProcessPaths2(IntegerVector srcuids1, IntegerVector trguids1, List uids_CountLoc1, IntegerVector joining_gene_sign1,
//                   IntegerVector srcuids1_2, IntegerVector trguids1_2, List uids_CountLoc1_2, IntegerVector joining_gene_sign1_2,
//                   IntegerVector srcuids2, IntegerVector trguids2, List uids_CountLoc2, IntegerVector joining_gene_sign2,
//                   IntegerVector srcuids3, IntegerVector trguids3, List uids_CountLoc3, IntegerVector joining_gene_sign3,
//                   IntegerVector srcuids4, IntegerVector trguids4, List uids_CountLoc4, IntegerVector joining_gene_sign4,
//                   IntegerVector srcuids5, IntegerVector trguids5, List uids_CountLoc5, IntegerVector joining_gene_sign5,
//                   IntegerVector data_inds1, IntegerVector data_inds1_2, IntegerVector data_inds2, IntegerVector data_inds3,
//                   IntegerMatrix data, IntegerMatrix data2, NumericMatrix ValueTable, int nCases, int nControls, int K,
//                   int iterations, IntegerMatrix CaseORControl, int method, int pathLength, int nthreads){
//
//   int vlen = (int) ceil(nCases/64.0);
//   int vlen2 = (int) ceil(nControls/64.0);
//
//   // Copy ValueTable into an aligned array
//   double ValueTable_scores[ValueTable.nrow()*ValueTable.ncol()] ALIGN;
//   copyValueTable2(ValueTable_scores, ValueTable);
//
//   // double ValueTable_scores_flipped[ValueTable.nrow()*ValueTable.ncol()] ALIGN;
//   // copyValueTable2(ValueTable_scores_flipped, ValueTable);
//
//   // // Parse CaseORControl matrix into a single aligned vector
//   // uint64_t CaseORControl2[iterations*(vlen+vlen2)] ALIGN;
//   // parseCaseORControl2(CaseORControl2, CaseORControl, nCases, nControls, vlen, vlen2);
//   //
//   // // Get the parsed data and store it into an aligned array
//   // uint64_t parsed_data[data.nrow()*(vlen + vlen2)] ALIGN;
//   // parseData2(parsed_data, data, nCases, nControls, vlen, vlen2);
//   //
//   // uint64_t parsed_data2[data2.nrow()*(vlen + vlen2)] ALIGN;
//   // parseData2(parsed_data2, data2, nCases, nControls, vlen, vlen2);
//
//
//   // // Process paths of length 1
//   // List lst1;
//   // int total_paths = getTotalPaths(trguids1, uids_CountLoc1);
//   // std::vector<std::vector<uint64_t> > paths_pos1 = getZeroMatrix(total_paths, vlen + vlen2);
//   // std::vector<std::vector<uint64_t> > paths_neg1;
//   // if(method == 2) paths_neg1 = getZeroMatrix(total_paths, vlen + vlen2);
//   // std::vector<std::vector<uint64_t> > paths_pos0 = getZeroMatrix(data_inds1.size(), vlen + vlen2);
//   // std::vector<std::vector<uint64_t> > paths_neg0;
//   // if(method == 2) paths_neg0 = getZeroMatrix(data_inds1.size(), vlen + vlen2);
//   // std::vector<std::vector<uint64_t> > paths_data_pos = matchData(parsed_data, data_inds1);
//   // std::vector<std::vector<uint64_t> > paths_data_neg;
//   // if(method == 2) paths_data_neg = getZeroMatrix(data_inds1.size(), vlen + vlen2);
//   // int nnodes1 = unique(srcuids1).size();
//   //
//   // std::cout << nnodes1 << " source nodes for paths of length " << 1 << " and their permutations will be processed!" << std::endl;
//   //
//   // lst1 = JoinPaths(paths_pos0, paths_neg0, paths_data_pos, paths_data_neg, paths_pos1, paths_neg1,
//   //                  srcuids1, trguids1, uids_CountLoc1, joining_gene_sign1, ValueTable2, CaseORControl22,
//   //                  nCases, nControls, K, iterations, method, 1, nthreads);
//   //
//   //
//   // List lst1_2;
//   // total_paths = getTotalPaths(trguids1_2, uids_CountLoc1_2);
//   // std::vector<std::vector<uint64_t> > paths_pos1_2 = getZeroMatrix(total_paths, vlen + vlen2);
//   // std::vector<std::vector<uint64_t> > paths_neg1_2;
//   // if(method == 2) paths_neg1_2 = getZeroMatrix(total_paths, vlen + vlen2);
//   // std::vector<std::vector<uint64_t> > paths_pos0_2 = getZeroMatrix(data_inds1_2.size(), vlen + vlen2);
//   // std::vector<std::vector<uint64_t> > paths_neg0_2;
//   // if(method == 2) paths_neg0_2 = getZeroMatrix(data_inds1_2.size(), vlen + vlen2);
//   // paths_data_pos = matchData(parsed_data2, data_inds1_2);
//   // if(method == 2) paths_data_neg = getZeroMatrix(data_inds1_2.size(), vlen + vlen2);
//   // nnodes1 = unique(srcuids1_2).size();
//   //
//   // std::cout << nnodes1 << " source nodes for paths of length " << 1 << " and their permutations will be processed!" << std::endl;
//   //
//   // lst1_2 = JoinPaths(paths_pos0_2, paths_neg0_2, paths_data_pos, paths_data_neg, paths_pos1_2, paths_neg1_2,
//   //                    srcuids1_2, trguids1_2, uids_CountLoc1_2, joining_gene_sign1_2, ValueTable2, CaseORControl22,
//   //                    nCases, nControls, K, iterations, method, 1, nthreads);
//   //
//   //
//   // // Process Paths of length 2
//   // List lst2;
//   // total_paths = getTotalPaths(trguids2, uids_CountLoc2);
//   // std::vector<std::vector<uint64_t> > paths_pos2 = getZeroMatrix(total_paths, vlen + vlen2);
//   // std::vector<std::vector<uint64_t> > paths_neg2;
//   // if(method == 2) paths_neg2 = getZeroMatrix(total_paths, vlen + vlen2);
//   // paths_data_pos = matchData(parsed_data, data_inds2);
//   // if(method == 2) paths_data_neg = getZeroMatrix(data_inds2.size(), vlen + vlen2);
//   // int nnodes2 = unique(srcuids2).size();
//   //
//   // std::cout << nnodes2 << " source nodes for paths of length " << 2 << " and their permutations will be processed!" << std::endl;
//   //
//   // lst2 = JoinPaths(paths_pos1, paths_neg1, paths_data_pos, paths_data_neg, paths_pos2, paths_neg2,
//   //                  srcuids2, trguids2, uids_CountLoc2, joining_gene_sign2, ValueTable2,
//   //                  CaseORControl22, nCases, nControls, K, iterations, method, 2, nthreads);
//   //
//   // // Process Paths of length 3
//   // List lst3;
//   // total_paths = getTotalPaths(trguids3, uids_CountLoc3);
//   // std::vector<std::vector<uint64_t> > paths_pos3 = getZeroMatrix(total_paths, vlen + vlen2);
//   // std::vector<std::vector<uint64_t> > paths_neg3;
//   // if(method == 2) paths_neg3 = getZeroMatrix(total_paths, vlen + vlen2);
//   // paths_data_pos = matchData(parsed_data, data_inds3);
//   // if(method == 2) paths_data_neg = getZeroMatrix(data_inds3.size(), vlen + vlen2);
//   // int nnodes3 = unique(srcuids3).size();
//   //
//   // std::cout << nnodes3 << " source nodes for paths of length " << 3 << " and their permutations will be processed!" << std::endl;
//   //
//   // lst3 = JoinPaths(paths_pos2, paths_neg2, paths_data_pos, paths_data_neg, paths_pos3, paths_neg3,
//   //                  srcuids3, trguids3, uids_CountLoc3, joining_gene_sign3, ValueTable2,
//   //                  CaseORControl22, nCases, nControls, K, iterations, method, 3, nthreads);
//   //
//   // // Process Paths of length 4
//   // List lst4;
//   // std::vector<std::vector<uint64_t> > paths_pos4;
//   // std::vector<std::vector<uint64_t> > paths_neg4;
//   // int nnodes4 = unique(srcuids4).size();
//   //
//   // std::cout << nnodes4 << " source nodes for paths of length " << 4 << " and their permutations will be processed!" << std::endl;
//   //
//   // lst4 = JoinPaths(paths_pos3, paths_neg3, paths_pos2, paths_neg2, paths_pos4, paths_neg4,
//   //                  srcuids4, trguids4, uids_CountLoc4, joining_gene_sign4, ValueTable2,
//   //                  CaseORControl22, nCases, nControls, K, iterations, method, 4, nthreads);
//   //
//   // // Process Paths of length 5
//   // List lst5;
//   // std::vector<std::vector<uint64_t> > paths_pos5;
//   // std::vector<std::vector<uint64_t> > paths_neg5;
//   // int nnodes5 = unique(srcuids5).size();
//   //
//   // std::cout << nnodes5 << " source nodes for paths of length " << 5 << " and their permutations will be processed!" << std::endl;
//   //
//   // lst5 = JoinPaths(paths_pos3, paths_neg3, paths_pos3, paths_neg3, paths_pos5, paths_neg5,
//   //                  srcuids5, trguids5, uids_CountLoc5, joining_gene_sign5, ValueTable2,
//   //                  CaseORControl22, nCases, nControls, K, iterations, method, 5, nthreads);
//   //
//   // return List::create(Named("lst1") = lst1_2,
//   //                     Named("lst2") = lst2,
//   //                     Named("lst3") = lst3,
//   //                     Named("lst4") = lst4,
//   //                     Named("lst5") = lst5);
//
//   return List::create();
//
// }

