#ifndef GCRE_H
#define GCRE_H

using namespace std;

vector<vector<uint64_t>> parseCaseORControl(IntegerMatrix CaseORControl, int nCases, int nControls);
void parsePaths(IntegerMatrix data, int nCases, int nControls, string file_path);
void StorePaths(vector<vector<uint64_t>> &paths, string file_path);
vector<vector<uint64_t>> readPaths(string file_path);
int getTotalPaths(IntegerVector trguids, List uids_CountLoc);
vector<vector<uint64_t>> getZeroMatrix(int dim1, int dim2);
int getTotalCountsCountLoc(List uids_CountLoc);

List JoinIndicesMethod1(IntegerVector srcuid, IntegerVector trguids2, List uids_CountLoc, IntegerVector joining_gene_sign,
  NumericMatrix ValueTable, int nCases, int nControls, int K,
  int iterations, IntegerMatrix CaseORControl, int pathLength, int nthreads, string pos_path1,
  string pos_path2, string dest_path_pos);

List JoinIndicesMethod2(IntegerVector srcuid, IntegerVector trguids2, List uids_CountLoc, IntegerVector joining_gene_sign,
  NumericMatrix ValueTable, int nCases, int nControls, int K,
  int iterations, IntegerMatrix CaseORControl, int pathLength, int nthreads, string pos_path1,
  string neg_path1, string conflict_path1, string pos_path2, string neg_path2, string conflict_path2,
  string dest_path_pos, string dest_path_neg, string dest_path_conflict);

#endif
