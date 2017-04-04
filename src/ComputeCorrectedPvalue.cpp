// #include "Utils.h"
//
//
//
// // [[Rcpp::export]]
// NumericVector computeCorrectedPvalues(IntegerMatrix path_data, IntegerVector sign,
//                                       NumericMatrix ValueTable, IntegerMatrix permutations,
//                                       int method, int pathLength){
//
//   NumericVector pvalues(pathLength*2);
//
//   for(int i = 0; i < path_data.nrow(); i++){
//
//     IntegerVector sub_path_pos1 = path_data.row(i);
//     IntegerVector sub_path_pos2(sub_path_pos1.size());
//
//     for(int j = i + 1; j < path_data.nrow(); j++){
//
//       IntegerVector temp = path_data.row(j);
//
//       for(int k = 0; k < path_data.nrow(); k++)
//         sub_path_pos2[k] |= temp[k];
//
//     }
//
//     int cases = 0;
//     int controls = 0;
//
//     for(int i = 0; i < sub_path_pos1.size(); i++){
//
//       if(sub_path_pos2[i] != sub_path_pos1[i] && sub_path_pos2[i] == 1){
//
//
//
//       }
//
//     }
//
//   }
//
//   return pvalues;
//
// }
