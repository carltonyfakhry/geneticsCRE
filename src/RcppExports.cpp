// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "geneticsCRE_types.h"
#include <Rcpp.h>

using namespace Rcpp;

// parsePaths
void parsePaths(IntegerMatrix data, int nCases, int nControls, std::string file_path);
RcppExport SEXP geneticsCRE_parsePaths(SEXP dataSEXP, SEXP nCasesSEXP, SEXP nControlsSEXP, SEXP file_pathSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type nCases(nCasesSEXP);
    Rcpp::traits::input_parameter< int >::type nControls(nControlsSEXP);
    Rcpp::traits::input_parameter< std::string >::type file_path(file_pathSEXP);
    parsePaths(data, nCases, nControls, file_path);
    return R_NilValue;
END_RCPP
}
// createPathSet
XPtr<paths_type> createPathSet(IntegerMatrix data, int num_cases, int num_controls, std::string method);
RcppExport SEXP geneticsCRE_createPathSet(SEXP dataSEXP, SEXP num_casesSEXP, SEXP num_controlsSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type num_cases(num_casesSEXP);
    Rcpp::traits::input_parameter< int >::type num_controls(num_controlsSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(createPathSet(data, num_cases, num_controls, method));
    return rcpp_result_gen;
END_RCPP
}
// getMatchingList
List getMatchingList(IntegerVector uids, IntegerVector counts, IntegerVector location);
RcppExport SEXP geneticsCRE_getMatchingList(SEXP uidsSEXP, SEXP countsSEXP, SEXP locationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type uids(uidsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type location(locationSEXP);
    rcpp_result_gen = Rcpp::wrap(getMatchingList(uids, counts, location));
    return rcpp_result_gen;
END_RCPP
}
// JoinIndices
List JoinIndices(IntegerVector r_src_uids, IntegerVector r_trg_uids, List uids_CountLoc, IntegerVector r_join_gene_signs, NumericMatrix r_value_table, int nCases, int nControls, int K, int iterations, IntegerMatrix CaseORControl, std::string method, int pathLength, int nthreads, SEXP xp_paths0, SEXP xp_paths1, bool keep_joined, std::string pos_path1, std::string neg_path1, std::string conflict_path1, std::string pos_path2, std::string neg_path2, std::string conflict_path2, std::string dest_path_pos, std::string dest_path_neg, std::string dest_path_conflict);
RcppExport SEXP geneticsCRE_JoinIndices(SEXP r_src_uidsSEXP, SEXP r_trg_uidsSEXP, SEXP uids_CountLocSEXP, SEXP r_join_gene_signsSEXP, SEXP r_value_tableSEXP, SEXP nCasesSEXP, SEXP nControlsSEXP, SEXP KSEXP, SEXP iterationsSEXP, SEXP CaseORControlSEXP, SEXP methodSEXP, SEXP pathLengthSEXP, SEXP nthreadsSEXP, SEXP xp_paths0SEXP, SEXP xp_paths1SEXP, SEXP keep_joinedSEXP, SEXP pos_path1SEXP, SEXP neg_path1SEXP, SEXP conflict_path1SEXP, SEXP pos_path2SEXP, SEXP neg_path2SEXP, SEXP conflict_path2SEXP, SEXP dest_path_posSEXP, SEXP dest_path_negSEXP, SEXP dest_path_conflictSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type r_src_uids(r_src_uidsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type r_trg_uids(r_trg_uidsSEXP);
    Rcpp::traits::input_parameter< List >::type uids_CountLoc(uids_CountLocSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type r_join_gene_signs(r_join_gene_signsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type r_value_table(r_value_tableSEXP);
    Rcpp::traits::input_parameter< int >::type nCases(nCasesSEXP);
    Rcpp::traits::input_parameter< int >::type nControls(nControlsSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type CaseORControl(CaseORControlSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type pathLength(pathLengthSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xp_paths0(xp_paths0SEXP);
    Rcpp::traits::input_parameter< SEXP >::type xp_paths1(xp_paths1SEXP);
    Rcpp::traits::input_parameter< bool >::type keep_joined(keep_joinedSEXP);
    Rcpp::traits::input_parameter< std::string >::type pos_path1(pos_path1SEXP);
    Rcpp::traits::input_parameter< std::string >::type neg_path1(neg_path1SEXP);
    Rcpp::traits::input_parameter< std::string >::type conflict_path1(conflict_path1SEXP);
    Rcpp::traits::input_parameter< std::string >::type pos_path2(pos_path2SEXP);
    Rcpp::traits::input_parameter< std::string >::type neg_path2(neg_path2SEXP);
    Rcpp::traits::input_parameter< std::string >::type conflict_path2(conflict_path2SEXP);
    Rcpp::traits::input_parameter< std::string >::type dest_path_pos(dest_path_posSEXP);
    Rcpp::traits::input_parameter< std::string >::type dest_path_neg(dest_path_negSEXP);
    Rcpp::traits::input_parameter< std::string >::type dest_path_conflict(dest_path_conflictSEXP);
    rcpp_result_gen = Rcpp::wrap(JoinIndices(r_src_uids, r_trg_uids, uids_CountLoc, r_join_gene_signs, r_value_table, nCases, nControls, K, iterations, CaseORControl, method, pathLength, nthreads, xp_paths0, xp_paths1, keep_joined, pos_path1, neg_path1, conflict_path1, pos_path2, neg_path2, conflict_path2, dest_path_pos, dest_path_neg, dest_path_conflict));
    return rcpp_result_gen;
END_RCPP
}
