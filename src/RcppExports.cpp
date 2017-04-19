// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getRels3
Rcpp::DataFrame getRels3(Rcpp::IntegerVector srcuid, Rcpp::IntegerVector trguid, Rcpp::IntegerVector sign, Rcpp::List uids_CountLoc3);
RcppExport SEXP geneticsCRE_getRels3(SEXP srcuidSEXP, SEXP trguidSEXP, SEXP signSEXP, SEXP uids_CountLoc3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type srcuid(srcuidSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type trguid(trguidSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type sign(signSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type uids_CountLoc3(uids_CountLoc3SEXP);
    rcpp_result_gen = Rcpp::wrap(getRels3(srcuid, trguid, sign, uids_CountLoc3));
    return rcpp_result_gen;
END_RCPP
}
// getMatchingList
Rcpp::List getMatchingList(Rcpp::IntegerVector uids, Rcpp::IntegerVector counts, Rcpp::IntegerVector location);
RcppExport SEXP geneticsCRE_getMatchingList(SEXP uidsSEXP, SEXP countsSEXP, SEXP locationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type uids(uidsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type location(locationSEXP);
    rcpp_result_gen = Rcpp::wrap(getMatchingList(uids, counts, location));
    return rcpp_result_gen;
END_RCPP
}
// ProcessPaths
Rcpp::List ProcessPaths(Rcpp::IntegerVector r_src_uids1, Rcpp::IntegerVector r_trg_uids1, Rcpp::List r_count_locs1, Rcpp::IntegerVector r_signs1, Rcpp::IntegerVector r_src_uids1_2, Rcpp::IntegerVector r_trg_uids1_2, Rcpp::List r_count_locs1_2, Rcpp::IntegerVector r_signs1_2, Rcpp::IntegerVector r_src_uids2, Rcpp::IntegerVector r_trg_uids2, Rcpp::List r_count_locs2, Rcpp::IntegerVector r_signs2, Rcpp::IntegerVector r_src_uids3, Rcpp::IntegerVector r_trg_uids3, Rcpp::List r_count_locs3, Rcpp::IntegerVector r_signs3, Rcpp::IntegerVector r_src_uids4, Rcpp::IntegerVector r_trg_uids4, Rcpp::List r_count_locs4, Rcpp::IntegerVector r_signs4, Rcpp::IntegerVector r_src_uids5, Rcpp::IntegerVector r_trg_uids5, Rcpp::List r_count_locs5, Rcpp::IntegerVector r_signs5, Rcpp::IntegerVector r_data_inds1, Rcpp::IntegerVector r_data_inds1_2, Rcpp::IntegerVector r_data_inds2, Rcpp::IntegerVector r_data_inds3, Rcpp::IntegerMatrix r_data1, Rcpp::IntegerMatrix r_data2, Rcpp::NumericMatrix r_value_table, int num_cases, int num_ctrls, int top_k, int iterations, Rcpp::IntegerMatrix r_perm_cases, std::string method, int path_length, int nthreads);
RcppExport SEXP geneticsCRE_ProcessPaths(SEXP r_src_uids1SEXP, SEXP r_trg_uids1SEXP, SEXP r_count_locs1SEXP, SEXP r_signs1SEXP, SEXP r_src_uids1_2SEXP, SEXP r_trg_uids1_2SEXP, SEXP r_count_locs1_2SEXP, SEXP r_signs1_2SEXP, SEXP r_src_uids2SEXP, SEXP r_trg_uids2SEXP, SEXP r_count_locs2SEXP, SEXP r_signs2SEXP, SEXP r_src_uids3SEXP, SEXP r_trg_uids3SEXP, SEXP r_count_locs3SEXP, SEXP r_signs3SEXP, SEXP r_src_uids4SEXP, SEXP r_trg_uids4SEXP, SEXP r_count_locs4SEXP, SEXP r_signs4SEXP, SEXP r_src_uids5SEXP, SEXP r_trg_uids5SEXP, SEXP r_count_locs5SEXP, SEXP r_signs5SEXP, SEXP r_data_inds1SEXP, SEXP r_data_inds1_2SEXP, SEXP r_data_inds2SEXP, SEXP r_data_inds3SEXP, SEXP r_data1SEXP, SEXP r_data2SEXP, SEXP r_value_tableSEXP, SEXP num_casesSEXP, SEXP num_ctrlsSEXP, SEXP top_kSEXP, SEXP iterationsSEXP, SEXP r_perm_casesSEXP, SEXP methodSEXP, SEXP path_lengthSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_src_uids1(r_src_uids1SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_trg_uids1(r_trg_uids1SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type r_count_locs1(r_count_locs1SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_signs1(r_signs1SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_src_uids1_2(r_src_uids1_2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_trg_uids1_2(r_trg_uids1_2SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type r_count_locs1_2(r_count_locs1_2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_signs1_2(r_signs1_2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_src_uids2(r_src_uids2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_trg_uids2(r_trg_uids2SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type r_count_locs2(r_count_locs2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_signs2(r_signs2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_src_uids3(r_src_uids3SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_trg_uids3(r_trg_uids3SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type r_count_locs3(r_count_locs3SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_signs3(r_signs3SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_src_uids4(r_src_uids4SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_trg_uids4(r_trg_uids4SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type r_count_locs4(r_count_locs4SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_signs4(r_signs4SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_src_uids5(r_src_uids5SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_trg_uids5(r_trg_uids5SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type r_count_locs5(r_count_locs5SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_signs5(r_signs5SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_data_inds1(r_data_inds1SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_data_inds1_2(r_data_inds1_2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_data_inds2(r_data_inds2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type r_data_inds3(r_data_inds3SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type r_data1(r_data1SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type r_data2(r_data2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type r_value_table(r_value_tableSEXP);
    Rcpp::traits::input_parameter< int >::type num_cases(num_casesSEXP);
    Rcpp::traits::input_parameter< int >::type num_ctrls(num_ctrlsSEXP);
    Rcpp::traits::input_parameter< int >::type top_k(top_kSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type r_perm_cases(r_perm_casesSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type path_length(path_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(ProcessPaths(r_src_uids1, r_trg_uids1, r_count_locs1, r_signs1, r_src_uids1_2, r_trg_uids1_2, r_count_locs1_2, r_signs1_2, r_src_uids2, r_trg_uids2, r_count_locs2, r_signs2, r_src_uids3, r_trg_uids3, r_count_locs3, r_signs3, r_src_uids4, r_trg_uids4, r_count_locs4, r_signs4, r_src_uids5, r_trg_uids5, r_count_locs5, r_signs5, r_data_inds1, r_data_inds1_2, r_data_inds2, r_data_inds3, r_data1, r_data2, r_value_table, num_cases, num_ctrls, top_k, iterations, r_perm_cases, method, path_length, nthreads));
    return rcpp_result_gen;
END_RCPP
}
