// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// EvalStatSlow
double EvalStatSlow(NumericMatrix X, int s, int t, String optimization);
RcppExport SEXP _OCBS_EvalStatSlow(SEXP XSEXP, SEXP sSEXP, SEXP tSEXP, SEXP optimizationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< String >::type optimization(optimizationSEXP);
    rcpp_result_gen = Rcpp::wrap(EvalStatSlow(X, s, t, optimization));
    return rcpp_result_gen;
END_RCPP
}
// MaxStats
List MaxStats(NumericMatrix X, String optimization, String method, bool circular, int min_seg);
RcppExport SEXP _OCBS_MaxStats(SEXP XSEXP, SEXP optimizationSEXP, SEXP methodSEXP, SEXP circularSEXP, SEXP min_segSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< String >::type optimization(optimizationSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    Rcpp::traits::input_parameter< bool >::type circular(circularSEXP);
    Rcpp::traits::input_parameter< int >::type min_seg(min_segSEXP);
    rcpp_result_gen = Rcpp::wrap(MaxStats(X, optimization, method, circular, min_seg));
    return rcpp_result_gen;
END_RCPP
}
// PermTest
bool PermTest(NumericMatrix X, IntegerVector boundary, double cand_stat, String optimization, String method, double alpha, int nr_perms, bool circular, int min_seg);
RcppExport SEXP _OCBS_PermTest(SEXP XSEXP, SEXP boundarySEXP, SEXP cand_statSEXP, SEXP optimizationSEXP, SEXP methodSEXP, SEXP alphaSEXP, SEXP nr_permsSEXP, SEXP circularSEXP, SEXP min_segSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< double >::type cand_stat(cand_statSEXP);
    Rcpp::traits::input_parameter< String >::type optimization(optimizationSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type nr_perms(nr_permsSEXP);
    Rcpp::traits::input_parameter< bool >::type circular(circularSEXP);
    Rcpp::traits::input_parameter< int >::type min_seg(min_segSEXP);
    rcpp_result_gen = Rcpp::wrap(PermTest(X, boundary, cand_stat, optimization, method, alpha, nr_perms, circular, min_seg));
    return rcpp_result_gen;
END_RCPP
}
// StoppingBoundary
IntegerVector StoppingBoundary(int nr_perms, double alpha, double eta_star);
RcppExport SEXP _OCBS_StoppingBoundary(SEXP nr_permsSEXP, SEXP alphaSEXP, SEXP eta_starSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nr_perms(nr_permsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type eta_star(eta_starSEXP);
    rcpp_result_gen = Rcpp::wrap(StoppingBoundary(nr_perms, alpha, eta_star));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_OCBS_EvalStatSlow", (DL_FUNC) &_OCBS_EvalStatSlow, 4},
    {"_OCBS_MaxStats", (DL_FUNC) &_OCBS_MaxStats, 5},
    {"_OCBS_PermTest", (DL_FUNC) &_OCBS_PermTest, 9},
    {"_OCBS_StoppingBoundary", (DL_FUNC) &_OCBS_StoppingBoundary, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_OCBS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}