// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ms_main
Rcpp::CharacterVector ms_main(NumericVector nsam, NumericVector nreps, NumericVector t, IntegerVector I_rcpp, NumericVector migr_rcpp, NumericVector migration, NumericMatrix en, NumericMatrix ej, NumericMatrix es);
RcppExport SEXP _msPac_ms_main(SEXP nsamSEXP, SEXP nrepsSEXP, SEXP tSEXP, SEXP I_rcppSEXP, SEXP migr_rcppSEXP, SEXP migrationSEXP, SEXP enSEXP, SEXP ejSEXP, SEXP esSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nsam(nsamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nreps(nrepsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type I_rcpp(I_rcppSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type migr_rcpp(migr_rcppSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type migration(migrationSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type en(enSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ej(ejSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type es(esSEXP);
    rcpp_result_gen = Rcpp::wrap(ms_main(nsam, nreps, t, I_rcpp, migr_rcpp, migration, en, ej, es));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_msPac_ms_main", (DL_FUNC) &_msPac_ms_main, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_msPac(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
