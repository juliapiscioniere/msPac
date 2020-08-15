//.h file for getnums.cpp
#include <Rcpp.h>

using namespace Rcpp;
void getnums(int nsam, int nreps, NumericVector t, IntegerVector I_rcpp, NumericVector migr_rcpp, NumericVector migration, NumericMatrix en, NumericMatrix ej, NumericMatrix es, NumericVector seeds, char *jtree);
