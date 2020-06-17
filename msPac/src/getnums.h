//.h file for getnums.cpp
#include <Rcpp.h>

using namespace Rcpp;
void getnums(int *howmany, NumericVector nsam, NumericVector nreps, NumericVector t,
        NumericVector variable_list_rcpp, IntegerVector I_rcpp, NumericVector migr_rcpp, NumericVector migration, NumericMatrix en, NumericMatrix ej, NumericVector seeds); 

//caseen(double *en, struct params pars);
//void caseej(double *ej);
//void caseI(int *I, double migr);
//casema(double *migmat_array, struct params pars);
