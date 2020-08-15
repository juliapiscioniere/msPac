/***** ms.c     ************************************************

 This is an R package that generates genetic trees. Modified from original ms software that also produced samples of gametes.
 
 See the complete testing document in order to see the package work with multiple different scenarios. Here is a brief overview of each parameter:

           npop:  Number of subpopulations which make up the total population - first number in I vector
 
           ni (rest of I vector):  the sample size from the i th subpopulation (all must be specified.) The output will have the gametes in order such that the first n1 gametes are from the first island, the next n2 are
           from the second island, etc.
 
           nsites: number of sites between which recombination can occur.
 
           theta: 4No times the neutral mutation rate
 
           track_len:  mean length of conversion track in units of sites.  The
		       total number of sites is nsites, specified with the -r option.
 
           mig_rate (migr) : migration rate: the fraction of each subpop made up of migrants times 4No.

	Note:  In the above definition, No is the total diploid population if
		npop is one, otherwise, No is the diploid population size of each
		subpopulation.


 Modified by Julia Piscioniere
 Code base from ms software written by Richard Hudson.
 Link to original source code: https://uchicago.app.box.com/s/13e5uf13tikfjm7e1il1eujitlsjdx13
 Written in collaboration with Dr. Allan Strand.
***************************************************************************/

#include <Rcpp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "getnums.h"

extern "C" {
 #include "ms.h"
}

using namespace Rcpp;

struct segl {
	int beg;
	struct node *ptree;
	int next;
	};

double *posit ;
double *agevec ;
double alleleage ;
double segfac ;
int nseeds ;

extern "C" {struct params pars;}
extern "C" {unsigned maxsites;}

// [[Rcpp::export]]
Rcpp::CharacterVector ms_main(NumericVector nsam, NumericVector nreps, NumericVector t, IntegerVector I_rcpp, NumericVector migr_rcpp, NumericVector migration, NumericMatrix en, NumericMatrix ej, NumericMatrix es)
{
    
    char *tree ;
   // int i, segsites;
    char **list, **cmatrix(int gsam, int len);
    double probss, tmrca, ttot;
    Rcpp::CharacterVector rcpp_tree (nreps[0]); /*initializing rcpp tree vector to be the same length as how many trees being produced*/
    Rcpp::NumericVector seeds(3);

    maxsites = SITESINC;

    for (int rep=0; rep<rcpp_tree.length();rep++) /*loop that uses nrep to produce that many trees - runs ms that many times*/
      {
        seeds = floor(runif(3)*1000000); /*seeds is a set of 3 random numbers*/
        
        tree = (char *)malloc(500 * nsam[0] * sizeof(char)); /*allocating memory to the one C tree being produced in this iteration of ms*/
        strcpy(&tree[0],"\0"); //string init

        getnums(nsam[0] , 1, t, I_rcpp, migr_rcpp, migration, en, ej, es, seeds, tree); /*getnums takes parameters and stores them in a global structure called pars*/
          
        list = cmatrix(pars.cp.nsam,maxsites+1);
        posit = (double *)malloc( (unsigned)( maxsites*sizeof( double)));
        agevec = (double *)malloc( (unsigned)( maxsites*sizeof(double)));

        gensam( list, &probss, &tmrca, &ttot, tree); /*gensam calls segtre_mig function as well as print tree functions in order to create the tree*/

        char *p = tree;
        rcpp_tree[rep] = p; /*adding this tree as another tree in the rcpp vector  */
        free(tree);

      }
    
    return rcpp_tree; /*return vector of trees */
}
