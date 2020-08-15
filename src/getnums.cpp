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

void getnums(int nsam, int nreps, NumericVector t, IntegerVector I_rcpp, NumericVector migr_rcpp, NumericVector migration, NumericMatrix en,
    NumericMatrix ej, NumericMatrix es, NumericVector seeds, char *jtree){
      
    int i, m, n, j, r, npop, pop, pop2, tempArg;
    struct devent *ptemp , *pt ;

    int nseeds;
    int seeds_c[3];
    int commandlineseed( int * ) ;
    void free_eventlist( struct devent *pt, int npop ), addtoelist( struct devent *pt, struct devent *elist );
    NumericVector rcpp_row_en, rcpp_row_ej, rcpp_row_es;
    
    /* Pars is the global variable structure that holds the call's variables as well as creates the past events. Check ms.h for the details of the structure. */
    pars.cp.nsamin = nsam;
    pars.cp.nsam = pars.cp.nsamin;
    if( pars.cp.nsam <= 0 ) { fprintf(stderr,"Error: nsam <= 0. \n"); }//usage();}
    pars.output_precision =  4;
    pars.cp.r = pars.mp.theta =  pars.cp.f = 0.0;
    pars.cp.track_len = 1;
    pars.cp.npop = npop = 1;
    pars.cp.mig_mat = (double **)malloc( (unsigned) sizeof( double *) );
    pars.cp.mig_mat[0] = (double *)malloc( (unsigned)sizeof(double ));
    pars.cp.mig_mat[0][0] =  0.0 ;
    pars.mp.segsitesin = 0;
    pars.mp.mfreq = 1;
    pars.cp.config = (int *) malloc( (unsigned)(( pars.cp.npop +1 ) *sizeof( int)) );
    (pars.cp.config)[0] = pars.cp.nsamin ;
    pars.cp.size= (double *) malloc( (unsigned)( pars.cp.npop *sizeof( double )) );
    (pars.cp.size)[0] = 1.0  ;
    pars.cp.alphag = (double *) malloc( (unsigned)(( pars.cp.npop ) *sizeof( double )) );
    (pars.cp.alphag)[0] = 0.0  ;
    pars.cp.nsites = 2;
    
    pars.cp.deventlist = NULL;

    /*End of intitialization*/
    
    /* Each case below exists to populate the global c structure pars with the rcpp getnums variables. Since it is possible to not include all of the cases in the same run, the code checks to make sure the variables are non-zero or have the correct structure in order to be considered.  */
    
    
    for(int i = 0; i < 3; i++){ // the seed (as three numbers) are passed in and changed to the commandlineseed parameter
      seeds_c[i] = seeds[i];
    }
    nseeds = commandlineseed(seeds_c);
    
    /* Case t */
    if(t[0] != 0){
        pars.mp.theta = t[0]; //populating global parameter theta
    }
    else{ fprintf(stderr,"Theta (t) cannot be 0");}
    
    // Case I
    if (I_rcpp.length() > 1){
        pars.cp.npop = I_rcpp[0]; /* I splits into subpopulations, with the first number being the amount of subpops */
        pars.cp.config = (int *) realloc(pars.cp.config, (unsigned)( pars.cp.npop*sizeof( int)));
        npop = pars.cp.npop;
        
        for(m = 1; m < I_rcpp.length(); m++){
            pars.cp.config[m-1] = I_rcpp[m];
        }
        
        /* Creating the migration matrix structure, and populating it according to the I vector as well as the migration parameter (migr). This only changes if the migration matrix variable is used */
        pars.cp.mig_mat = (double **)realloc(pars.cp.mig_mat, (unsigned)(pars.cp.npop*sizeof(double *) )) ;
        pars.cp.mig_mat[0] = (double *)realloc(pars.cp.mig_mat[0], (unsigned)( pars.cp.npop*sizeof(double)));
        
        for(i=1; i<pars.cp.npop; i++){
            pars.cp.mig_mat[i] = (double *)malloc( (unsigned)( pars.cp.npop*sizeof(double)));
        }
        
        pars.cp.size = (double *)realloc( pars.cp.size, (unsigned)( pars.cp.npop*sizeof( double )));
        pars.cp.alphag = (double *) realloc( pars.cp.alphag, (unsigned)( pars.cp.npop*sizeof( double )));
        
        for( i=1; i< pars.cp.npop ; i++) {
            (pars.cp.size)[i] = (pars.cp.size)[0]  ;
            (pars.cp.alphag)[i] = (pars.cp.alphag)[0] ;
        }

        for( i=0; i<pars.cp.npop; i++){
            for( j=0; j<pars.cp.npop; j++) {
                pars.cp.mig_mat[i][j] = migr_rcpp[0]/(pars.cp.npop - 1) ;
            }
        }

        for( i=0; i< pars.cp.npop; i++){
                pars.cp.mig_mat[i][i] = migr_rcpp[0] ;
        }
    }
    
    if (migration.length() > 1){ /*repopulates the migration matrix based off of the vector*/
        double migmat_array[migration.length()];
        for(m = 0; m < migration.length(); m++){
            migmat_array[m] = migration[m];
        }
        tempArg = -1;
        for( pop = 0; pop <npop; pop++){
            for( pop2 = 0; pop2 <npop; pop2++){
                tempArg = tempArg + 1;
                pars.cp.mig_mat[pop][pop2]= migmat_array[tempArg];
            }
        }
        for( pop = 0; pop < npop; pop++) {
            pars.cp.mig_mat[pop][pop] = 0.0 ;
            for( pop2 = 0; pop2 < npop; pop2++){
                if( pop2 != pop ) pars.cp.mig_mat[pop][pop] += pars.cp.mig_mat[pop][pop2] ;
                }
            }
    }
    
    /* The e cases follow the same structure - going through each row of the matrix to turn them int an event structure and then populating the last two numbers of each row into their specific pars variable. Most of the code is the same for every e case; the unique code of each case comes in the last few lines */
    
    if (en.ncol() == 3){
        double row_en[3];
        NumericVector numrow = en.nrow();
        for(n = 0; n < numrow.length(); n++){ /* loops through the amount of en instances there are - the number of rows in the matrix */
            NumericVector rcpp_row_en = en(n,_);
            for(r = 0; r < 3; r++){
                row_en[r] = rcpp_row_en[r]; /*creating c variable for each row*/
            }
            pt = (struct devent *)malloc( sizeof( struct devent) ) ; /*creating an event*/
            pt->detype = 'n';
            
            pt->time = row_en[0]; /* time is the first number of every row*/
            pt->nextde = NULL ;

            /*adds the event to the event list at the correct time: */
            if( pars.cp.deventlist == NULL ){
                pars.cp.deventlist = pt ;
            }
            else if ( pt->time < pars.cp.deventlist->time ) {
                ptemp = pars.cp.deventlist ;
                pars.cp.deventlist = pt ;
                pt->nextde = ptemp ;
            }
            else {
                addtoelist( pt, pars.cp.deventlist ) ;
            }
            pt->popi =  row_en[1] -1 ; /*popi = subpopulation, the second number in the row */
            pt->paramv = row_en[2] ; /*paramv = size of subpopulation, third number in the row*/
        }
    }

    if (ej.ncol() == 3){
        double row_ej[3];
        NumericVector numrow_ej = ej.nrow();
        for(j = 0; j < numrow_ej.length(); j++){
            Rcpp::NumericVector rcpp_row_ej = ej(j,_);
            for(r = 0; r < 3; r++){
                row_ej[r] = rcpp_row_ej[r];
                }
            pt = (struct devent *)malloc( sizeof( struct devent) ) ;
            pt->detype = 'j';
            
            pt->time = row_ej[0];
            pt->nextde = NULL ;
            if( pars.cp.deventlist == NULL ){
                pars.cp.deventlist = pt ;
		        }
            else if ( pt->time < pars.cp.deventlist->time ) {
                ptemp = pars.cp.deventlist ;
                pars.cp.deventlist = pt ;
                pt->nextde = ptemp ;
            }
            else{ addtoelist( pt, pars.cp.deventlist ) ;}
            pt ->popi = row_ej[1] - 1; /*subpopulation 1*/
            pt ->popj = row_ej[2] - 1; /*subpopulation 2*/
            }
    }
      
      if (es.ncol() == 3){
          double row_es[3]; //vector for each row
          NumericVector numrow_es = es.nrow();
        for(j = 0; j < numrow_es.length(); j++) { //loop to go down each row
              Rcpp::NumericVector rcpp_row_es = es(j,_);
              for(r = 0; r < 3; r++){
                  row_es[r] = rcpp_row_es[r];
                  }
              pt = (struct devent *)malloc( sizeof( struct devent) ) ;
              pt->detype = 's';
              
              pt->time = row_es[0];
              pt->nextde = NULL ;
              if( pars.cp.deventlist == NULL ){
                  pars.cp.deventlist = pt ;
                  }
              else if ( pt->time < pars.cp.deventlist->time ) {
                  ptemp = pars.cp.deventlist ;
                  pars.cp.deventlist = pt ;
                  pt->nextde = ptemp ;
              }
              else{ addtoelist( pt, pars.cp.deventlist ) ; }
              
              pt->popi = row_es[1] - 1 ; /* subpopulation */
              pt->paramv = row_es[2] - 1 ; /*paramv = probability*/
        }
      }
}

