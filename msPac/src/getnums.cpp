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
//initialization, the case I initialiZation, and then the migration inirialization

void getnums(int *phowmany, NumericVector nsam, NumericVector nreps, NumericVector t, NumericVector variable_list,
             IntegerVector I_rcpp, NumericVector migr_rcpp, NumericVector migration, NumericMatrix en, 
             NumericMatrix ej, NumericVector seeds){
    int i, m, n, j, r, npop, I_length, pop, pop2, tempArg;
    struct devent *ptemp , *pt ;

    int nseeds;
    int seeds_c[3]; //seeds vector will be 3 numbers
    int commandlineseed( int * ) ; //need this after i changed it?
    void free_eventlist( struct devent *pt, int npop ), addtoelist( struct devent *pt, struct devent *elist );
    NumericVector rcpp_row_en, rcpp_row_ej;
    
    pars.commandlineseedflag =  variable_list[0];
    pars.cp.nsamin = nsam[0]; //one number vectors 
    pars.cp.nsam = pars.cp.nsamin;
    if( pars.cp.nsam <= 0 ) { fprintf(stderr,"Error: nsam <= 0. \n"); usage();}
    *phowmany = nreps[0];
    if( *phowmany  <= 0 ) { fprintf(stderr,"Error: howmany <= 0. \n"); usage();}
    pars.output_precision =  variable_list[1];
    pars.cp.r = pars.mp.theta =  pars.cp.f = variable_list[2]; //should be 0.0
    pars.cp.track_len = variable_list[3] ;
    pars.cp.npop = npop = variable_list[4] ;
    pars.cp.mig_mat = (double **)malloc( (unsigned) sizeof( double *) );
    pars.cp.mig_mat[0] = (double *)malloc( (unsigned)sizeof(double ));
    pars.cp.mig_mat[0][0] =  0.0 ;
    pars.mp.segsitesin = variable_list[5] ;
    pars.mp.treeflag = variable_list[6] ;
    pars.mp.timeflag =  variable_list[7];
    pars.mp.ageflag =  variable_list[8];
    pars.mp.mfreq = variable_list[9];
    pars.cp.config = (int *) malloc( (unsigned)(( pars.cp.npop +1 ) *sizeof( int)) );
    (pars.cp.config)[0] = pars.cp.nsamin ;
    pars.cp.size= (double *) malloc( (unsigned)( pars.cp.npop *sizeof( double )) );
    (pars.cp.size)[0] = 1.0  ;
    pars.cp.alphag = (double *) malloc( (unsigned)(( pars.cp.npop ) *sizeof( double )) );
    (pars.cp.alphag)[0] = 0.0  ;
    pars.cp.nsites = variable_list[10] ; //*
    
    pars.cp.deventlist = NULL;

    /*End of intitialization*/
    
    /* Case t */
    if(t[0] != 0){ //making sure that case t is initialized correctly, then stores the number for t into the global variable
        pars.mp.theta = t[0];
    }
    
    for(int i = 0; i < 3; i++){ // the seed (as three numbers) are passed in and changed to the commandlineseed parameter
      seeds_c[i] = seeds[i];
    }
    nseeds = commandlineseed(seeds_c);
    
    /* Case T - Needs to be 1 to return tree */
    pars.mp.treeflag = 1 ;
    
    // Case I
    I_length = I_rcpp.length(); //              of the rcpp vector = how many numbers in array
    if (I_length > 1){ //if length is greater than 1, then the I vector is there
        pars.cp.npop = I_rcpp[0];
        pars.cp.config = (int *) realloc(pars.cp.config, (unsigned)( pars.cp.npop*sizeof( int)));
        npop = pars.cp.npop;
        
        for(m = 1; m < I_length; m++){
            pars.cp.config[m-1] = I_rcpp[m];
        }
        
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
                pars.cp.mig_mat[i][j] = migr_rcpp[0]/(pars.cp.npop - 1) ; //changed it so it uses migr from rcpp
            }
        }
        for( i=0; i< pars.cp.npop; i++){
                pars.cp.mig_mat[i][i] = migr_rcpp[0] ;
        }
    }
    
    if (migration.length() > 1){ //if migration vector is there
        double migmat_array[migration.length()];
        for(m = 0; m < migration.length(); m++){
            //Rprintf("migration[%d] = %f", m, migration[m]);
            //Rprintf("\n");
            migmat_array[m] = migration[m];
            //Rprintf("migmat_array[%d] = %f", m, migmat_array[m]);
            //Rprintf("\n");
        }
        tempArg = -1;
        for( pop = 0; pop <npop; pop++){ //going through the correct number, up to 3 (0, 1, 2) for a 3 x 3 matrix
            for( pop2 = 0; pop2 <npop; pop2++){
//                Rprintf("tempArg before = %d", tempArg);
//                Rprintf("\n");
                tempArg = tempArg + 1;
                pars.cp.mig_mat[pop][pop2]= migmat_array[tempArg];
                
//                Rprintf("First Loop - migmat_array[%d] = %f", tempArg, migmat_array[tempArg]);
//                Rprintf("\n");
//                Rprintf("First Loop - temp arg = %d -- pars.cp.migmat[%d][%d] = %f", tempArg, pop, pop2, pars.cp.mig_mat[pop][pop2]);
//                Rprintf("\n");
            }
        }
        for( pop = 0; pop < npop; pop++) {
            pars.cp.mig_mat[pop][pop] = 0.0 ;
            for( pop2 = 0; pop2 < npop; pop2++){
//                Rprintf("\n");
//                Rprintf("pars.cp.migmat[%d][%d] = %f", pop, pop2, pars.cp.mig_mat[pop][pop2]);
//                Rprintf("\n");
                if( pop2 != pop ) pars.cp.mig_mat[pop][pop] += pars.cp.mig_mat[pop][pop2] ;
                }
            }
//        Rprintf("\n");
//    for(i =0; i< pars.cp.npop; i++){
//        for(j = 0; j < pars.cp.npop; j++){
////            Rprintf("\n");
////            Rprintf("pars.cp.npop = %d", pars.cp.npop);
////            Rprintf("\n");
////            Rprintf("npop = %d", npop);
////            Rprintf("\n");
////            Rprintf("i = %d", i);
////            Rprintf("\n");
////            Rprintf("j = %d", j);
////            Rprintf("\n");
////            Rprintf("migmat = %f", pars.cp.mig_mat[i][j]);
//        }
//        Rprintf("\n");
//        }
    }
    
    if (en.ncol() == 3){ // for no en, pass matrix that does not have 3 columns
        double row_en[3];
        NumericVector numrow = en.nrow();
        for(n = 0; n < numrow.length(); n++){
            NumericVector rcpp_row_en = en(n,_);
            for(r = 0; r < 3; r++){ //hard coded as three because the en and ej should only have 3 numbers
                row_en[r] = rcpp_row_en[r];
            }
            pt = (struct devent *)malloc( sizeof( struct devent) ) ;
            pt->detype = 'n';
            
            pt->time = row_en[0];
            pt->nextde = NULL ;

            if( pars.cp.deventlist == NULL ){
                pars.cp.deventlist = pt ; 
		}	
            else if ( pt->time < pars.cp.deventlist->time ) {
                ptemp = pars.cp.deventlist ;
                pars.cp.deventlist = pt ;
                pt->nextde = ptemp ;
            }
            else {
                addtoelist( pt, pars.cp.deventlist ) ;}
            pt->popi =  row_en[1] -1 ;
            pt->paramv = row_en[2] ;
            }
        }
    
    if (ej.ncol() == 3){
        double row_ej[3]; //vector for each row
        NumericVector numrow_ej = ej.nrow();
        for(j = 0; j < numrow_ej.length(); j++){ //loop to go down each row
            Rcpp::NumericVector rcpp_row_ej = ej(j,_);
            for(r = 0; r < 3; r++){
                row_ej[r] = rcpp_row_ej[r]; //populatng the solo vector for the row
                }
            pt = (struct devent *)malloc( sizeof( struct devent) ) ;
            pt->detype = 'j';
            
            pt->time = row_ej[0]; //time is the first element in the row ??
            pt->nextde = NULL ;
            if( pars.cp.deventlist == NULL ){
                pars.cp.deventlist = pt ; // do i need to populate pt itself?
		}
            else if ( pt->time < pars.cp.deventlist->time ) {
                ptemp = pars.cp.deventlist ;
                pars.cp.deventlist = pt ;
                pt->nextde = ptemp ;
            }
            else{
                addtoelist( pt, pars.cp.deventlist ) ;}
            pt ->popi = row_ej[1] - 1; // do these need to go before or here?
            pt ->popj = row_ej[2] - 1;
            }
        }
}

