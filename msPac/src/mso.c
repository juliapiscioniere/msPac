#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "ms.h"
#include <R.h>

#define NL putchar('\n')
#define size_t unsigned

#define MIN(x, y) ( (x)<(y) ? (x) : (y) )

//#define ERROR(message) error(message)

#define SEGINC 80 


struct segl {
	int beg;
	struct node *ptree;
	int next;
};

double *posit ;
double *agevec ;
double alleleage ;
double segfac ;
int count, nseeds ; //ntbs


void gensam( char **list, double *pprobss, double *ptmrca, double *pttot, char *jtree)
{
	int nsegs, k, seg, ns, start, end, len, segsit ;
	struct segl *seglst, *segtre_mig(struct c_params *p, int *nsegs ) ;
	double nsinv,  tseg, tt, ttime(struct node *, int nsam), ttimemf(struct node *, int nsam, int mfreq) ;
	int segsitesin,nsites;
	double theta;
	int nsam, mfreq ;
	void prtree( struct node *ptree, int nsam, char *jtree);
	void make_gametes(int nsam, int mfreq,  struct node *ptree, double tt, int newsites, int ns, char **list );
 	void ndes_setup( struct node *, int nsam );

    
	nsites = pars.cp.nsites ;
	nsinv = 1./nsites;
	seglst = segtre_mig(&(pars.cp),  &nsegs ) ;
	nsam = pars.cp.nsam;
	segsitesin = pars.mp.segsitesin ;
    theta = pars.mp.theta ;
	mfreq = pars.mp.mfreq ;
    
    ns = 0 ;
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
      if( (pars.cp.r > 0.0 ) || (pars.cp.f > 0.0) ){
         end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
         start = seglst[seg].beg ;
         len = end - start + 1 ;
      }
      prtree( seglst[seg].ptree, nsam, jtree) ;
    }
    
    ns = 0 ;
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
        if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
        end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
        start = seglst[seg].beg ;
        len = end - start + 1 ;
        tseg = len*(theta/nsites) ;
        if( mfreq == 1) tt = ttime(seglst[seg].ptree, nsam);
        else tt = ttimemf(seglst[seg].ptree, nsam, mfreq );
        segsit = poisso( tseg*tt );
        if( (segsit + ns) >= maxsites ) {
            maxsites = segsit + ns + SITESINC ;
            posit = (double *)realloc(posit, maxsites*sizeof(double) ) ;
            agevec = (double *)realloc(agevec, maxsites*sizeof(double) ) ;
              biggerlist(nsam, list) ;
        }
        make_gametes(nsam,mfreq,seglst[seg].ptree,tt, segsit, ns, list );
        free(seglst[seg].ptree) ;
        locate(segsit,start*nsinv, len*nsinv,posit+ns);
        ns += segsit;
    }
}

	void 
ndes_setup(struct node *ptree, int nsam )
{
	int i ;

	for( i=0; i<nsam; i++) (ptree+i)->ndes = 1 ;
	for( i = nsam; i< 2*nsam -1; i++) (ptree+i)->ndes = 0 ;
	for( i= 0; i< 2*nsam -2 ; i++)  (ptree+((ptree+i)->abv))->ndes += (ptree+i)->ndes ;

}

	void
biggerlist(int nsam,  char **list)
{
	int i;
	for( i=0; i<nsam; i++){
	   list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
	   if( list[i] == NULL ) perror( "realloc error. bigger");
	   }
}
	   


/* allocates space for gametes (character strings) */
	char **
cmatrix(nsam,len)
	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned) nsam*sizeof( char* ) ) ) )
	   perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) len*sizeof( char ) )))
			perror("alloc error in cmatric. 2");
		}
	return( m );
}



	void
locate(int n,double beg, double len,double *ptr)
{
	int i;

	ordran(n,ptr);
	for(i=0; i<n; i++)
		ptr[i] = beg + ptr[i]*len ;

}

	void
addtoelist( struct devent *pt, struct devent *elist ) 
{
	struct devent *plast, *pevent, *ptemp  ;

	pevent = elist ;
	while(  (pevent != NULL ) && ( pevent->time <= pt->time ) )  {
		plast = pevent ;
		pevent = pevent->nextde ;
		}
	ptemp = plast->nextde ;
	plast->nextde = pt ;
	pt->nextde = ptemp ;
}

	void 
free_eventlist( struct devent *pt, int npop )
{
   struct devent *next ;
   int pop ;
   while( pt != NULL){
	  next = pt->nextde ;
	  if( pt->detype == 'a' ) {
	     for( pop = 0; pop < npop; pop++) free( (pt->mat)[pop] );
		 free( pt->mat );
	  }
	  free(pt);
    
	  pt = next ;
   }
}

#define STATE1 '1'
#define STATE2 '0'

	void
make_gametes(int nsam, int mfreq, struct node *ptree, double tt, int newsites, int ns, char **list )
{
	int  tip, j,  node ;
        int pickb(int nsam, struct node *ptree, double tt), 
            pickbmf(int nsam, int mfreq, struct node *ptree, double tt) ;

	for(  j=ns; j< ns+newsites ;  j++ ) {
		if( mfreq == 1 ) node = pickb(  nsam, ptree, tt);
		else node = pickbmf(  nsam, mfreq, ptree, tt);
        agevec[j] = alleleage ;
		for( tip=0; tip < nsam ; tip++) {
		   if( tdesn(ptree, tip, node) ) list[tip][j] = STATE1 ;
		   else list[tip][j] = STATE2 ;
        }
    }
}


/***  ttime.c : Returns the total time in the tree, *ptree, with nsam tips. **/

	double
ttime( ptree, nsam)
	struct node *ptree;
	int nsam;
{
	double t;
	int i;

	t = (ptree + 2*nsam-2) -> time ;
	for( i=nsam; i< 2*nsam-1 ; i++)
		t += (ptree + i)-> time ;
    /* adna */
        for( i=0; i<nsam; i++)
            t -= (ptree+i)->time ;
	return(t);
}


	double
ttimemf( ptree, nsam, mfreq)
	struct node *ptree;
	int nsam, mfreq;
{
	double t;
	int i;

	t = 0. ;
	for( i=0;  i< 2*nsam-2  ; i++)
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
		t += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
	return(t);
}


	void
prtree( ptree, nsam, jtree) //pointer or &?
	struct node *ptree;
	int nsam;
    char *jtree;
{
	int i, *descl, *descr ;
	void parens( struct node *ptree, int *descl, int *descr, int noden, char *jtree );

	descl = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	descr = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	for( i=0; i<2*nsam-1; i++) descl[i] = descr[i] = -1 ;
	for( i = 0; i< 2*nsam-2; i++){
	  if( descl[ (ptree+i)->abv ] == -1 ) descl[(ptree+i)->abv] = i ;
	  else descr[ (ptree+i)->abv] = i ;
	 }
	parens( ptree, descl, descr, 2*nsam-2, jtree);
	free( descl ) ;
	free( descr ) ;
}

	void
parens( struct node *ptree, int *descl, int *descr,  int noden, char *jtree)
{
	double time ;
   
	if( descl[noden] == -1 ) {
        char *main_parens = (char*)malloc(100 * sizeof(char));
        sprintf(main_parens, "%d:%5.3lf", noden+1, (ptree+ ((ptree+noden)->abv))->time  - (ptree+noden )->time ); /* adna */
        strcat(jtree, main_parens);
	}
 	else{
        strcat(jtree, "(");
        parens( ptree, descl,descr, descl[noden], jtree) ;
        strcat(jtree, ",");
        parens(ptree, descl, descr, descr[noden], jtree) ;
        if((ptree+noden)->abv == 0 ){
            strcat(jtree, ");");
        }
        else {
            time = (ptree + (ptree+noden)->abv )->time - (ptree+noden)->time ;
            char *time_char = (char*)malloc(100 * sizeof(char));
            sprintf(time_char, "):%5.3lf", time );
            strcat(jtree, time_char);
            free(time_char);
        }
    }
}

/***  pickb : returns a random branch from the tree. The probability of picking
              a particular branch is proportional to its duration. tt is total
	      time in tree.   ****/

	int
pickb(nsam, ptree, tt)
	int nsam;
	struct node *ptree;
	double tt;
{
	double x, y, ran1();
	int i;

	x = ran1()*tt;
	for( i=0, y=0; i < 2*nsam-2 ; i++) {
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
        if( y >= x ) {
            alleleage = (ptree+i)->time + ran1()* ( (ptree + (ptree+i)->abv )->time - (ptree+i)->time ) ;
            return( i ) ;
        }
        
    }
    i = 2*nsam - 3 ;
    alleleage = (ptree+i)->time + ran1()* ( (ptree + (ptree+i)->abv )->time - (ptree+i)->time ) ;
	return( 2*nsam - 3  );  /* changed 4 Feb 2010 */
}

	int
pickbmf(nsam, mfreq, ptree, tt )
	int nsam, mfreq;
	struct node *ptree;
	double tt;
{
	double x, y, ran1();
	int i, lastbranch = 0 ;

	x = ran1()*tt;
	for( i=0, y=0; i < 2*nsam-2 ; i++) {
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ){
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
		lastbranch = i ;    /* changed 4 Feb 2010 */
	  }
     if( y >= x ) {
            alleleage = (ptree+i)->time + ran1()* ( (ptree + (ptree+i)->abv )->time - (ptree+i)->time ) ;
            return( i ) ;
     }
	}
    i = lastbranch;
    alleleage = (ptree+i)->time + ran1()* ( (ptree + (ptree+i)->abv )->time - (ptree+i)->time ) ;
	return( lastbranch );   /*  changed 4 Feb 2010 */
}

/****  tdesn : returns 1 if tip is a descendant of node in *ptree, otherwise 0. **/

	int
tdesn(struct node *ptree, int tip, int node )
{
	int k;

	for( k= tip ; k < node ; k = (ptree+k)->abv ) ;
	if( k==node ) return(1);
	else return(0);
}


/* pick2()  */

	int
pick2(int n, int *i, int *j)
{
	double ran1();

	*i = n * ran1() ;
	while( ( *j = n * ran1() ) == *i )
		;
	return(0) ;
}

/**** ordran.c  ***/

	void
ordran(int n,double pbuf[])
{
	ranvec(n,pbuf);
	order(n,pbuf);
	return;
}


	void
mnmial(int n, int nclass, double p[], int rv[])
{
	double ran1();
	double x, s;
	int i, j;

	for(i=0; i<nclass; i++) rv[i]=0;
	for(i=0; i<n ; i++) {
	   x = ran1();
	   j=0;
	   s = p[0];
	   while( (x > s) && ( j<(nclass-1) ) )  s += p[++j];
	   rv[j]++;
	   }
	return;
}

        void
order(int n,double pbuf[])
{
        int gap, i, j;
        double temp;

        for( gap= n/2; gap>0; gap /= 2)
           for( i=gap; i<n; i++)
                for( j=i-gap; j>=0 && pbuf[j]>pbuf[j+gap]; j -=gap) {
                   temp = pbuf[j];
                   pbuf[j] = pbuf[j+gap];
                   pbuf[j+gap] = temp;
                   }
    return;
}


	void
ranvec(int n,double pbuf[])
{
	int i;
	double ran1();

	for(i=0; i<n; i++)
		pbuf[i] = ran1();

	return;
}



	int
poisso(double u)
{
	double  cump, ru, ran1(), p, gasdev(double, double) ;
	int i=1;

	if( u > 30. ){
	    i =  (int)(0.5 + gasdev(u,u)) ;
	    if( i < 0 ) return( 0 ) ;
	    else return( i ) ;
	  }
	 
	ru = ran1();
	p = exp(-u);
	if( ru < p) return(0);
	cump = p;
	
	while( ru > ( cump += (p *= u/i ) ) )
		i++;
	return(i);
}


/* a slight modification of crecipes version */

double gasdev(m,v)
	double m, v;
{
	static int iset=0;
	static float gset;
	float fac,r,v1,v2;
	double ran1();

	if  (iset == 0) {
		do {
			v1=2.0*ran1()-1.0;
			v2=2.0*ran1()-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset= v1*fac;
		iset=1;
		return( m + sqrt(v)*v2*fac);
	} else {
		iset=0;
		return( m + sqrt(v)*gset ) ;
	}
}
