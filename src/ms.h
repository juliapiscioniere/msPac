struct devent {
	double time;
	int popi;
	int popj;
	double paramv;
    int num ;  /* adna */
	double **mat ;
	char detype ;
	struct devent *nextde;
	} ;
struct c_params {
	int npop;
	int nsam; /* total sample size including ancient samples */
    int nsamin;  /*  sample size at present time,  does not include ancient samples */
	int *config;
	double **mig_mat;
	double r;
	int nsites;
	double f;
	double track_len;
	double *size;
	double *alphag;
	struct devent *deventlist ;
	} ;
struct m_params {
    double theta;
	int segsitesin;
	//int treeflag;
	//int timeflag;
	int mfreq;
    //int ageflag ;
	 } ;
struct params { 
	struct c_params cp;
	struct m_params mp;
	//int commandlineseedflag ;
	int output_precision;
	};

struct node{
	int abv;
	int ndes;
	float time;
};

#include <R.h>
#define SITESINC 10

#ifdef __cplusplus
extern "C" struct params pars;
extern "C" unsigned maxsites;
#else
extern struct params pars;
extern unsigned maxsites;
#endif

/*KRT -- prototypes added*/
void ordran(int n, double pbuf[]);
void ranvec(int n, double pbuf[]);
void order(int n, double pbuf[]);

void biggerlist(int nsam,  char **list);
int poisso(double u);
void locate(int n,double beg, double len,double *ptr);
void mnmial(int n, int nclass, double p[], int rv[]);
int tdesn(struct node *ptree, int tip, int node );
int pick2(int n, int *i, int *j);
int xover(int nsam,int ic, int is);
int links(int c);

void gensam( char **list, double *probss, double *ptmrca, double *pttot, char *tree) ;
/* void seedit( const char * ) ; */
void addtoelist( struct devent *pt, struct devent *elist );
char ** cmatrix(int nsam, int len);
int commandlineseed(int *seeds);
