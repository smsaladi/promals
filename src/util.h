#ifndef __util_pei_
#define __util_pei_

double *dvector(int d1);
int *ivector(int d1);
char *cvector(int d1);
double **dmatrix(int d1, int d2);
int **imatrix(int d1, int d2);
char **cmatrix(int d1, int d2);
void free_dmatrix(double **mat, int m, int n);
void free_imatrix(int **mat, int m, int n);
void free_cmatrix(char **mat, int m, int n);
void printinfo(char *info, int doprint);


// generic vector
template<class kind>
kind *gvector(int d1) {
	kind *v;
	v = new kind [d1+1];
	return v;
}


// generic matrix
template<class kind>
kind **gmatrix(int d1, int d2) {
	int i;
	kind **m;
	m = new kind * [d1+1];
	for(i=0;i<=d1;i++) {
	   m[i] = new kind [d2+1];
	}
	return m;
}

// generic matrix free memory
template<class kind>
void free_gmatrix(kind **gmat, int d1, int d2) {

	int i;
	if(!gmat) return;
	for(i=0;i<=d1;i++) {
		delete [] gmat[i];
	}
	delete [] gmat;
	gmat = 0;
}

// generic matrix free memory
template<class kind>
void free_gvector(kind *gvec) {

	if(!gvec) return;
	delete [] gvec;
	gvec = 0;
}

/*
//static float sqrarg;
//#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

//static float maxarg1,maxarg2;
//#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

//static float minarg1,minarg2;
//#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ? (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

*/

#endif

