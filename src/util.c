#include "header_cpp.h"
#include "util.h"
//#include "matrixalgebra.c"

double *dvector(int d1) {
	int i;
	double *v;
	v = new double [d1+1];
	for(i=0;i<=d1;i++) v[i] = 0;
	return v;
}

int *ivector(int d1) {
	int i;
	int *v;
	v = new int [d1+1];
	for(i=0;i<=d1;i++) v[i] = 0;
	return v;
}

char *cvector(int d1) {
	char *v;
	v = new char [d1+1];
	return v;
}

double **dmatrix(int d1, int d2) {
	int i,j;
	double **m;
	// cout << "==========-----"<<endl;
	// cout << d1 << "\t" << d2 << endl;
	m = new double * [d1+1];
	// cout << "==========-----"<<endl;
	for(i=0;i<=d1;i++) {
	   m[i] = new double [d2+1];
	}
	// cout << "==========-----"<<endl;
	for(i=0;i<=d1;i++) {
	   for(j=0;j<=d2;j++) {
		m[i][j] = 0;
	   }
	}
	return m;
}

int **imatrix(int d1, int d2) {
	int i,j;
	int **m;
	m = new int * [d1+1];
        for(i=0;i<=d1;i++) {
           m[i] = new int [d2+1];
        }
	for(i=0;i<=d1;i++) {
	    for(j=0;j<=d2;j++) {
		m[i][j] = 0;
	    }
	}
        return m;
} 

char **cmatrix(int d1, int d2) {
        int i;
        char **m;
        m = new char * [d1+1];
        for(i=0;i<=d1;i++) {
           m[i] = new char [d2+1];
        }
        return m;
}

void free_dmatrix(double **mat, int m, int n) {

	int i;
	if(mat==0) return;
	// cout << m << "          " << n << endl;
	for(i=0;i<=m;i++) {
		delete [] mat[i];
		// fprintf(stdout, "=---------------\n"); fflush(stdout);
	}
	delete [] mat;
	mat = 0; // The setting to NULL does not work sometimes
}


void free_imatrix(int **mat, int m, int n) {

	int i;
	if(mat==0) return;
	for(i=0;i<=m;i++) {
		delete [] mat[i];
	}
	delete [] mat;
	mat = 0;
}


void free_cmatrix(char **mat, int m, int n) {

	int i;
	
	if(mat==0) return;

	for(i=0;i<=m;i++) {
		delete [] mat[i];
	}
	delete [] mat;
	mat = 0;
}

void printinfo(char *info, int doprint) {
        if(doprint) {
                cout << info << endl;
        }
}
