#ifndef __math_function_
#define __math_function_

#include "mathnr.h"

void nrerror(char error_text[]);
double erfcc(double x);
double gammln(double xx);
double betacf(double a, double b, double x);
double betai(double a, double b, double x);
void pearsn(double x[], double y[], unsigned long n, double *r, double *prob, double *z);
double ran3(long *idum);
void moment(double data[], int n, double *ave, double *adev, double *sdev, double *var, double *skew, double *curt);
void moment2(double data[], int n, double *ave, double *sdev);

void balanc(double **a, int n);
void elmhes(double **a, int n);
void hqr(double **a, int n, double wr[], double wi[]);
void svdcmp(double **a, int m, int n, double w[], double **v);
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);

void sort(unsigned long n, double arr[]);
template<class kind> void sort(unsigned long n, kind arr[]);


#define NRANSI
//#include "nrutil.h"
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

template<class kind>
void sort(unsigned long n, kind arr[])
{
	unsigned long i,ir=n,j,k,l=1;
	int jstack=0,*istack;
	kind a,temp;

	//istack=ivector(1,NSTACK);
	istack = ivector(NSTACK);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				for (i=j-1;i>=1;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
				}
				arr[i+1]=a;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);
			}
			arr[l]=arr[j];
			arr[j]=a;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	//free_ivector(istack,1,NSTACK);
	delete [] istack;
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 21].,t45316. */

#endif
