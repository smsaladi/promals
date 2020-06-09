#ifndef __consv_
#define __consv_

// this is the class that calculate various types of sequence frequencies and conservations
// fq[][]:  raw frequencies
// hfq[][]: henikoff weighted frequencies
// icfq[][]: independent count frequencies

#include "header_cpp.h"
#include "util.h"
#include "subalign.h"

class consv
{
public:
	int nal;
	int alilen;
	int **alignment;
	subalign *align;

	double gap_thr;

	double **fq;  
	double **hfq;
	double **icfq;
	double *h_oafq;
	double **effcount;
	double *gap_fraction;
	int *goodpos;
	double *hwt_all;
	int goodposnum;
	
	void freq();
	void h_freq();
	void ic_freq();

	consv();
	consv(const subalign &);
	consv(const consv &);
	~consv();
	const consv &operator=( const consv &);
	void h_weight_all();

private:
	void setgap_thr(double t);

	double *oaf_ip;
	double *h_oaf_ip; 
	double **u_oaf,**h_oaf;

	void h_weight(int **ali, int ip, double *hwt);
	void overall_freq(int **ali, int startp, int endp, int *mark, double *oaf);
	void overall_freq_wgt(int **ali,int startp,int endp,int *mark,double *wgt, double *oaf);
	double effective_number_nogaps(int **ali, int *marks, int n, int start, int end);
	double *entro_conv(double **f, int **ali, double *econv);
	double *pairs_conv(double **f,int **ali,int **matrix1,int indx,double *pconv);
	double *variance_conv(double **f, int **ali, double **oaf, double *vconv);
	
};

#endif

