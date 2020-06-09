#ifndef __subalign_
#define __subalign_

#include "util.h"
#include "amino.h"

class subalign {

    public:
        subalign();
	subalign(char *filename);
	subalign(string filename);
	// copy constructor
	subalign(const subalign &);
        ~subalign();

	// assign
	const subalign &operator=( const subalign &);

	void readali(char *filename);
	void printali(int blocksize);
	void printali(char *filename, int blocksize);

	int getNal();
        int getAlilen();
        int **getAlignment();
        char **getAseq();
        char **getAname();
	double *getEff_num_seq();

        void setNal(int n);
        void setAlilen(int len);
        void setAlignment(int **ali);
        void setAseq(char **seq);
        void setAname(char **name);
	void convertAseq2Alignment();

	//subalign *sub2align(int *mark);
	//subalign *sub2align(int *mark, int *mark1);
	
	subalign sub2align(int *mark);
	subalign sub2align(int *mark, int *mark1);

	static int getsubaligncount();
	void add_sequence(char *name, char *seq);


    public:
        char **aname, **aseq; /* aseq: [0,nal-1][0,alilen-1] */
        int **alignment; /* alignment: [1,nal][1,alilen] */
	int nal, alilen, *astart, *alen;
	int mnamelen;

	//double *eff_num_seq;
        //int *ngap;

	
    private:
	void *mymalloc(int size);
	char * strsave(char *str);
	char * strnsave(char *str, int l);
	char ** incbuf(int n, char **was);
	int * incibuf(int n, int *was);
	// convert the allocations by malloc to allocations by new
	void reassign();
	
	static int subaligncount;

    public:
        int *apos_filtr;
        double **n_effAa;
	double **pseudoCnt;
        double *sum_eff_let;
        int *maskgapRegion;
	double qmatrix[21][21];

	int alilen_mat;
	double n_eff;
	double gap_threshold;
	double gapt_threshold;
	double beta;

	void profile(); int done_profile;
        void neffsForEachCol_maskGapReg(int **ali, int n, int len, double effgapmax, double effgapRegionMin, double **n_effAa, double *sum_eff_let, int *maskgapRegion, int *apos_filtr, int *len_lowgaps, double *nef);
	double effective_number_nogaps(int **ali, int *marks, int n, int start, int end);
	void pseudoCounts(double **matrix, double n_eff, int len, double **pseudoCnt);

    public:
	double **distMat;
};


subalign * oneSeq2subalign(char *seq, char *name);

#endif
