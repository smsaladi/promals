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
	void printali(int blocksize, int header);
	void printali(char *filename, int blocksize);
	void printProfile();

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

	void set_gapt(double gapt);

	//subalign *sub2align(int *mark);
	//subalign *sub2align(int *mark, int *mark1);
	
	subalign *sub2align(int *mark);
	subalign *sub2align(int *mark, int *mark1);

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
	double average_sum_eff_let;
	
	double *gap_content;

	int alilen_mat;
	double n_eff;
	double gap_threshold;
	double gapt_threshold;
	double beta;

	void profile(); int done_profile;
        void neffsForEachCol_maskGapReg(int **ali, int n, int len, double effgapmax, double effgapRegionMin, double **n_effAa, double *sum_eff_let, int *maskgapRegion, int *apos_filtr, int *len_lowgaps, double *nef);
	double effective_number_nogaps(int **ali, int *marks, int n, int start, int end);
	void pseudoCounts(double **matrix, double n_eff, int len, double **pseudoCnt);

	void h_weight_all();
	double *hwt_all;

    public:
	double **distMat;


    public:
	// purge an alignment by
	// 1. removing highly similar sequences (id > high_id_thr)
	// 2. removing highly divergent sequences (id < low_id_thr)
	// 3. removing sequence fragment (gap_fraction > gap_fraction_thr)
	// 4. selecting only a subset if remaining num > max_num_kept
	subalign * purge_align(double low_id_thr, double high_id_thr, int max_num_kept, double gap_fraction_thr);

    // my new version of profile - effective position determined by Henikoff gap weight
    public:
	int *prof_pos;
	double **prof_effn;
	double **prof_freq;
	double *prof_sum_eff;
	double *prof_hwt_all;
	double *prof_gap_content;
	double prof_nef;
	double prof_average_sum_eff_let;
	int prof_len;
	double prof_raw_gap_threshold;
	double prof_gap_threshold;
	void prof_h_weight_all(double raw_gap_fraction_cutoff);
	void prof_positions(double prof_gap_threshold_here);
	void set_prof_raw_gap_threshold(double gapthr);
	void set_prof_gap_threshold(double gapthr);
	void prof();
	void prof_get_effn(int **ali, int n, int len, int prof_len, double **n_effAa, double *sum_eff_let, int *prof_pos, double *nef);
	int done_prof;

	
};


subalign * oneSeq2subalign(char *seq, char *name);

#endif
