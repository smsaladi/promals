#ifndef _profilehmm__
#define _profilehmm__

#include "header_cpp.h"
#include "subalign.h"
#include "smallnumber.h"
#include "mathfunc.h"
#include "ScoreType.h"
#include "ss_prof.h"

// potentially important parameters: 
// 1: gap_threshold; 2: beta for profile calculation

// generalization of sum-of-pairs scoring
class  profilehmm {

    public:

	profilehmm(subalign *x1);

	~profilehmm();

	// build profile hmm
	void build_hmm(subalign *x1);

	subalign *x;  // the model is derived from x
	int nal, alilen_mat, alilen;
	int **alignment;
	subalign *y;  // the alignment or sequence to be evaluated
	int lenx, leny;
	void set_align(subalign *y1);

	int prof_len;
	int *prof_pos;
	double *prof_sum_eff;
	double *hwt_all;

	//consv *csv;

	// transition probabilities among states
	// states named according to Figure 1 in HMMER user guide
	float t_sn, t_sb, t_nn, t_nb;
	float *t_bm, t_bd1;
	float **t_m, **t_i, **t_d;
	float t_ec, t_et, t_ct, t_cc;

	// convert the transition probabilities to log scale
	void log_convert();

	// weight of the regularizer
	float reg_weight;
	
	int use_position_specific_regularizer;

	// number of alphabets of the second track
	int N2; // fixed to be 9

	// secondary structure profiles
	//ss_prof *ssx;
	//ss_prof *ssy;

	// set up parameters 
	void set_parameters(char *file_name, char *ss_dir_name, int use_ss);
	void print_transitions();

	// parameter estimation from an alignment

	// scores
	float **score_matrix; // for match states
	float *score_bg;      // for insertion states, N,C-terminal unaligned regions
	// weight for secondary structure scores
	float ss_weight;
	int use_picasso;  // picasso type scoring function
			  // divide the scores by the total number of effective counts
	void get_scores(int bg_type);

	// viterbi
	float **Vm, **Vi, **Vd;
	float *Vn, *Vb, *Ve, *Vc;
    
	// forward
	float **Fm, **Fi, **Fd;
	float *Fn, *Fb, *Fe, *Fc;
	float Ft;
	void forward();
	
	// backward
	float **Bm, **Bi, **Bd;
	float *Bn, *Bb, *Be, *Bc;
	float Bs;
	void backward();

	
	// putting together from beginning to end
	void get_match_prob();
	float **probMat;

};

#endif

