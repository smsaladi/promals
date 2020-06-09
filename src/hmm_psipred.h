#ifndef _hmm_psipred__
#define _hmm_psipred__

#include "subalign.h"
#include "header_cpp.h"
#include "smallnumber.h"
#include "mathfunc.h"
#include "ScoreType.h"

// potentially important parameters: 
// 1: gap_threshold; 2: beta for profile calculation

// hmm parameters
class hmm_psipred_parameters {

    public:

	int env_number;
	float ***aa_pair;
	float **aa_loop;
	float **ss_pair;
	float **ss_loop;
	float ***aa_pair1;
	float **aa_bg1; // background frequencies estimated from aa_pair
	float **aa_loop1; // every name ending with a 1 is raw frequences; 
			  // otherwise logarithms are taken
	float **ss_pair1;
	float **ss_loop1;
	float *tran_begin;
	float **tran_match;
	float **tran_x;
	float **tran_y;

	hmm_psipred_parameters(int num) { env_number = num; }
	hmm_psipred_parameters(char *filename, int num, int log_mark_emission) { 
		read_parameters(filename, num, log_mark_emission);		
	}
	void read_parameters(char *filename, int num, int log_mark_emission);
	void print_parameters();
	float mylog(float a);
	float mylog0(float a);
};
		

// generalization of sum-of-pairs scoring
class  hmm_psipred {

    public:
	hmm_psipred(subalign *x1, subalign *y1) {
		x=x1; y=y1; lenx = x1->alilen; leny = y1->alilen; 
		Fm = NULL; Fx = NULL; Fy = NULL;
		Bm = NULL; Bx = NULL; By = NULL;
		score_matrix = NULL;
		score_bg_x = NULL;
		score_bg_y = NULL;
		probMat = NULL;
	}
	~hmm_psipred();

	// parameters
	hmm_psipred_parameters *p;

	// two alignments
        subalign *x, *y;
	int lenx, leny;
	
	void set_parameters(hmm_psipred_parameters *pointer1);

	int env_number;
	float ***aa_pair;
	float **aa_loop;
	float **ss_pair;
	float ***aa_pair1;
	float **aa_loop1;
	float **ss_pair1;
	float **ss_loop;
	float **ss_loop1;
	float *tran_begin;
	float **tran_match;
	float **tran_x;
	float **tran_y;
	int *x_alphabet;
	int *y_alphabet;

	// viterbi
	//float ***V;
	//float **Vm, **Vx, **Vy, VE; // maximum logrithm probability ending with m, x, or y or the last one E
	//int ***T;
	//int TE;
	//float VE;
	//int **Tm, **Tx, **Ty, TE; // tracking 
	//int *path, *path1, *path2; // path: 0, aligned position; 1, x to -; -1: - to y
	//int v_alilen; // profile alignment length
	//float log_odds_score(int xi, int yj);
	//float log_odds_x_d(int xi);
	//float log_odds_y_d(int yj);
	//float max3(float a,float b, float c, int &t); // t stores the index of the max
	//float max2(float a, float b, int &t); // t stores the index of the max
	//float max(float *list, int &t);
	//void viterbi(int lor);
	//void printHmmAlign(int blocksize); // print the alignment of the viterbi
	//subalign *productViterbiAlign(); // generate an alignment of the two alignments

	float **score_matrix;
	float *score_bg_x;
	float *score_bg_y;
	void get_scores(double ss_weight);
	void get_scores(double ss_weight, double score_weight);
	void get_scores_sum_of_pairs(double ss_weight);
	void get_scores_sum_of_pairs(double ss_weight, double score_weight);
	

	// forward
	ScoreType **Fm, **Fx, **Fy, FE;
	//ScoreType ***F, FE;
	void forward();
	void forward1();
	void forward_no_end_penalty();
	//ScoreType log_odds(int xi, int yj);
	//ScoreType log_odds_x(int xi);
	//ScoreType log_odds_y(int yj);

	// backward
	ScoreType **Bm, **Bx, **By, BE;
	//ScoreType ***B, BE;
	float **probMat;
	void backward();
	void backward1();
	void backward_no_end_penalty();

	//
	void forward2();
	void backward2();

	// posterior

	// putting together from beginning to end
	void getPosterior();

	// posterior alignment
	int *posteriorAlignment();

};

#endif

