#ifndef _hmm_multim__
#define _hmm_multim__

#include "header_cpp.h"
#include "subalign.h"
#include "smallnumber.h"
#include "mathfunc.h"
#include "ScoreType.h"

// potentially important parameters: 
// 1: gap_threshold; 2: beta for profile calculation

// hmm parameters
class hmm_parameters {

    public:
	int num_solv;
	int num_ss;
	int num_states;
	int num_match_states;
	int num_indel_states;
	int unaligned;
	float **transition_matrix;
	float ***match_emission_matrix;
	float *indel_emission_vector1;
	float *indel_emission_vector2;
	float *background_frequencies;

	float **t_log;
	float ***m_lor; // match_emission matrix log-odss ratio
	float *indel_lor1;
	float *indel_lor2;
	float ***m_log; // match_emission matrix log-odss ratio
	float *indel_log1;
	float *indel_log2;

	hmm_parameters(int solv, int ss, int unaligned_match) {
		num_solv = solv; num_ss = ss;
		num_match_states = solv * ss + unaligned_match;
		num_states = num_match_states + 4;
		num_indel_states = 2;
		unaligned = unaligned_match;
	}
	void read_parameters(char *filename);
	void print_parameters();
	
};
		


// generalization of sum-of-pairs scoring
class  hmm_multim {

    public:

	hmm_multim(subalign *x1, subalign *y1) {
		x=x1; y=y1; lenx = x1->alilen;leny = y1->alilen; 
		F = NULL; B = NULL; probMat = NULL;
	}
	~hmm_multim();

	// parameters
	hmm_parameters *p;

	void set_parameters(hmm_parameters *pointer1) {p=pointer1;
		num_match_states = p->num_match_states;
		p1 = num_match_states+1;
		p2 = num_match_states+2;
		p_end = num_match_states+3;
		unaligned = p->unaligned;
	}

	// transition probabilities
	//float tau, delta, epsilon;
	
	// two alignments
        subalign *x, *y;
	int lenx, leny;

	int num_match_states, p1, p2, p_end; // p1: insertion state 1; p2: insertion state 2; p_end: ending state
	int unaligned;

	int gapopen_end;

	// viterbi
	float ***V;
	//float **Vm, **Vx, **Vy, VE; // maximum logrithm probability ending with m, x, or y or the last one E
	int ***T;
	int TE;
	float VE;
	//int **Tm, **Tx, **Ty, TE; // tracking 
	int *path, *path1, *path2; // path: 0, aligned position; 1, x to -; -1: - to y
	int v_alilen; // profile alignment length
	//float log_odds_score(int xi, int yj);
	//float log_odds_x_d(int xi);
	//float log_odds_y_d(int yj);
	//float max3(float a,float b, float c, int &t); // t stores the index of the max
	//float max2(float a, float b, int &t); // t stores the index of the max
	float max(float *list, int &t);
	void viterbi(int lor);
	void printHmmAlign(int blocksize); // print the alignment of the viterbi
	subalign *productViterbiAlign(); // generate an alignment of the two alignments

	// forward
	//ScoreType **Fm, **Fx, **Fy, FE;
	ScoreType ***F, FE;
	void forward();
	//ScoreType log_odds(int xi, int yj);
	//ScoreType log_odds_x(int xi);
	//ScoreType log_odds_y(int yj);

	// backward
	//ScoreType **Bm, **Bx, **By, BE;
	ScoreType ***B, BE;
	float **probMat;
	void backward();

	// posterior

	// putting together from beginning to end
	void getPosterior();

	// posterior alignment
	int *posteriorAlignment();

};

#endif

