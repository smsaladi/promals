#ifndef _hmm_local__
#define _hmm_local__

#include "header_cpp.h"
#include "subalign.h"
#include "smallnumber.h"
#include "mathfunc.h"
#include "ScoreType.h"

// potentially important parameters: 
// 1: gap_threshold; 2: beta for profile calculation

// generalization of sum-of-pairs scoring
class  hmm_local {

    public:

	hmm_local(subalign *x1, subalign *y1);
	~hmm_local();
	void deleteSubalign();

	// transition probabilities
	float tau, delta, epsilon, theta, eta;
	
	// two alignments
        subalign *x, *y;
	int lenx, leny;

	int gapopen_end;

	// viterbi
	float **Vm, **Vx, **Vy, VE; // maximum logrithm probability ending with m, x, or y or the last one E
	int **Tm, **Tx, **Ty, TE; // tracking 
	int *path, *path1, *path2; // path: 0, aligned position; 1, x to -; -1: - to y
	int v_alilen; // profile alignment length
	float log_odds_score(int xi, int yj);
	float log_odds_x_d(int xi);
	float log_odds_y_d(int yj);
	float max3(float a,float b, float c, int &t); // t stores the index of the max
	float max2(float a, float b, int &t); // t stores the index of the max
	void viterbi();
	void printHmmAlign(int blocksize); // print the alignment of the viterbi
	subalign *productViterbiAlign(); // generate an alignment of the two alignments

	// forward
	ScoreType **Fm, **Fx, **Fy, FE;
	ScoreType *Frx1, **Fry1, **Frx2, **Fry2, *Fn1, **Fn2, **Fn3, **Fn4, **Fe;
	void forward();
	ScoreType log_odds(int xi, int yj);
	ScoreType log_odds_x(int xi);
	ScoreType log_odds_y(int yj);
	void getTransitions();

	// backward
	ScoreType **Bm, **Bx, **By, BE;
	ScoreType **Brx1, **Bry1, **Brx2, *Bry2, **Bn1, **Bn2, **Bn3, *Bn4, **Be;
	float **probMat;
	void backward();

	// posterior

	// putting together from beginning to end
	void getPosterior();

	// Baum-Welch
	void  baumwelch();
	float estimatedE;
	float estimatedD;
	float estimatedEP;
	float estimatedT;
	float countE, countE1, countD, countD1, countEP, countEP1, countT, countT1;

	ScoreType L_eta, L1_eta, L_delta, L1_delta2_tau, L_tau, L_epsilon, L1_epsilon_tau;
	ScoreType e, e1, d, d2t1, t, ep, ept1;
};

extern float ne;
extern float nep;
extern float nt;
extern float nd;
void trainingLocalHMM(char *fastaFileList);

#endif

