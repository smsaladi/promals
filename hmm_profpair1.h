#ifndef _hmm_pairaln1__
#define _hmm_pairaln1__

#include "header_cpp.h"
#include "subalign.h"
#include "smallnumber.h"
#include "mathfunc.h"
#include "ScoreType.h"

// potentially important parameters: 
// 1: gap_threshold; 2: beta for profile calculation

// generalization of sum-of-pairs scoring
class  hmm_profpair1 {

    public:

	hmm_profpair1(subalign *x1, subalign *y1);
	~hmm_profpair1();

	// transition probabilities
	float tau, delta, epsilon, theta;
	
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
	void forward();
	ScoreType log_odds(int xi, int yj);
	ScoreType log_odds_x(int xi);
	ScoreType log_odds_y(int yj);

	// backward
	ScoreType **Bm, **Bx, **By, BE;
	float **probMat;
	void backward();

	// posterior
	int *posteriorAlignment();

	// putting together from beginning to end
	void getPosterior();
	
	void baumwelch();
	void getTransitions();
	float ep, ep1, d, d2,th, th2;
	float countEPx, countEP1x, countEPy, countEP1y, countDx, countDy, countD2m;
	float countEP, countEP1, countD, countD2;

	// 

};

extern void trainingHMM_profpair1(char *fastaFileList);
extern float newEP, newD;

#endif

