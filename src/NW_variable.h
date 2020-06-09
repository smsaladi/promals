#ifndef __NW_variable_
#define __NW_variable_

// NW_variable: Needleman-Wunch

extern int w(int x, int y);

class NW_variable {

    public:

	NW_variable();
	//NW_variable(int m, int n, int _g, int _h);
	NW_variable(int m, int n, double *_go1, double *_ge1, double *_go2, double *_ge2);
        ~NW_variable();

        void setM(int x);
	void setN(int x);
	void set_g(double x);
	void set_h(double x);
	void set_use_end_penalty(int x);
	void set_go1(double *x);
	void set_go2(double *x);
	void set_ge1(double *x);
	void set_ge2(double *x);
	void set_w(double **sm);

	void dp();
	void traceback();
	int *tr; // traceback path record [pt-1..0] from the end to the beginning
	int *r_tr; // reverse traceback path record [1..pt] from the beginning to the end
		   // -1: deletion in the second dimension (x -)
		   // 0: substitution (x y)
		   // 1: insertion in the second dimension (- y)
	double maxscore;
	int pt; // alignment length; number of positions in the traceback path
 
    private:

	int M;  // length of the first dimension
	int N;  // length of the second dimension
	double g;  // gap opening penalty
	double h;  // gap extension penalty

	double *go1;
	double *go2;
	double *ge1;
	double *ge2;
	
	double **w;

	double **D;  // best deletion 
	double **C;  // best cumulative
	double **I;  // best insertion
	int **T; // traceback direction
	int **T_I;  // trace back insert
	int **T_D; // trace back delete

	int use_end_penalty;

	int debug;

	int tracedirection; // trace direction -1 0 1
	int tracedirection_I; // trace direction insert
	int tracedirection_D; // trace direction delete
	double MAX_I(double x, double y); 
	double MAX_D(double x, double y);
	double MAX3(double x, double y, double z);
	void warning(char *s);

	
};

#endif
