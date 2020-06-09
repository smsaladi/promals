#ifndef __MM_
#define __MM_

// MM: global alignment by Myers and Miller algorithm

//extern int w(int x, int y);

class MM {

    public:

	MM();
	MM(int m, int n, int _g, int _gh);
        ~MM();

        void setM(int x);
	void setN(int x);
	void set_g(int x);
	void set_h(int x);

	void dp();
	void dp(int **smat);
	int *displ;
	int maxscore;
	int print_ptr;
	int endgappenalties;

	int **scoreMat;
 
    private:

	int M;  // length of the first dimension
	int N;  // length of the second dimension
	int g;  // gap opening penalty
	int gh;  // gap extension penalty

	int *HH;  // best cumulative 
	int *DD;  // best deletion 
	int *RR;  // 
	int *SS; // traceback direction
	int *gS;
	int last_print;

	int w(int x, int y) {return scoreMat[x][y];}

	int diff(int A,int B,int M,int N,int tb,int te); // simple Myers-Miller with affine gap penalties
	int diff1(int A,int B,int M,int N,int tb,int te); // also consider the situation of no end gap penalties
	int open_penalty1(int, int);
	int ext_penalty1(int, int);
	int gap_penalty1(int, int, int);
	int open_penalty2(int, int);
        int ext_penalty2(int, int);
        int gap_penalty2(int, int, int);
	void del(int k);
	void add(int k);
	void palign();

	int debug;

	void warning(char *s);

	
};

#endif
