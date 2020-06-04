#ifndef __gap_refining__
#define __gap_refining__
#include "all.h"

extern double gap_matrix[21][21];

class group {

   public:
	int length;
	vector<int> seqnum;
	vector<string> gseq;

};

extern void gap_align_two_groups(group &a, group &b);

class gap_refine {

    public: 
	gap_refine();
	~gap_refine();

	// parameters
	int num_end_letters;
	int min_cb_size;
	int problem_gr_size;

	// 
	int nseqs;
	int alilen;
	string *seq;
	//string **cb;
	//string **gr;

	int *gr_begin, *gr_end, num_gr;

	// 
	void setup_seq(char **aseq, int n, int len);
	void setup_seq(subalign *x);
	void treat_end_letters();
	void define_cb_gr();
	void treat_continuous_gaps();

	void printseq(int blocksize);
	void printseq(int blocksize, int *mark);

	//
	double pb62[21][21];
	int done_get_pb62;
	void get_pb62();
	static void define_v(double v);

	// get groups
	vector<group> gr_group;

	void refine_group(int begin, int end); // do_refine mode 2
	void refine_group_by_seqnum(int begin, int end); // do_refine mode 1
	void refine_group_push_aside(int begin, int end); // do_refine mode 0
	void gap_align_two_groups(group &a, group &b);

	void do_refine(int mode);

	// 
	char **outalign();
	char **batch_refine(subalign *x, int mode);

};

#endif
