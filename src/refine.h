#ifndef __refine__
#define __refine__

#include "header_cpp.h"
#include <algorithm>
#include "subalign.h"
#include <map>
#include "consv1.h"

struct gapinfo {

	char aa;  // amino acid
	int lg;   // number of gaps to the left
	int rg;   // number of gaps to the right
	int single;  // is a single letter: -P-
	int doublet; // is the first of a double letters: -PX-
};

class alignrefine {
	
    public:
	subalign *a;
	char **aseq;
	int nal;
	int alilen;
	gapinfo **gi;  // gap information
	double *sw; // [1..nal] sequence weight for sequences

	void get_sw(); // by henikoff weight scheme in consv
		       // sum of weights for each sequence should be 1

	double *aaw; // [1..alilen] Sum-of-weight for all non-gapped letters for positions
	double gappy_threshold;
	int CORE_POS_NUMBER;
	void set_gappy_threshold(double g);
	void set_CORE_POS_NUMBER(int n);
	double *sdw; // [1..alilen] Sum-of-weight for letters that are single letters or doublets for positions

	double totalweight;

	alignrefine(subalign *x, int NC_added_gaps);
	~alignrefine();

	void assign_gi();
	void print_gi_html();
	void assign_weight(double *weight);
	void calculate_weights();

	// deal with single and double letters
	double single_fraction;
	void operation();
	void operation_together(int p);
	void operation_separate(int p);
	void shift_left(int n, int p); // shift a signle letter 
	void shift_right(int n, int p); // shift a single letter
	void treat_single_and_doublet(); // batch operations
	void gi2aseq();

	// deal with gappy regions
	// get rid of positions with all gaps; and reassign nal, alilen, ..
	void refresh_alignment();
	// for n_terminal gappy regions; shift letters to the right
	void deal_with_N();
	void deal_with_C();
	void shift_gappy_right(int gb, int ge, int n);
	void shift_gappy_left(int gb, int ge, int n);

	vector <double *> w_aa;
	vector <int> isgappy;
	vector <double> vaaw;

	void deal_with_gappy(int filter_length);
	void deal_with_gappy(int gb, int ge, int filter_length);
	void delete_all_gap(int b, int e);
	void delete_NC_terminal_gaps();
	void deal_with_gappy_individual(int gb, int ge, int n, int lmr);
	string add_middle_gaps(string mstr, int length);
	string add_right_gaps(string mstr, int length);
	string add_left_gaps(string mstr, int length);
	double compare_score(string a1, string a2, vector<int> b, int n);

};

#endif
