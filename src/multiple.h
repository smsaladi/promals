#ifndef __multiple__
#define __multiple__

#include "btree_template.h"
#include "subalign.h"
#include "sequences.h"
#include "header_cpp.h"

extern FILE *logfp;

class multiple {

    public:

	multiple(char *filename);

	char inputfileName[200];

	double distance_cutoff_similar;
	double distance_cutoff_div;

	sequences allseqs;

	btree<tnode> alltree;

	void alignSimilar();

	void alignDivergent();
	void alignDivergent(float divergent_cutoff);

	void store_similar(tnode *r);

	void map_allseqs_pos_to_tnode();

	int N_small;
	double **dist_matrix_preAligned;
	void get_distance_matrix_for_preAligned(int N_smallest);

	sequences aligned_seqs;

        void set_distance_cutoff_similar(double dist_cutoff);
	void set_distance_cutoff_similar(double dist_cutoff, int Ngroup);
	void distance2leaf(tnode *r, double * array);
        void output_alignment();

	void align_profilehmm(int use_ss, char *ss_dir_name);

	void addSimilar();

	void alignDivergent_psipred(int x);
	void alignDivergent_psipred_sum_of_pairs(int x);

};

subalign * merge_align_by_one_sequence(subalign *a, subalign *b, char *seq_name);

extern char ssint2ss(int i);

#endif
