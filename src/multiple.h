#ifndef __multiple__
#define __multiple__

#include "btree_template.h"
#include "subalign.h"
#include "sequences.h"
#include "header_cpp.h"

class multiple {

    public:

	multiple(char *filename);

	char inputfileName[200];

	double distance_cutoff_similar;
	double distance_cutoff_div;

	sequences allseqs;

	btree<tnode> alltree;

	void set_distance_cutoff_similar(double dist_cutoff);

	void alignSimilar();

	void output_alignment();

	void alignDivergent();

	void store_similar(tnode *r);

	void map_allseqs_pos_to_tnode();

	int N_small;
	double **dist_matrix_preAligned;
	void get_distance_matrix_for_preAligned(int N_smallest);

	sequences aligned_seqs;

	void addSimilar();

};

subalign * merge_align_by_one_sequence(subalign *a, subalign *b, char *seq_name);

#endif
