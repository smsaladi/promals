#include "all.h"
#include "sequences.h"
#include "btree_template.h"
#include "progressiveAlignHMM.h"
#include "multiple.h"

//static int debug = 1;

int main(int argc, char **argv) {

	int i,j,k;

	getParameter(argc, argv, 1);
	printParameters();

	get_log_robinson_freq();
	get_log_q_blosum62_ratio();
	get_log_q_blosum62();

	if(probconsBLOSUM) useProbconsFrequencies();

	multiple *my_align = new multiple(argv[1]);

	//cout << id_thr << endl;
	my_align->set_distance_cutoff_similar(1-id_thr);
	
	my_align->alignSimilar();
	my_align->alignDivergent();
}

