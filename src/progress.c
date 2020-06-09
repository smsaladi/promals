#include "all.h"
#include "sequences.h"
#include "btree_template.h"
#include "progressiveAlignHMM.h"
#include "multiple.h"

//static int debug = 1;
void configure_parameters();

int main(int argc, char **argv) {

	int i,j,k;

	//sequences tmp1;
	//tmp1.readFasta(argv[1], 1);

	getParameter(argc, argv, 1);
	//printParameters();
	configure_parameters();

	printParameters();
	printParameters1();

	//exit(0);
	
	hmm_psipred_parameters *psi_param = new hmm_psipred_parameters(psipred_parameter_file, psipred_env_number, 1);

	//psi_param->print_parameters();

	get_log_robinson_freq();
	get_log_q_blosum62_ratio();
	get_log_q_blosum62();

	//strcpy(blast_dir, ".");
	//strcpy(psipred_dir, ".");

	if(probconsBLOSUM) useProbconsFrequencies();

	multiple *my_align = new multiple(argv[1]);

	//my_align->set_distance_cutoff_similar(1-id_thr);
	my_align->set_distance_cutoff_similar(1-id_thr, max_group_number);
	//my_align->alltree.writeTree("tmp.tre");
	//exit(0);
	
	my_align->alignSimilar();
	//my_align->alignDivergent();
	my_align->alignDivergent_psipred(2);
}

void configure_parameters() {

	int i, j, k;

	// psipred dir
	//strcpy(psipred_dir, "/home/jpei/astral/dp_1.69/psipred/psipred");

	// uniref_dir
	//set_uniref90("/home/jpei/uniref90/uniref90.fasta");
	//set_uniref90("/home/jpei/uniref90/uniref90_filt");

}
