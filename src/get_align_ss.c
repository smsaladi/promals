#include "all.h"
#include "sequences.h"
#include "btree_template.h"
#include "progressiveAlignHMM.h"
#include "multiple.h"

//static int debug = 1;
void configure_parameters();

int main(int argc, char **argv) {

	int i,j,k;

	getParameter(argc, argv, 1);
	//printParameters();
	configure_parameters();

	printParameters();
	printParameters1();

	//exit(0);
	
	//hmm_psipred_parameters *psi_param = new hmm_psipred_parameters(psipred_parameter_file, psipred_env_number, 1);
	//get_log_robinson_freq();
	//get_log_q_blosum62_ratio();
	//get_log_q_blosum62();

	//strcpy(blast_dir, ".");
	//strcpy(psipred_dir, ".");

	//if(probconsBLOSUM) useProbconsFrequencies();

	subalign *x = new subalign(argv[1]);

	x->select_representative();

	subalign *aux_align = get_blastpgp_alignment(x->aname[0], x->aname[0], x->aseq[0]);

	x->get_ss_prof1(blast_dir, x->aname[0], runpsipred1_command);

}

void configure_parameters() {

	int i, j, k;

	// psipred dir
	//strcpy(psipred_dir, "/home/jpei/astral/dp_1.69/psipred/psipred");

	// uniref_dir
	//set_uniref90("/home/jpei/uniref90/uniref90.fasta");
	//set_uniref90("/home/jpei/uniref90/uniref90_filt");

}
