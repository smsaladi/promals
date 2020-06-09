#include "all.h"
#include "sequences.h"
#include "btree_template.h"
#include "progressiveAlignHMM.h"
#include "multiple.h"

//static int debug = 1;
void configure_parameters();

int main(int argc, char **argv) {

	//int i,j,k;

	//sequences tmp1;
	//tmp1.readFasta(argv[1], 1);
        time(&timestart);
        timecount=0;

	getParameter(argc, argv, 1);
	//printParameters();
	configure_parameters();

	printParameters();
	printParameters1();

	//exit(0);
	
	//hmm_psipred_parameters *psi_param = new hmm_psipred_parameters(psipred_parameter_file, psipred_env_number, 1);

	//psi_param->print_parameters();

	get_log_robinson_freq();
	get_log_q_blosum62_ratio();
	get_log_q_blosum62();

	//strcpy(blast_dir, ".");
	//strcpy(psipred_dir, ".");

	if(probconsBLOSUM) useProbconsFrequencies();

	multiple *my_align = new multiple(argv[1]);
        //cout << "Here" << endl;

	if(my_align->allseqs.nseqs <= my_align->max_cluster_elem) {
                my_align->build_tree();
                my_align->set_distance_cutoff_similar(1-id_thr, max_group_number);
                my_align->alignSimilar();
        }
        else {
                my_align->build_cluster();
        }

        // remove filter files
        char command[500];
        sprintf(command, "rm -f %s_fiLtEr*", argv[1]);
        system(command);

	my_align->alignDivergent_psipred(2);
}

void configure_parameters() {

	//int i, j, k;

	// psipred dir
	//strcpy(psipred_dir, "/home/jpei/astral/dp_1.69/psipred/psipred");

	// uniref_dir
	//set_uniref90("/home/jpei/uniref90/uniref90.fasta");
	//set_uniref90("/home/jpei/uniref90/uniref90_filt");
	//set_uniref90("/home/jpei/uniref90_02282008/uniref90_filt");

}
