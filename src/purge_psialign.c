#include "all.h"

// this program purges psiblast alignment
/*
      - initial step of puring the alignment
            high identity cutoff: 0.90
            low identity cutoff: 0.20
            fragment removal: gap to the first sequence >0.5
     - if the number of sequences is >1000, select the top 1000 sequences, ignore the rest.
     - then calculate the profile and determine the average effective number of sequences in the core block positions (eff_gap_content <0.5)
            if ave_eff_num <5: purge the initial alignment using the following criteria:
            high identity cutoff: 0.95
            low identity cutoff: 0.20
            fragment removal: gap to the first sequence >0.5

      output is a file that is named in the format of "argv[1].20_high_aveeffnum.aln", high is 90 or 95
*/

int main(int argc, char **argv) {

    int i,j;

    char align_name[300];

    strcpy(align_name, argv[1]);
    strcat(align_name, ".20_90");

    //cout << align_name << endl;

    subalign a(argv[1]);
    
    subalign *b;
    b = a.purge_align(0.2, 0.9, 1000, 0.5);
    b->printali(80, 1);

    delete_complete_gap_positions(b);
    b->profile();
	
    cout << b->average_sum_eff_let << endl;

    if(b->average_sum_eff_let<10) {
	strcpy(align_name, argv[1]);
	strcat(align_name, ".20_95");
	b = a.purge_align(0.2, 0.95, 1000, 0.5);
	delete_complete_gap_positions(b);
	//b->printali(80, 1);
	//b->printali(80, 0);
	b->profile();
    }

    sprintf(align_name, "%s_%1.1f.aln", align_name, b->average_sum_eff_let );
    //cout << align_name << endl;
    
    b->printali(align_name, 80);

    return 0;
}
