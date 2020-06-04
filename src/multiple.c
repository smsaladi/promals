#include "subalign.h"
#include "sequences.h"
#include "btree_template.h"
#include "progressiveAlignHMM.h"
#include "multiple.h"

static int debug_here = 11;

FILE *logfp;
static int arrayindex = 0;

multiple::multiple(char *filename) {

	strcpy(inputfileName, filename);
	distance_cutoff_similar = 0.4;
	distance_cutoff_div = 0.7;
	N_small = 20;
	
	if(relax_count>0) N_small = relax_count;

	allseqs.readFasta(filename, 1);

	sequences allseqs1(allseqs);

	/* Old way of getting distances
	// k-mer dist matrix calculation
	allseqs1.toDayhoff6();
	cout << "finished toDayhoff6" << endl;
	allseqs1.generateD6t(6);
	cout << "finished generateD6t" << endl;
	allseqs1.diffCountD2t(allseqs1.d6t[1], allseqs1.d6t[2]);
	cout << "finished diffCountD2t" << endl;
	allseqs1.d6t2DistMat(6);
	cout << "finished d6t2DistMat" << endl;
	*/

	// a much simpler and much faster way of getting distances
	allseqs1.get_kmer_array(6);
	allseqs1.get_kmer_distance(6);

	allseqs.distMat = allseqs1.distMat; // adapt the distance matrix

	// build k-mer dist based tree
	alltree.UPGMA(allseqs1.distMat, allseqs.seq, allseqs1.name, allseqs1.nseqs);
	//cout << "finished UPGMA" << endl;
	//exit(0);
}

void multiple::set_distance_cutoff_similar(double dist_cutoff) {

      distance_cutoff_similar = dist_cutoff;

}

// for a large number of sequences, find the cutoff that results in a fixed number of pre-aligned groups
void multiple::set_distance_cutoff_similar(double dist_cutoff, int Ngroup) {

	char logfile[500];
	strcpy(logfile, inputfileName);
	strcat(logfile, ".promals.logfile");
	logfp = fopen(logfile, "w");

	int i;
	if(allseqs.nseqs<=Ngroup) {
		distance_cutoff_similar = dist_cutoff;
		return;
	}

	// find the distance to the leaf node for each node
	double dist2leaf[allseqs.nseqs];

	arrayindex = 1;
	distance2leaf(alltree.root, dist2leaf);
	//for(i=1;i<=allseqs.nseqs-1;i++) { cout << dist2leaf[i] << endl; }
	sort(allseqs.nseqs-1, dist2leaf);
	//for(i=1;i<=allseqs.nseqs-1;i++) { cout << dist2leaf[i] << endl; }
	double tmp_dist_cutoff = (dist2leaf[allseqs.nseqs-1-(Ngroup-1)]+dist2leaf[allseqs.nseqs-1-(Ngroup-2)])/2 * 2;
	if(tmp_dist_cutoff > dist_cutoff) {
		distance_cutoff_similar = tmp_dist_cutoff;
		fprintf(logfp, "The original identity threshold %3.2f results in the number of pre-aligned groups > %d.\n", 1-dist_cutoff, Ngroup);
		fprintf(logfp, "The new identity threshold is adjusted to %3.2f to set the number of pre-aligned groups to %d.\n\n", 1-tmp_dist_cutoff, Ngroup);
		fprintf(stdout, "  identity threshold adjusted to: %3.2f (to save time)\n\n", 1-distance_cutoff_similar);
	}
	else distance_cutoff_similar = dist_cutoff;
	//cout << dist_cutoff << "\t"<< tmp_dist_cutoff << endl;
	//cout << "\t: " << distance_cutoff_similar  << endl;
	//exit(0);
}

void multiple::distance2leaf(tnode *r, double * array) {
	
	double dist;
	if(r->childL == NULL) {
		return;
	}
	dist = 0;
	tnode *tmp = r;
	while(tmp->childL!=NULL) {
		dist += tmp->childL->branchlen;
		tmp = tmp->childL;
	}
	array[arrayindex] = dist;
	//cout << arrayindex << " " << dist << endl;
	arrayindex++;
	distance2leaf(r->childL, array);
	distance2leaf(r->childR, array);
}

void multiple::alignSimilar() {

	int i, j;

	// progressively align similar sequences using general substitution matrix
	fprintf(logfp, "Start aligning similar sequences ...");
    	fflush(logfp);
	alltree.progressiveAlignHMM_FastStage(alltree.root, distance_cutoff_similar/2);
	fprintf(logfp, " Done.\n");
    	fflush(logfp);
        // if the root is aligned, stop
        if(alltree.root->aligned)  {
                //output_alignment();
		//cout << "NUMBER OF SEQUENCES: " << allseqs.nseqs << " NUMBER OF GROUPS: " << 1  << endl;
		//fprintf(logfp, "Number of input sequences: %d\n",  allseqs.nseqs);
		//fprintf(logfp, "Number of pre-aligned groups: 1\n");
		//fprintf(logfp, "PROMALS is now finished\n");
		//fclose(logfp);
                //exit(0);
        }


	// store pre-aligned groups in the stopped nodes and select one representative from each group
	store_similar(alltree.root);
	//cout << "here........" << endl;

	//exit(0);

	alltree.obtainPreAligned(alltree.root);
	//cout << "here........" << endl;

	fprintf(logfp, "      Number of input sequences: %d\n",  allseqs.nseqs);
	fprintf(logfp, "      Number of pre-aligned groups: %d\n", alltree.preAligned.size());
	fprintf(logfp, "\n");
    	fflush(logfp);

	map_allseqs_pos_to_tnode();
	//cout << "here........" << endl;

	get_distance_matrix_for_preAligned(N_small);
	//cout << "here........" << endl;

	// output distance matrix for representatives
	if(debug>1) {
		for(i=0;i<alltree.preAligned.size();i++) {
			for(j=0;j<alltree.preAligned.size();j++) {
				cout << dist_matrix_preAligned[i][j] << " ";
			}
			cout << endl;
		}
	}
	//cout << "here........" << endl;

	if(debug_here>1) {
	    cout << "NUMBER OF SEQUENCES: " << allseqs.nseqs << " NUMBER OF GROUPS: " << alltree.preAligned.size() << endl;
	    //exit(0);
	}
}

void multiple::alignDivergent() {

	int i, j;
		
	// right now, just the option of multim - for probablistic consistency 
        hmm_parameters params(solv,ss,unaligned);
        params.read_parameters(parameter_file);
	if(debug_here>11) { cout << "Before consistency" << endl; }
        alltree.profileConsistency_multim(&params, dist_matrix_preAligned, 1);

	if(debug_here>11) { cout << "After consistency" << endl; }

	alltree.computeConsistencyAlignment(alltree.root);

	output_alignment();
}

void multiple::alignDivergent_psipred(int use_homologs) {

	int i, j;

	if(alltree.preAligned.size()==1) {
		fprintf(logfp, "Start getting structural information ...\n");
	}
	else fprintf(logfp, "Start aligning divergent groups ...\n");
    	fflush(logfp);
	// right now, just the option of multim - for probablistic consistency 
        hmm_psipred_parameters params(psipred_env_number);
        params.read_parameters(psipred_parameter_file, psipred_env_number, 1);

	//cout << "Here: " << endl;

	//cout << psipred_dir << endl;

	// select representatives
	for(i=0;i<alltree.preAligned.size();i++) {
		alltree.preAligned[i]->aln->select_representative();
	}

	// use_homologs==1: use profile of pre-aligned group
	if(use_homologs==1) {
	     for(i=0;i<alltree.preAligned.size();i++)
		alltree.preAligned[i]->aux_align = alltree.preAligned[i]->similarSet->purge_align_one_seq_name(alltree.preAligned[i]->aln->aname[0]);
	}
	// use_homologs==2: use database homologs
	else if(use_homologs==2) {
	    fprintf(logfp, "      - Running PSI-BLAST and PSIPRED ...\n");
	    fprintf(logfp, "        ");
    	    fflush(logfp);
	    for(i=0;i<alltree.preAligned.size();i++) {
		alltree.preAligned[i]->aux_align = get_blastpgp_alignment(alltree.preAligned[i]->aln->aname[0], alltree.preAligned[i]->aln->aname[0], alltree.preAligned[i]->aln->aseq[0]);
		//cout << alltree.preAligned[i]->aln->aname[0] << " " <<  alltree.preAligned[i]->aln->aname[0] << endl; exit(0);
		alltree.preAligned[i]->aln->get_ss_prof1(blast_dir, alltree.preAligned[i]->aln->aname[0], runpsipred1_command);
		if(!alltree.preAligned[i]->aux_align) {
			cout << "reading blastpgp alignment failed " << alltree.preAligned[i]->aln->aname[0] << endl;
			alltree.preAligned[i]->aux_align = alltree.preAligned[i]->aln;
		}
		//fprintf(logfp, "             Repres. sequence %d\r", i+1);
		if(clean_blast_after) clean_blast_psipred(alltree.preAligned[i]->aln->aname[0]);
		fprintf(logfp, "*");
    	    	fflush(logfp);
	    }
	    fprintf(logfp, "\n");
	}
	// use_homologs==0: use single sequence as profile
	else {
	     for(i=0;i<alltree.preAligned.size();i++)
	    	alltree.preAligned[i]->aux_align = alltree.preAligned[i]->aln;
	}
	//cout << "Here: " << endl;
	
	subalign *taln;
	fprintf(logfp, "      - Calculating profiles ...\n");
	fprintf(logfp, "        ");
    	fflush(logfp);
	for(i=0;i<alltree.preAligned.size();i++) {
		alltree.preAligned[i]->aux_align->ss = alltree.preAligned[i]->aln->ss;
		taln = alltree.preAligned[i]->aux_align;
		taln->prof(1);
		//prof_freq allocation 
		taln->prof_freq = dmatrix(taln->prof_len, 20);
		// use_ss_freq:
		// 0 - use blosum62 matrix to get pseudocount
		// 1 - use aa_pair1[0] and aa_bg1[0] to get pseudocount
		// 2 - use aa_pair1 and aa_bg1 and ss->alphabet1 to get pseudocount
		if(use_ss_freq==0) {
			taln->pseudoCounts(taln->prof_effn, taln->prof_nef, taln->prof_len, taln->prof_freq);
		}
		else if(use_ss_freq==1) {
			taln->pseudoCounts(taln->prof_effn, taln->prof_nef, taln->prof_len, taln->prof_freq, params.aa_pair1[0], params.aa_bg1[0]);
		}
		else if(use_ss_freq==2) {
			if(psipred_env_number==9) taln->pseudoCounts(taln->prof_effn, taln->prof_nef, taln->prof_len, taln->prof_freq, params.aa_pair1, params.aa_bg1, taln->ss->alphabet1);
			if(psipred_env_number==3) taln->pseudoCounts(taln->prof_effn, taln->prof_nef, taln->prof_len, taln->prof_freq, params.aa_pair1, params.aa_bg1, taln->ss->sstype);
		}
		taln->log_pseudoCounts();
		
		//alltree.preAligned[i]->aux_align->get_prof_freq(use_ss_freq);
		//alltree.preAligned[i]->aux_align->get_prof_map_ss(alltree.preAligned[i]->aln->aname[0]);
		//a1->get_prof_alphabet1();
		//a1->get_prof_freq(0);
		//alltree.preAligned[i]->aux_align->get_score_bg_mine(params.aa_loop, use_ss_freq);
		//fprintf(logfp, "             Repres. sequence %d\r", i+1);
		fprintf(logfp, "*");
    	    	fflush(logfp);
	}
	fprintf(logfp, "\n");
	//cout << "Here: " << endl;
		
	if(debug_here>11) { cout << "Before consistency" << endl; }
	//cout << "Here: " << endl;
	fprintf(logfp, "      - Making consistency scoring fuction ...\n");
	fprintf(logfp, "        ");
    	fflush(logfp);
        alltree.profileConsistency_psipred(&params, dist_matrix_preAligned, 1, use_homologs);

	if(debug_here>11) { cout << "After consistency" << endl; }
	//cout << "Here: " << endl;

	fprintf(logfp, "      - Making progressive alignments ...\n");
    	fflush(logfp);
	alltree.computeConsistencyAlignment(alltree.root);
	//cout << "Here: " << endl;

	output_alignment();
	fprintf(logfp, "\nPROMALS is now finished\n");
    	fflush(logfp);
	fclose(logfp);
}

void multiple::alignDivergent_psipred_sum_of_pairs(int use_homologs) {

	int i, j;

	// right now, just the option of multim - for probablistic consistency 
        hmm_psipred_parameters params(psipred_env_number);
        params.read_parameters(psipred_parameter_file, psipred_env_number, 1);

	//cout << "Here: " << endl;

	//cout << psipred_dir << endl;

	for(i=0;i<alltree.preAligned.size();i++) {
		// psipred_dir is global variable
		alltree.preAligned[i]->aln->select_representative();
		// DO NOT run psipred from scratch
		//alltree.preAligned[i]->aln->get_ss_prof(psipred_dir, runpsipred_command);
	}
	//cout << "Here: " << endl;

	if(use_homologs==1) {
	     for(i=0;i<alltree.preAligned.size();i++)
		alltree.preAligned[i]->aux_align = alltree.preAligned[i]->similarSet->purge_align_one_seq_name(alltree.preAligned[i]->aln->aname[0]);
	}
	else if(use_homologs==2) {
	    for(i=0;i<alltree.preAligned.size();i++) {
		alltree.preAligned[i]->aux_align = get_blastpgp_alignment(alltree.preAligned[i]->aln->aname[0], alltree.preAligned[i]->aln->aname[0], alltree.preAligned[i]->aln->aseq[0]);
		alltree.preAligned[i]->aln->get_ss_prof1(blast_dir, alltree.preAligned[i]->aln->aname[0], runpsipred1_command);

		if(!alltree.preAligned[i]->aux_align) {
			cout << "reading blastpgp alignment failed " << alltree.preAligned[i]->aln->aname[0] << endl;
			alltree.preAligned[i]->aux_align = alltree.preAligned[i]->aln;
		}
	    }
	}
	else {
	     for(i=0;i<alltree.preAligned.size();i++)
	    	alltree.preAligned[i]->aux_align = alltree.preAligned[i]->aln;
	}
	//cout << "Here: " << endl;
	
	for(i=0;i<alltree.preAligned.size();i++) {
		alltree.preAligned[i]->aux_align->ss = alltree.preAligned[i]->aln->ss;
		//alltree.preAligned[i]->aux_align->ss->print_ss_info();
		alltree.preAligned[i]->aux_align->prof(1);
		//alltree.preAligned[i]->aux_align->get_prof_freq(use_ss_freq);
		//alltree.preAligned[i]->aux_align->get_prof_map_ss(alltree.preAligned[i]->aln->aname[0]);
		//a1->get_prof_alphabet1();
		//a1->get_prof_freq(0);
		//alltree.preAligned[i]->aux_align->get_score_bg_mine(params.aa_loop, use_ss_freq);
	}
	//cout << "Here: " << endl;
		
	if(debug_here>11) { cout << "Before consistency" << endl; }
	//cout << "Here: " << endl;
        alltree.profileConsistency_psipred_sum_of_pairs(&params, dist_matrix_preAligned, 1, use_homologs);

	if(debug_here>11) { cout << "After consistency" << endl; }
	//cout << "Here: " << endl;

	alltree.computeConsistencyAlignment(alltree.root);
	//cout << "Here: " << endl;

	output_alignment();
}

void multiple::output_alignment() {

        char outFileName[200];
        strcpy(outFileName, inputfileName);
        int tmpLen = strlen(outFileName);
        if( (outFileName[tmpLen-1]=='a')&&(outFileName[tmpLen-2]=='f')&&(outFileName[tmpLen-3]=='.')){
            outFileName[tmpLen-3] = '\0';
        }
        strcat(outFileName, ".promals.aln");
        if(!outFile.empty()) {
                strcpy(outFileName, outFile.c_str() );
        }
        cout << endl << "  output file Name: " << outFileName << endl << endl;
        alltree.printAlignmentFromAbs(alltree.root, outFileName);
        cout << "  program finished"  << endl << endl;

}

void multiple::align_profilehmm(int use_ss, char *ss_dir_name) {

	int i, j, k;
	subalign *x;

	// get the sequence profile and secondary structure profile for each tnode
	for(i=0;i<alltree.preAligned.size();i++) {
		x = alltree.preAligned[i]->similarSet;
		x->prof();
	        if( (use_ss==1)||(use_ss==0) ) {
       	        	x->get_prof_freq(use_ss, 1);
        	}
        	else if(use_ss==2) {
                	x->select_representative();
                	x->get_ss_prof(ss_dir_name, runpsipred_command);
                	//cout << x->repres_name << endl;
                	x->get_prof_map_ss(x->repres_name);
                	x->get_prof_alphabet1();
                	x->get_prof_freq(use_ss, 1);
                	//use_position_specific_regularizer = 1;
        	}
        	else {
                	cout << "Error: use_ss must be 0, 1, or 2" << endl;
                	exit(0);
        	}
		x->get_score_bg(use_ss);
	}

}

/*
// align divergent sequences to form pre-aligned groups, until distance between any two
// neighboring groups is larger than divergent_cutoff
void multiple::alignDivergent(float divergent_cutoff) {

	int i, j;
		
	// right now, just the option of multim - for probablistic consistency 
        hmm_parameters params(solv,ss,unaligned);
        params.read_parameters(parameter_file);
	if(debug_here>11) { cout << "Before consistency" << endl; }
        alltree.profileConsistency_multim(&params, dist_matrix_preAligned);

	if(debug_here>11) { cout << "After consistency" << endl; }

	alltree.computeConsistencyAlignment(alltree.root, divergent_cutoff/2);

	if(!alltree.root->aligned) return;

        char outFileName[200];
        strcpy(outFileName, inputfileName);
        int tmpLen = strlen(outFileName);
        if( (outFileName[tmpLen-1]=='a')&&(outFileName[tmpLen-2]=='f')&&(outFileName[tmpLen-3]=='.')){
            outFileName[tmpLen-3] = '\0';
        }
        strcat(outFileName, ".mummals.aln");
        if(!outFile.empty()) {
                strcpy(outFileName, outFile.c_str() );
        }
        cout << "  output file Name: " << outFileName << endl << endl;
        alltree.printAlignmentFromAbs(alltree.root, outFileName);
}

void alignExtremelyDivergent() {

	int i, j;
	
	alltree.obtainPreAligned(alltree.root);

}

void multiple::addSimilar() {

	int i,j;
}

*/

// store pre-aligned groups in the stopped nodes and select one representative from each group
void multiple::store_similar(tnode *r) {

	int i, j;

	if(!r->aligned) {
		store_similar(r->childL);
		store_similar(r->childR);
		return;
	}
	
	if(r->aln->nal == 1) return;
	
	// store the original subalign to "similarSet"
	r->similarSet = new subalign(*(r->aln));
	if(debug_here>11) r->similarSet->printali(80);

	// find a representative sequence for the group, make it the "aln"
	// right now, the representative is the longest sequence (excluding gaps)
	int tmp_index = 0;
	int tmp_count_aa = 0;
	int max_count_aa = 0;
	int find_target = 0;
	int target_index = -1;
	int target_count = 0;
	for(i=0;i<r->aln->nal;i++) {
		/*
		if(strlen(r->aln->aname[i])<=6) {
			target_index = i;
			find_target = 1;
			target_count++;
			tmp_index = target_index;
			tmp_count_aa = 0;
			for(j=0;j<r->aln->alilen;j++) {
				if(r->aln->aseq[i][j]!='-') tmp_count_aa++;
			}
			max_count_aa = tmp_count_aa;
			break;
		}
		*/
		tmp_count_aa = 0;
		for(j=0;j<r->aln->alilen;j++) {
			if(r->aln->aseq[i][j]!='-') tmp_count_aa++;
		}
		if(tmp_count_aa>max_count_aa) {
			max_count_aa = tmp_count_aa;
			tmp_index = i;
		}
	}
	/*if(target_count>1) { cout << "Not a good case: target count larger than one" << endl; exit(0);
	}
	if(find_target) {
		if(target_index!=tmp_index) {
			cout << "Not a good case" << endl;
			exit(0);
		}
	}
	*/
	char *tmp_name = new char [strlen(r->aln->aname[tmp_index])+1];
	char *tmp_seq = new char[max_count_aa+1];
	strcpy(tmp_name, r->aln->aname[tmp_index]);
	int tmp_array_index = 0;
	for(i=0;i<r->aln->alilen;i++) {
		if(r->aln->aseq[tmp_index][i]!='-') {
			tmp_seq[tmp_array_index] = r->aln->aseq[tmp_index][i];
			tmp_array_index++;
		}
	}
	tmp_seq[tmp_array_index] = '\0';
	r->aln = oneSeq2subalign(tmp_seq, tmp_name);
	if(debug_here>11) r->aln->printali(60);

}

// determine the position in the allseqs for any tnode in preAligned vector
//                p_seq
void multiple::map_allseqs_pos_to_tnode() {

	int i, j;
	subalign *tmpaln;

	if(debug_here>11) cout << alltree.preAligned.size() << "  " << allseqs.nseqs << endl;

	for(i=0;i<alltree.preAligned.size();i++) {
		tmpaln = alltree.preAligned[i]->aln;
		for(j=1;j<=allseqs.nseqs;j++) {
			//cout << tmpaln->aname[0] << "  " <<  allseqs.name[j].c_str() << endl;
			if(strcmp(tmpaln->aname[0], allseqs.name[j].c_str() )==0) {
				alltree.preAligned[i]->p_seq = j;
				if(debug_here>11) cout << i << " " << j << " " << tmpaln->aname[0] << endl;
				continue;
			}
		}
	}
}

// find N sequences with smallest distances to a sequence in the preAligned vector
// store them in a distance matrix
void multiple::get_distance_matrix_for_preAligned(int N_smallest) {

	int i, j;

	double tmp_dist_array[alltree.preAligned.size()+1];
	int auxilary[alltree.preAligned.size()+1];

	dist_matrix_preAligned = dmatrix(alltree.preAligned.size(), alltree.preAligned.size());
	
	for(i=0;i<alltree.preAligned.size();i++) {
		for(j=0;j<alltree.preAligned.size();j++) {
			tmp_dist_array[j+1] = 0;
			tmp_dist_array[j+1] = allseqs.distMat[alltree.preAligned[i]->p_seq][alltree.preAligned[j]->p_seq];
			auxilary[j+1] = j;
		}

		sort2(alltree.preAligned.size(), tmp_dist_array, auxilary);

		for(j=1;j<=((N_smallest+1<alltree.preAligned.size())?N_smallest+1:alltree.preAligned.size());j++) {
			dist_matrix_preAligned[i][auxilary[j]] = tmp_dist_array[j];
		}
	}

}

// merge a and b according to the sequence with seq_name
// b is added to a
subalign * merge_align_by_one_sequence(subalign *a, subalign *b, char *seq_name) {

	int i, j, k, l;

	int ia=-1, ib=-1;
	
	for(i=0;i<a->nal;i++) {
		if(strcmp(seq_name, a->aname[i])==0) { ia = i; break; }	
	}
	if(ia<0) { cout << "Name " << seq_name << " is not present in a" << endl; exit(0); }
	for(i=0;i<b->nal;i++) {
		if(strcmp(seq_name, b->aname[i])==0) { ib = i; break; }	
	}
	if(ib<0) { cout << "Name " << seq_name << " is not present in b" << endl; exit(0); }


	// determine the number of gaps of the linking sequence in a, and in b
	int ngapa=0, ngapb = 0;
	char *no_gap_seqa = cvector(a->alilen);
	char *no_gap_seqb = cvector(b->alilen);
	int tmp_index = 0;
	for(i=0;i<a->alilen;i++) {
		if(a->aseq[ia][i]!='-') {
			no_gap_seqa[tmp_index] = a->aseq[ia][i];
			tmp_index++;	
		}
		else { ngapa+=1; }
	} 	
	no_gap_seqa[tmp_index] = '\0';

	if(debug_here>11) cout << "ngapa: " << ngapa << endl;
	
	tmp_index = 0;
	for(i=0;i<b->alilen;i++) {
		if(b->aseq[ib][i]!='-') {
			no_gap_seqb[tmp_index] = b->aseq[ib][i];
			tmp_index++;	
		}
		else { ngapb+=1; }
	} 	
	no_gap_seqb[tmp_index] = '\0';

	if(debug_here>11) cout << "ngapb: " << ngapb << endl;
	if(debug_here>11) cout << "no_gap_seqb: " << no_gap_seqb << endl;

	int none_gap_len = strlen(no_gap_seqb);

	if(strcmp(no_gap_seqa, no_gap_seqb)!=0) {
		cout << "None-gapped sequences are different in subalign a and subalign b" << endl;
		cout << no_gap_seqa << endl;
		cout << no_gap_seqb << endl;
		exit(0);
	}

	// find the gap pattern arrays for a and b
	int *gap_pattern_a = ivector(none_gap_len);
	int *gap_pattern_b = ivector(none_gap_len);

	int tmp_gap_count = 0;
	tmp_index = 0;
	for(i=0;i<a->alilen;i++) {
		if(a->aseq[ia][i]!='-') {
			gap_pattern_a[tmp_index] = tmp_gap_count;
			if(debug>1) cout << tmp_index << " " << gap_pattern_a[tmp_index] << endl;
			tmp_index+=1;
			tmp_gap_count = 0;
		}
		else tmp_gap_count += 1;
	}
	gap_pattern_a[tmp_index] = tmp_gap_count;
	if(debug>1) cout << tmp_index << " " << gap_pattern_a[tmp_index] << endl;

	tmp_gap_count = 0;
	tmp_index = 0;
	for(i=0;i<b->alilen;i++) {
		if(b->aseq[ib][i]!='-') {
			gap_pattern_b[tmp_index] = tmp_gap_count;
			if(debug>1) cout << tmp_index << " " << gap_pattern_b[tmp_index] << endl;
			tmp_index+=1;
			tmp_gap_count = 0;
		}
		else tmp_gap_count += 1;
	}
	gap_pattern_b[tmp_index] = tmp_gap_count;
	if(debug>1) cout << tmp_index << " " << gap_pattern_b[tmp_index] << endl;

	// set up a new align
	subalign *new_aln = new subalign();
	new_aln->nal = a->nal + b->nal - 1;
	new_aln->alilen = strlen(no_gap_seqb) + ngapa + ngapb;
	new_aln->mnamelen = 0;
	for(i=0;i<a->nal;i++) { 
		if(strlen(a->aname[i])> new_aln->mnamelen) new_aln->mnamelen = strlen(a->aname[i]); 
	}
	for(i=0;i<b->nal;i++) { 
		if(strlen(b->aname[i])> new_aln->mnamelen) new_aln->mnamelen = strlen(b->aname[i]); 
	}
	if(debug_here>11) {
		cout << "new_aln: " << new_aln->nal << " " << new_aln->alilen << " " <<  new_aln->mnamelen << endl;
	}

	new_aln->aseq = cmatrix(new_aln->nal, new_aln->alilen+1);
	new_aln->aname = cmatrix(new_aln->nal, new_aln->mnamelen+1);

	for(i=0;i<a->nal;i++) {
		strcpy(new_aln->aname[i], a->aname[i]);
	}
	for(i=0;i<b->nal;i++) {
		if(i==ib) continue;
		if(i<ib) strcpy(new_aln->aname[i+a->nal], b->aname[i]);
		if(i>ib) strcpy(new_aln->aname[i+a->nal-1], b->aname[i]);
	}

	// generate new sequences
	// atmp_len: temporary length of the new sequences
	// tmp_index: the index of non-gapped residues in the linking sequence
	int atmp_len = 0;
	tmp_index = 0;
	for(i=0;i<a->nal;i++) {
		if(i!=ia) continue;
		for(j=0;j<a->alilen;j++) {
			if(a->aseq[i][j]!='-') {
				for(l=1;l<=gap_pattern_b[tmp_index];l++) {
					for(k=0;k<a->nal;k++) {
						new_aln->aseq[k][atmp_len] = '-';
					}
					atmp_len++;
				}
				for(k=0;k<a->nal;k++) {
					new_aln->aseq[k][atmp_len] = a->aseq[k][j];
				}
				atmp_len++;
				tmp_index+=1;
			}
			else {
				for(k=0;k<a->nal;k++) {
					new_aln->aseq[k][atmp_len] = a->aseq[k][j];
				}
				atmp_len++;
			}
		}
		// gaps at the C-terminal
		for(l=1;l<=gap_pattern_b[tmp_index];l++) {
			for(k=0;k<a->nal;k++) {
				new_aln->aseq[k][atmp_len] = '-';
			}
			atmp_len++;
		}
		// closing the sequences
		for(k=0;k<a->nal;k++) {
			new_aln->aseq[k][atmp_len] = '\0';
		}
	}

	atmp_len = 0;
	tmp_index = 0;
	for(i=0;i<b->nal;i++) {
		if(i!=ib) continue;
		for(j=0;j<b->alilen;j++) {
			if(b->aseq[i][j]!='-') {
				for(l=1;l<=gap_pattern_a[tmp_index];l++) {
					for(k=0;k<b->nal;k++) {
						// skipping the linking sequence
						if(k==ib) continue;
						if(k<ib) new_aln->aseq[k+a->nal][atmp_len] = '-';
						if(k>ib) new_aln->aseq[k+a->nal-1][atmp_len] = '-';
					}
					atmp_len++;
				}
				for(k=0;k<b->nal;k++) {
					if(k==ib) continue;
					if(k<ib) new_aln->aseq[k+a->nal][atmp_len] = b->aseq[k][j];
					if(k>ib) new_aln->aseq[k+a->nal-1][atmp_len] = b->aseq[k][j];
				}
				atmp_len++;
				tmp_index+=1;
			}
			else {
				for(k=0;k<b->nal;k++) {
					if(k==ib) continue;
					if(k<ib) new_aln->aseq[k+a->nal][atmp_len] = b->aseq[k][j];
					if(k>ib) new_aln->aseq[k+a->nal-1][atmp_len] = b->aseq[k][j];
				}
				atmp_len++;
			}
		}
		for(l=1;l<=gap_pattern_a[tmp_index];l++) {
			for(k=0;k<b->nal;k++) {
				if(k==ib) continue;
				if(k<ib) new_aln->aseq[k+a->nal][atmp_len] = '-';
				if(k>ib) new_aln->aseq[k+a->nal-1][atmp_len] = '-';
			}
			atmp_len++;
		}
		for(k=0;k<b->nal;k++) {
			if(k==ib) continue;
			if(k<ib) new_aln->aseq[k+a->nal][atmp_len] = '\0';
			if(k>ib) new_aln->aseq[k+a->nal-1][atmp_len] = '\0';
		}
	}
	
	if(debug>1) {
	a->printali(80);
	b->printali(80);
	new_aln->printali(80);
	}

	delete [] gap_pattern_a;
	delete [] gap_pattern_b;

	return new_aln;

}
char ssint2ss(int i) {

	if(i==1) return 'H';
	if(i==2) return 'E';
	if(i==3) return 'C';
}

