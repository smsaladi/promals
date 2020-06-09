#include "subalign.h"
#include "sequences.h"
#include "btree_template.h"
#include "progressiveAlignHMM.h"
#include "multiple.h"

static int debug_here = 1;

multiple::multiple(char *filename) {

	strcpy(inputfileName, filename);
	distance_cutoff_similar = 0.4;
	distance_cutoff_div = 0.7;
	N_small = 20;

	allseqs.readFasta(filename, 1);

	sequences allseqs1(allseqs);

	// k-mer dist matrix calculation
	allseqs1.toDayhoff6();
	allseqs1.generateD6t(6);
	allseqs1.diffCountD2t(allseqs1.d6t[1], allseqs1.d6t[2]);
	allseqs1.d6t2DistMat(6);

	allseqs.distMat = allseqs1.distMat; // adapt the distance matrix

	// build k-mer dist based tree
	alltree.UPGMA(allseqs1.distMat, allseqs.seq, allseqs1.name, allseqs1.nseqs);
}

void multiple::set_distance_cutoff_similar(double dist_cutoff) {

	distance_cutoff_similar = dist_cutoff;

}

void multiple::alignSimilar() {

	int i, j;

	// progressively align similar sequences using general substitution matrix
	alltree.progressiveAlignHMM_FastStage(alltree.root, distance_cutoff_similar/2);

       // if the root is aligned, stop
       if(alltree.root->aligned)  {
               output_alignment();
               exit(0);
       }

	// store pre-aligned groups in the stopped nodes and select one representative from each group
	store_similar(alltree.root);

	alltree.obtainPreAligned(alltree.root);

	map_allseqs_pos_to_tnode();

	get_distance_matrix_for_preAligned(N_small);

	// output distance matrix for representatives
	if(debug>1) {
		for(i=0;i<alltree.preAligned.size();i++) {
			for(j=0;j<alltree.preAligned.size();j++) {
				cout << dist_matrix_preAligned[i][j] << " ";
			}
			cout << endl;
		}
	}

	if(debug_here>-1) {
	    cout << endl;
	    cout << "Aligning..." << endl;
	    cout << "  NUMBER OF SEQUENCES: " << allseqs.nseqs;
	    cout << "  NUMBER OF PRE-ALIGNED GROUPS: " << alltree.preAligned.size() << endl;
	    cout << endl;
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

        output_alignment();

}

void multiple::output_alignment() {


	alltree.computeConsistencyAlignment(alltree.root);
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
        cout << "  program finished"  << endl << endl;

}

void multiple::addSimilar() {

	int i,j;

	

}

// store pre-aligned groups in the stopped nodes and select one representative from each group
void multiple::store_similar(tnode *r) {

	int i, j;

	if(!r->aligned) {
		store_similar(r->childL);
		store_similar(r->childR);
		return;
	}
	
	if(r->aln->nal == 1) return;

	//cout << "r->aln->nal: " << r->aln->nal << endl;
	
	// store the original subalign to "similarSet"
	r->similarSet = new subalign(*(r->aln));
	//cout << "=============" << endl;
	if(debug_here>1) r->similarSet->printali(80);

	// find a representative sequence for the group, make it the "aln"
	// right now, the representative is the longest sequence (excluding gaps)
	int tmp_index = 0;
	int tmp_count_aa = 0;
	int max_count_aa = 0;
	for(i=0;i<r->aln->nal;i++) {
		tmp_count_aa = 0;
		for(j=0;j<r->aln->alilen;j++) {
			if(r->aln->aseq[i][j]!='-') tmp_count_aa++;
		}
		if(tmp_count_aa>max_count_aa) {
			max_count_aa = tmp_count_aa;
			tmp_index = i;
		}
	}
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
	if(debug_here>1) r->aln->printali(60);

}

// determine the position in the allseqs for any tnode in preAligned vector
//                p_seq
void multiple::map_allseqs_pos_to_tnode() {

	int i, j;
	subalign *tmpaln;

	if(debug_here>1) cout << alltree.preAligned.size() << "  " << allseqs.nseqs << endl;

	for(i=0;i<alltree.preAligned.size();i++) {
		tmpaln = alltree.preAligned[i]->aln;
		for(j=1;j<=allseqs.nseqs;j++) {
			//cout << tmpaln->aname[0] << "  " <<  allseqs.name[j].c_str() << endl;
			if(strcmp(tmpaln->aname[0], allseqs.name[j].c_str() )==0) {
				alltree.preAligned[i]->p_seq = j;
				if(debug_here>1) cout << i << " " << j << " " << tmpaln->aname[0] << endl;
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

	if(debug_here>1) cout << "ngapa: " << ngapa << endl;
	
	tmp_index = 0;
	for(i=0;i<b->alilen;i++) {
		if(b->aseq[ib][i]!='-') {
			no_gap_seqb[tmp_index] = b->aseq[ib][i];
			tmp_index++;	
		}
		else { ngapb+=1; }
	} 	
	no_gap_seqb[tmp_index] = '\0';

	if(debug_here>1) cout << "ngapb: " << ngapb << endl;
	if(debug_here>1) cout << "no_gap_seqb: " << no_gap_seqb << endl;

	int none_gap_len = strlen(no_gap_seqb);

	if(strcmp(no_gap_seqa, no_gap_seqb)!=0) {
		cout << "None-gapped sequences are different in subalign a and subalign b" << endl;
		cout << no_gap_seqa << endl;
		cout << no_gap_seqb << endl;
		//exit(0);
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
	if(debug_here>1) {
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

	return new_aln;

}
