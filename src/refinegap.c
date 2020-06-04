#include "consv1.h"
#include "mathfunc.h"

void refinegap(subalign *aln, double gap_thr_lowest, int use_hwt, int remove_frag, int extend_gap_region);

void adjust_gap(subalign *aln, int start, int end);

void delete_complete_gap_positions(subalign *aln);

/*
int main (int argc, char **argv) {

	int i,j,k;
	char outalnfile[200];
	subalign *aln = new subalign (argv[1]);
	aln->printali("tmp1.aln", 70);
	double gap_thr_lowest = 0.5;
	if(argc>2) {
		gap_thr_lowest = atof(argv[2]);
	}
	int use_hwt = 1;
	int remove_frag = 1;
	int extend_gap_region = 1;

	refinegap(aln, gap_thr_lowest, use_hwt, remove_frag, extend_gap_region);
	// refinegap(aln, 0.5, 1, 1, 1);

	delete_complete_gap_positions(aln);

	sprintf(outalnfile, "%s_%2.2lf_%d_%d.clustal", argv[1], gap_thr_lowest, use_hwt, extend_gap_region);
	aln->printali(outalnfile, 60);

	return 0;

}
*/

void refinegap(subalign *aln, double gap_thr_lowest, int use_hwt, int remove_frag, int extend_gap_region) {

   	int i,j,k;
	int nal, alilen;
	double gap_thr;
	int debug=0;
	// subalign *aln = new subalign (argv[1]);
	consv  *csv = new consv(*aln);
	csv->freq();
	//csv->h_freq();
	//csv->ic_freq();
	//cout << "========"<< endl;
	csv->h_weight_all();

	nal = csv->nal; alilen=csv->alilen;
	
	// calculate gap fraction
	double *gap_fraction = dvector(csv->alilen);
	double *gap_fraction1 = dvector(csv->alilen);
	if(!use_hwt) {
		for(i=1;i<=csv->alilen;i++) gap_fraction[i] = csv->gap_fraction[i];
	}
	else {
		double total_hwt = 0;
		for(i=1;i<=csv->nal;i++) {
			total_hwt += csv->hwt_all[i];
		}
		//cout << "total_hwt: "<< total_hwt << endl;
		for(i=1;i<=csv->alilen;i++) {
			gap_fraction[i] = 0;
			for(j=1;j<=csv->nal;j++) {
				if(csv->alignment[j][i] == 0) gap_fraction[i] += csv->hwt_all[j];
			}
			gap_fraction[i] /= total_hwt;
			gap_fraction1[i] = gap_fraction[i];
		}
	}
	
	sort(csv->alilen, gap_fraction);
	//for(i=1;i<=csv->alilen;i++) { cout << i << "\t" << gap_fraction[i] << endl; }
	// if(gap_fraction[csv->alilen/5]<0.1) gap_thr = 0.1;
	if(gap_fraction[csv->alilen/5]<gap_thr_lowest) gap_thr = gap_thr_lowest;
	else gap_thr = gap_fraction[csv->alilen/5];
	//cout << "gap_thr: " << gap_thr << endl;

	// number of effective positions where gap_fraction is less than 0.5
	int *num_eff_pos = ivector(csv->nal);
	int ave_num_eff_pos = 0;
	for(i=1;i<=csv->nal;i++) {
		num_eff_pos[i] = 0;
		for(j=1;j<=csv->alilen;j++) {
			if(csv->alignment[i][j] && (gap_fraction1[j]<0.5) ) {
				num_eff_pos[i]++;
			}
		}
		ave_num_eff_pos += num_eff_pos[i];
	}
	ave_num_eff_pos /= csv->nal;
	//for(i=1;i<=csv->nal;i++) { cout << i << "\t" << num_eff_pos[i] << endl; }

	/* if(remove_frag) {
	   if(!use_hwt) {
		for(i=1;i<=csv->alilen;i++) gap_fraction[i] = csv->gap_fraction[i];
	   }
	   else {
		double total_hwt = 0;
		for(i=1;i<=csv->nal;i++) {
			total_hwt += csv->hwt_all[i];
		}
		cout << "total_hwt: "<< total_hwt << endl;
		for(i=1;i<=csv->alilen;i++) {
			gap_fraction[i] = 0;
			for(j=1;j<=csv->nal;j++) {
				if(csv->alignment[j][i] == 0) gap_fraction[i] += csv->hwt_all[j];
			}
			gap_fraction[i] /= total_hwt;
			gap_fraction1[i] = gap_fraction[i];
		}
	   }
	} */

	for(i=1;i<=csv->alilen;i++) gap_fraction[i] = gap_fraction1[i];
	int *gap_flag = ivector(csv->alilen);
	for(i=1;i<=csv->alilen;i++) {
		gap_flag[i] = 0;
		if(csv->gap_fraction[i]>gap_thr) {
			gap_flag[i] = 1;
		}
		if(extend_gap_region) {
			if(i>1) if(csv->gap_fraction[i-1]>gap_thr) gap_flag[i] = 1;
			if(i<csv->alilen) if(csv->gap_fraction[i+1]> gap_thr) gap_flag[i] = 1;
		}
	}
	//for(i=1;i<=csv->alilen;i++) cout << gap_flag[i]; cout << endl;

	int starts[alilen], ends[alilen];
	int start_count=0, end_count=0;
	for(i=1;i<=alilen;i++) {
		if(i==1) {
			if(gap_flag[i]==1) {
				starts[start_count] = i;
				start_count++;
			}
		}
		if(i==alilen) {
			if(gap_flag[i] ==1) {
				ends[end_count] = i;
				end_count++;
			}
			continue;
		}
		if( (!gap_flag[i]) && (gap_flag[i+1]) ){
			starts[start_count] = i+1;
			// cout << "start i: " << i << endl;
			start_count++;
		}
		if( (!gap_flag[i+1]) && (gap_flag[i]) ) {
			ends[end_count] = i;
			// cout << "end i: " << i << endl;
			end_count++;
		}
	}
	if(start_count!=end_count) 
	 	cout << "start_count: " << start_count << "\t endcount " << end_count << endl;
	if(debug>1) for(i=0;i<start_count;i++) {
		cout << i << "\t" << starts[i] << "\t" << ends[i] << endl;
	}

	for(i=0;i<start_count;i++) {
		//cout << starts[i] << " " << ends[i] << endl;
		adjust_gap(aln, starts[i], ends[i]);
	}		
	// aln->printali("tmp.aln", 60);

	delete csv;
	delete [] gap_fraction;
	delete [] gap_fraction1;
	delete [] gap_flag;	
	delete [] num_eff_pos;
			
}

void adjust_gap(subalign *aln, int start, int end) {
	
	int i,j;
	int nal = aln->nal;
	int alilen = aln->alilen;
	int gap_length = end - start +1;
	int tmparray[gap_length+1];
	char aarray[gap_length+1];
	int num_letters_in_gap;
	int **alignment = aln->alignment;
	char **aseq = aln->aseq;

	for(i=1;i<=nal;i++) {
		num_letters_in_gap = 0;
		for(j=start;j<=end;j++) {
			if(aln->alignment[i][j]!=0) {
				num_letters_in_gap++;
				tmparray[num_letters_in_gap-1] = aln->alignment[i][j];
				aarray[num_letters_in_gap-1] = aln->aseq[i-1][j-1];
			}
			//cout << "num_letters_in_gap-1: " << tmparray[num_letters_in_gap-1] << " |\n";
		}
		//cout << endl;
		if(num_letters_in_gap==0) {
			//cout << i<< endl;
			continue;
		}
	    if(start==1) { // push to the right
		for(j=start;j<=end-num_letters_in_gap;j++) {
			alignment[i][j] = 0;
			aseq[i-1][j-1] = '-';
		}
		for(j=end-num_letters_in_gap+1;j<=end;j++) {
			alignment[i][j] = tmparray[j-end+num_letters_in_gap-1]+00000;
			aseq[i-1][j-1] = aarray[j-end+num_letters_in_gap-1];
		}
	    }
	    else if(end==alilen) { // push to the left
		for(j=start;j<=start+num_letters_in_gap-1;j++) {
			alignment[i][j] = tmparray[j-start]+00000;
			aseq[i-1][j-1] = aarray[j-start];
		}
		for(j=start+num_letters_in_gap;j<=end;j++) {
			alignment[i][j] = 0;
			aseq[i-1][j-1] = '-';
		}
	    } else {
		for(j=start;j<=start+num_letters_in_gap/2-1;j++) {
			alignment[i][j] = tmparray[j-start]+00000;
			aseq[i-1][j-1] = aarray[j-start];
		}
		for(j=start+num_letters_in_gap/2;j<=end-(num_letters_in_gap-num_letters_in_gap/2);j++) {
			alignment[i][j] = 0;
			aseq[i-1][j-1] = '-';
		}
		for(j=end-(num_letters_in_gap-num_letters_in_gap/2)+1;j<=end;j++) {
			alignment[i][j] = tmparray[j-end+num_letters_in_gap-1]+00000;
			aseq[i-1][j-1] = aarray[j-end+num_letters_in_gap-1];
		}
	    }
		//for(j=start;j<=end;j++) { aseq[i-1][j-1] = am[alignment[i][j]]; }	
	}

}

void delete_complete_gap_positions(subalign *aln) {

    int i, j, k;

    int nal=aln->nal;
    int alilen = aln->alilen;

	//cout << nal << " " << alilen << endl;
    int all_gap_flag;

    char **tmp_seq = new char *[nal];
    for(i=0;i<nal;i++) tmp_seq[i] = new char[alilen+1];
	//cout << nal << " " << alilen << endl;


    int tmp_len = 0;
    for(i=0;i<alilen;i++) {
	all_gap_flag = 1;
	for(j=0;j<nal;j++) {
	    if(aln->aseq[j][i]!='-') all_gap_flag = 0;
	}
	if(all_gap_flag) continue;
	for(j=0;j<nal;j++) {
	    tmp_seq[j][tmp_len] = aln->aseq[j][i];
	}
	tmp_len += 1;
    }
    for(j=0;j<nal;j++) {
	tmp_seq[j][tmp_len] = '\0';
    }

    for(j=0;j<nal;j++) {
	strcpy(aln->aseq[j], tmp_seq[j]);
    }
    aln->alilen = tmp_len;
	
    free_cmatrix(tmp_seq, nal-1, alilen);

}
	

void treat_single_residues(subalign *aln, double position_fraction) {

	int i, j, k;
	int nal = aln->nal;
	int alilen = aln->alilen;
	int gap_count;

	int **sr_record = imatrix(nal, alilen);
	double *gap_fraction = dvector(alilen);

	for(j=1;j<=alilen;j++) {
		gap_count = 0;
		for(i=1;i<=nal;i++) {
			if(j==1) { if(aln->alignment[i][j+1] == 0) { sr_record[i][j] = 1; } }
			else if(j==alilen) { if(aln->alignment[i][j-1]==0) {sr_record[i][j] = 1; } }
			else if((aln->alignment[i][j+1] == 0)&&(aln->alignment[i][j-1]==0)) {sr_record[i][j] = 1; }

			if(aln->alignment[i][j]==0) gap_count++;
		}
		gap_fraction[j] = 1.0 * gap_count/nal;
	}
	int sr_count = 0;
	int n_pos, c_pos, t_pos;
	for(j=1;j<=alilen;j++) {
		sr_count = 0;
		// find out how many residues are single residues in a position
		for(i=1;i<=nal;i++) {
			if(sr_record[i][j] == 1) sr_count++;
		}
		// if the fraction of single residues exceeds position_fraction; do not consider them
		if(1.0 * sr_count/nal > position_fraction) {
			for(i=1;i<=nal;i++) sr_record[i][j] = 0;
		}
		// otherwise make the corrections --> find N, C adjacent residues, switch with the gap
		else {
			for(i=1;i<=nal;i++) {
				if(sr_record[i][j] == 1) {
					n_pos = c_pos = 0;
					for(k=1;k<alilen;k++) {
						if(j+k>alilen) break;
						if(aln->alignment[i][j+k]!=0) {break;} 
					}
					c_pos = j+k-1;
					for(k=1;k<alilen;k++) {
						if(j-k<0) break;
						if(aln->alignment[i][j-k]!=0) {break;}
					}
					n_pos = j-k+1;
					if(j==1) t_pos = c_pos;
					else if(j==alilen) t_pos = n_pos;
					else if(gap_fraction[n_pos]>gap_fraction[c_pos]) t_pos = c_pos;
					else t_pos = n_pos;
					aln->alignment[i][t_pos] = aln->alignment[i][j];
					aln->alignment[i][j] = 0;
					aln->aseq[i-1][t_pos-1] = aln->aseq[i-1][j-1];
					aln->aseq[i-1][j-1] = '-';
				}
			}
		}
	}

	free_imatrix(sr_record, nal, alilen);
	delete [] gap_fraction;
}


