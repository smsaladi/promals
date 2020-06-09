#include "refine.h"
#include <algorithm>

static int debug_here = 1;
void GetLower(string& str);

alignrefine::alignrefine(subalign *x, int NC_added_gaps) {

	int i, j;

	a = x;
	nal = a->nal;
	alilen = a->alilen;

	// add 50 gaps to both N and C terminal of the sequences
	aseq = new char * [ nal];
	for(i=0;i<nal;i++) aseq[i] = new char [alilen+NC_added_gaps+NC_added_gaps+1];
	for(i=0;i<nal;i++) {
		for(j=0;j<NC_added_gaps;j++) aseq[i][j] = '-';
		for(j=NC_added_gaps;j<NC_added_gaps+alilen;j++) aseq[i][j] = a->aseq[i][j-NC_added_gaps];
		for(j=NC_added_gaps+alilen;j<NC_added_gaps+NC_added_gaps+alilen;j++) aseq[i][j] = '-';
		aseq[i][alilen+NC_added_gaps+NC_added_gaps] = '\0';
		if(debug_here>1)cout << aseq[i] << endl;
	}
	alilen += (NC_added_gaps+NC_added_gaps);

	assign_gi();

	single_fraction = 0.5;

	sw = new double[nal+1];
	aaw = new double[alilen+1];
	sdw = new double[alilen+1];
	gappy_threshold = 0.1;

	CORE_POS_NUMBER = 3;
}

void alignrefine::set_gappy_threshold(double g) {

	gappy_threshold = g;

}

void alignrefine::set_CORE_POS_NUMBER(int n) {

	CORE_POS_NUMBER = n;

}

void alignrefine::assign_gi() {

	int i, j, k;

	gi = new gapinfo *[nal+1];
	for(i=1;i<=nal;i++) {
		gi[i] = new gapinfo[alilen+1];
	}

	for(i=1;i<=nal;i++) {
		for(j=1;j<=alilen;j++) {
			gi[i][j].aa = aseq[i-1][j-1];
			gi[i][j].lg = 0;
			gi[i][j].rg = 0;

			for(k=j-1;k>=1;k--) {
				if(aseq[i-1][k-1]=='-') gi[i][j].lg++;
				else break;
			}
			for(k=j+1;k<=alilen;k++) {
				if(aseq[i-1][k-1]=='-') gi[i][j].rg++;
				else break;
			}
			gi[i][j].single = 0;
			gi[i][j].doublet = 0;
			if( (gi[i][j].aa!='-')&& (gi[i][j].lg) && (gi[i][j].rg) ) {
				gi[i][j].single = 1;
			}
			if( (gi[i][j].aa!='-')&& (j==1) && (gi[i][j].rg) ) {
				gi[i][j].single = 1;
			}
			if( (gi[i][j].aa!='-')&& (j==alilen) && (gi[i][j].lg) ) {
				gi[i][j].single = 1;
			}
			if( (j<=alilen-2) && (gi[i][j].aa!='-')&&(aseq[i-1][j]!='-')&&(gi[i][j].lg)&&(aseq[i-1][j+1]=='-') ){
				gi[i][j].doublet = 1;
			}
			if( (j==1) && (j<=alilen-2) && (gi[i][j].aa!='-')&&(aseq[i-1][j]!='-')&&(aseq[i-1][j+1]=='-') ){
				gi[i][j].doublet = 1;
			}
			if( (j==alilen-1) && (gi[i][j].aa!='-')&&(aseq[i-1][j]!='-')&&(gi[i][j].lg) ){
				gi[i][j].doublet = 1;
			}

			if(debug_here > 1) cout << i << " " << j << " " << gi[i][j].aa << " " << gi[i][j].lg << " " << gi[i][j].rg << " " << gi[i][j].single << " " << gi[i][j].doublet << endl;
		}
	}
}

void alignrefine::print_gi_html() {

	string tmpstr;
	int i, j;

	for(i=1;i<=nal;i++) {
		tmpstr = "";
		for(j=1;j<=alilen;j++) {
			if(gi[i][j].single) {
				tmpstr += "<font color=red>";
				tmpstr += gi[i][j].aa;
				tmpstr += "</font>";
			}
			else if(gi[i][j].doublet) {
				tmpstr += "<font color=blue>";
				tmpstr += gi[i][j].aa;
				tmpstr += "</font>";
			}
			else {
				tmpstr += gi[i][j].aa;
			}
		}
		cout << tmpstr << endl;
	}
}

void alignrefine::assign_weight(double *weight) {
	
	int i;
	
	for(i=1;i<=nal;i++) { sw[i] = weight[i]; }

}

// weight at each position; gappy content;
// sw: sequence weight
// aaw: weight of amino acids at each position
// sdw: singlet/doublet weights for positions
// vaaw: a vector of weights of amino acids for positions
// w_aa: a vector of lists of 21 weight values for amino acids + gap character for positions
void alignrefine::calculate_weights() {

	int i, j, k;
	totalweight = 0;
	for(i=1;i<=nal;i++) {
		totalweight += sw[i];
		if(debug_here > 1) cout << "sw: " << sw[i] << endl;
	}
	if(debug_here > 1) cout << "totalweight: " << totalweight << endl;

	isgappy.clear();
	for(i=0;i<w_aa.size();i++) delete [] w_aa[i];
	w_aa.clear();
	vaaw.clear();
	vaaw.push_back(-1);
	w_aa.push_back(new double[1]);
	isgappy.push_back(-1);
	if(debug_here > 1) cout << "weight" << endl;
	for(j=1;j<=alilen;j++) {
		aaw[j] = 0;
		sdw[j] = 0;
		double *tmpw = new double[21];
		for(k=0;k<=20;k++) tmpw[k] = 0;
		w_aa.push_back(tmpw);
		for(i=1;i<=nal;i++) {
			//cout << "i: " << i <<" " << gi[i][j].aa<< endl;
			if(gi[i][j].aa!='-') {
				aaw[j] += sw[i];
				w_aa[j][am2num(gi[i][j].aa)]+=sw[i];
				if(gi[i][j].single || gi[i][j].doublet) {
					sdw[j] += sw[i];
				}
			}
			else {
				w_aa[j-1][0] += sw[i];
			}
		}

		if(aaw[j]<gappy_threshold) isgappy.push_back(1);
		else isgappy.push_back(0);
		vaaw.push_back(aaw[j]);
		if(debug_here > 1) cout << j << ": " << aaw[j] << " " << sdw[j] << " " << isgappy[j] << endl;
	}

}

void alignrefine::operation() {

	int i, j, k;

	for(i=1;i<=alilen;i++) {
		//cout << sdw[i] << " " << totalweight << " " << sdw[i]/totalweight << " " << a->aseq[0][i-1]<< endl;
		if(sdw[i]/totalweight > single_fraction) {
			operation_together(i);
		}
		else {
			operation_separate(i);
		}
	}
        for(i=1;i<=nal;i++) {
                for(j=1;j<=alilen;j++) {
                        a->aseq[i-1][j-1] = gi[i][j].aa; 
			a->alignment[i][j] = am2num(a->aseq[i-1][j-1]);
                } 
        } 

}

void alignrefine::operation_together(int p) {

	int i, j, k;

	int found;
	//cout << "Together: " << p << endl;

	// if the position or the next position contains a conserved 
	// amino acid, do not do the operation (frequency above 0.5)
	for(i=1;i<=20;i++) {
		if(w_aa[p][i]>=0.5) return;
		if(p<alilen) if(w_aa[p+1][i]>=0.5) return;
	}
	//cout << "here is operation together" << endl;

	// if p is the last position
	/*
	if(p==alilen) {
		for(i=1;i<=nal;i++) { shift_left(i, p); }
		return;
	}
	*/
		

	// it->first: lg   it->second: sum-of-weight with a certain lg
	map<int, double> gap_pattern;
	map<int, double>::const_iterator it;
	for(i=1;i<=nal;i++) {
		if(gi[i][p].aa=='-') continue; // ignore gaps
		if(gi[i][p].lg==0) continue; // ignore positions where the left position is not a gap
		if(p - gi[i][p].lg==1) continue; // ignore positions with a N-temrinal gap
		found = 0;
		for(it=gap_pattern.begin(); it != gap_pattern.end(); ++it) {
			// if there is a gap before position P
			// -PX
			if(it->first == gi[i][p].lg) {
				gap_pattern[gi[i][p].lg] += sw[i];
				found = 1;
				break;
			}
		}
		if(found == 0) gap_pattern[gi[i][p].lg] = sw[i];
	}
	if(debug_here>1) cout << "operation_together: " << p << endl;
	double max_left = 0;
	for(it=gap_pattern.begin(); it != gap_pattern.end(); ++it) {
		if(it->second>max_left) max_left = it->second;
		if(debug_here>1) cout << it->first << " " << it->second << endl;
	}
	// max_left: lg with the maximum sum-of-weight
	if(debug_here>1) cout << "max_left: " << max_left << endl;
	gap_pattern.clear();

	for(i=1;i<=nal;i++) {
		if(gi[i][p].aa=='-') continue; // ignore gaps
		if(p + gi[i][p].rg==alilen) continue; // ignore positions with a C-temrinal gap
		found = 0;
		for(it=gap_pattern.begin(); it != gap_pattern.end(); ++it) {
			// if there is a gap following position p
			// XP-
			if( (gi[i][p].rg ) && (it->first == gi[i][p].rg) ) {
				gap_pattern[gi[i][p].rg] += sw[i];
				found = 1;
				if(debug_here>1) cout << "single: " << i << " " << sw[i] << endl;
				break;
			}
			// otherwise if there is a gap following position p+1
			// XPA-
			else if( (gi[i][p+1].rg) && (it->first == gi[i][p+1].rg) ) {
				gap_pattern[gi[i][p+1].rg] += sw[i];
				found = 1;
				if(debug_here>1) cout << "doublet: " << i << " " << sw[i] << endl;
				break;
			}
		}
		if(found == 0) {
			if(gi[i][p].single) gap_pattern[gi[i][p].rg] = sw[i];
			if(gi[i][p].doublet) gap_pattern[gi[i][p+1].rg] = sw[i];
		}
	}
	double max_right = 0;
	for(it=gap_pattern.begin(); it != gap_pattern.end(); ++it) {
		if(it->second>max_right) max_right = it->second;
		if(debug_here>1) cout << it->first << " " << it->second << endl;
	}
	if(debug_here>1) cout << "max_right: " << max_right << endl;
	gap_pattern.clear();

	if(max_right <  max_left) {
		if(debug_here>1) cout << "operation_together_left: " << p << endl;
		for(i=1;i<=nal;i++) {
			if(gi[i][p].doublet) {
				shift_left(i, p);
				shift_left(i, p+1);
			}
			else {
				shift_left(i, p);
			}
		}
	}	
	else {
		if(debug_here>1) cout << "operation_together_right: " << p << endl;
		for(i=1;i<=nal;i++) {
			shift_right(i, p+1);
			shift_right(i, p);
			//if(gi[i][p].left!=0) if(gi[i][p].aa!='-') shift_right(i, p);
				
		}
	}		

}

void alignrefine::operation_separate(int p) {

	int i, j, k;

	for(i=1;i<=nal;i++) {
		if(gi[i][p].single) {
			// if the residue is highly conserved in the position, do not move
			//cout << i << " " << am2num(gi[i][p].aa) << " " << w_aa[p][am2num(gi[i][p].aa)] << endl;
			if(w_aa[p][am2num(gi[i][p].aa)]>=0.5) continue;

			if(debug_here>1) cout << "operation separate - single: " << i << " " << p;
			if(debug_here>1) cout << " " << gi[i][p].lg << " " << alilen << endl;
			if(p+gi[i][p].rg<=alilen) if( (aaw[p-gi[i][p].lg]<aaw[p]-0.5) && (aaw[p+gi[i][p].rg]<aaw[p]-0.5) ) continue;
			if(p-gi[i][p].lg==1) {                     // P-------AAAAA
				if(debug_here>1) cout << " right - Nterminal" << endl;
				shift_right(i, p);
				continue;
			}
			if(p+gi[i][p].rg==alilen) {               // AAAAA-------P
				if(debug_here>1) cout << " left - Cterminal" << endl;
				shift_left(i, p);
				continue;
			}
			if(aaw[p-gi[i][p].lg]>aaw[p+gi[i][p].rg]) {
				if(debug_here>1) cout << " left" << endl;
				shift_left(i, p);
			}
			else {
				if(debug_here>1) cout << " right" << endl;
				shift_right(i, p);
			}
		}
		if(gi[i][p].doublet) {
			// if the residue is highly conserved in the position, do not move
			if(w_aa[p][am2num(gi[i][p].aa)]>=0.5) continue;
			if(w_aa[p+1][am2num(gi[i][p+1].aa)]>=0.5) continue;

			//cout << "operation separate - single: " << i << " " << p;
			if(debug_here>1) cout << "operation separate - doublet: " << i << " " << p;
			if( (aaw[p-gi[i][p].lg]<aaw[p]-0.5) && (aaw[p+gi[i][p+1].rg]<aaw[p]-0.5) ) continue;
			if(p-gi[i][p].lg==1) {       // PP------AAAAAAA
				if(debug_here>1) cout << " right - Nterminal" << endl;
				shift_right(i, p+1);
				shift_right(i, p);
				continue;
			}
			if(p+1+gi[i][p+1].rg==alilen) { // AAAAAAAAA-------PP
				if(debug_here>1) cout << " left - Cterminal" << endl;
				shift_left(i, p);
				shift_left(i, p+1);
				continue;
			}
			if(aaw[p-gi[i][p].lg]>aaw[p+1+gi[i][p+1].rg]) {
				if(debug_here>1) cout << " left" << endl;
				shift_left(i, p);
				shift_left(i, p+1);
			}
			else {
				if(debug_here>1) cout << " right" << endl;
				shift_right(i, p+1);
				shift_right(i, p);
			}
		}
	}

}

//         3333333334
//         1234567890
//         A-----A--A
//          s    p t
void alignrefine::shift_left(int n, int p) {

	int i, j, k;
	
	if(gi[n][p].lg==0) return;
	if(gi[n][p].aa=='-') return;

	int s = p-gi[n][p].lg;
	int t = p+gi[n][p].rg;

	// for position s-1 
	if(s-1 > 0) {
		gi[n][s-1].rg = 0;
		gi[n][s-1].single = 0;
		if( (gi[n][s-1].lg) && (gi[n][s+1].aa=='-') ) {
			gi[n][s-1].doublet = 1;
		}
	}

	// for positon s
	gi[n][s].lg = 0;
	gi[n][s].aa = gi[n][p].aa;
	gi[n][s].single = 0;
	gi[n][s].doublet = 0;
	gi[n][s].rg = gi[n][p].rg + gi[n][p].lg;
	aaw[s] += sw[n];
	
	// for positon s+1 to positon p-1
	for(j=s+1;j<=p-1;j++) {
		gi[n][j].lg -= 1;
		gi[n][j].rg = gi[n][j].rg + 1 + gi[n][p].rg;
	}

	// for position p
	gi[n][p].aa = '-';
	gi[n][p].lg -= 1;
	gi[n][p].single = 0;
	gi[n][p].doublet = 0;
	aaw[p] -= sw[n];
	
	// for position p+1 to position t
	for(j=p+1;j<=t;j++) {
		gi[n][j].lg += gi[n][p].lg;
		gi[n][j].lg++;
	}
	
	// for position t+1
	if(t+1<=alilen) {
		gi[n][t+1].lg += gi[n][p].lg;
		gi[n][t+1].lg++;
		if( (gi[n][t].aa=='-') && (gi[n][t+1].rg) ) {
			gi[n][t+1].single = 1;
		}
		if(t+2<=alilen) if( (gi[n][t].aa=='-') && (gi[n][t+1].rg==0) && (gi[n][t+2].rg) ){
			gi[n][t+1].doublet = 1;
		}
	}
}

//         33333333344
//         12345678901
//         A-----A---A
//          s    p  t
void alignrefine::shift_right(int n, int p) {

	int i, j, k;

	if(gi[n][p].aa=='-') return;

	if(gi[n][p].rg==0) return;

	int s = p-gi[n][p].lg;
	int t = p+gi[n][p].rg;
	if(debug_here>1) cout << n << " " << p << " " << t << endl;

	// for position t+1;  aa, rg
	if(t+1<=alilen) {
		gi[n][t+1].lg = 0;
		if(gi[n][t+1].single) {
			gi[n][t+1].single = 0;
			sdw[t+1] -= sw[n];
		}
		if(gi[n][t+1].doublet) {
			gi[n][t+1].doublet = 0;
			sdw[t+1] -= sw[n];
		}
	}

	// for position t; 
	gi[n][t].aa = gi[n][p].aa;
	gi[n][t].single = 0;
	if(t+2<=alilen) if(gi[n][t+2].aa=='-') {
		gi[n][t].doublet = 1;
		sdw[t] += sw[n];
	}
	gi[n][t].rg = 0;
	gi[n][t].lg = gi[n][p].lg + gi[n][p].rg;
	aaw[t] += sw[n];

	// for position p+1 to position t-1; aa
	for(j=p+1;j<=t-1;j++) {
		gi[n][j].rg -= 1;
		gi[n][j].lg += 1;
		gi[n][j].lg += gi[n][p].lg;
	}
	
	// for position p
	gi[n][p].aa = '-';
	gi[n][p].rg -= 1;
	gi[n][p].lg;
	gi[n][p].single = gi[n][p].doublet = 0;
	aaw[p] -= sw[n];

	// for position s to position p-1;
	for(j=s;j<=p-1;j++) {
		gi[n][j].rg += 1;
		gi[n][j].rg += gi[n][p].rg;
	}
	
	// for position s-1
	gi[n][s-1].rg += 1;
	gi[n][s-1].rg += gi[n][p].rg;
	if(s-2>0) if( (gi[n][s].aa=='-') && (gi[n][s-2].aa=='-') ) gi[n][s-1].single = 1;
	gi[n][s-1].doublet = 0;

}

// henikoff weighting in the consv1 class
void alignrefine::get_sw() {

	int i, j;

	consv  *csv = new consv(*a);
	if(debug_here>1) cout << "========+++++++++" << endl;
	csv->freq();
	if(debug_here>1) cout << "========+++++++++" << endl;
	csv->h_weight_all();
	if(debug_here>1) cout << "========+++++++++" << endl;

	double total_hwt = 0;
        for(i=1;i<=csv->nal;i++) {
                total_hwt += csv->hwt_all[i];
		sw[i] = csv->hwt_all[i];
        }
	if(debug_here>1) cout << "========+++++++++" << endl;

	if(debug_here>1) cout << "total_hwt: " << total_hwt << endl;

}

void alignrefine::gi2aseq() {

	int i, j;

	for(i=1;i<=nal;i++) {
		for(j=1;j<=alilen;j++) {
			aseq[i-1][j-1] = gi[i][j].aa;
		}
	}

}

void alignrefine::treat_single_and_doublet(){

	int i, j;

	//cout << "+++++++++++++++++" << endl;
	assign_gi();
	if(debug_here>1) cout << "========" << endl;
	get_sw();
	if(debug_here>1) cout << "========" << endl;
	calculate_weights();
	if(debug_here>1) cout << "========" << endl;
	operation();
	if(debug_here>1) cout << "========" << endl;

}

// get rid of positions with all gaps
void alignrefine::refresh_alignment() {

	int i, j;

	if(debug_here>1) cout << "===========" << endl;
	delete_complete_gap_positions(a);
	if(debug_here>1) cout << "===========" << endl;

	//a->printali(80);

	alilen = a->alilen;
	aseq = a->aseq;

	delete [] aaw;
	delete [] sdw;
	isgappy.clear();
	for(i=1;i<=nal;i++) delete [] gi[i];
	delete [] gi;
	aaw = new double[alilen+1];
	sdw = new double[alilen+1];
	
	assign_gi();
	get_sw();
	calculate_weights();
}

// deal with gappy regions
// for n_terminal gappy regions
void alignrefine::deal_with_N() {

	int i, j;

	// deterimine where the N-terminal gappy region ends
	for(i=1;i<=alilen;i++) {
		if(isgappy[i]==0) break;
	}
	int N_bound = i-1;
	if(debug_here>1) cout << "N_bound: " << N_bound << endl;
	
	if(N_bound == 0) return;

	for(i=1;i<=nal;i++) {
		shift_gappy_right(1, N_bound, i);
	}
}
			
// for c terminal gappy region
void alignrefine::deal_with_C() {

	int i, j;

	// determine where the C_terminal gappy region begins
	for(i=alilen;i>=1;i--) {
		if(isgappy[i]==0) break;
	}
	int C_bound = i+1;
	if(debug_here>1) cout << "C_bound: " << C_bound << endl;
	if(debug_here>1) cout << "alilen: " << alilen << endl;
	if(C_bound>alilen) return;
	for(i=1;i<=nal;i++) {
		shift_gappy_left(C_bound, alilen, i);
	}
}

void alignrefine::shift_gappy_right(int gb, int ge, int n) {

	int i, j, k;
	
	assert(gb<=ge);
	deal_with_gappy_individual(gb, ge, n, -1);
	/*
	for(i=ge;i>=gb;i--) {
		if(aseq[n-1][i-1]!='-') {
			for(j=i+1;j<=alilen;j++) {
				if(aseq[n-1][j-1]!='-') break;
			}
			if(debug_here>1) cout << i << " " << j-1 << " " << n << endl;
			if(j-1==i) continue;
			aseq[n-1][j-2] = aseq[n-1][i-1];
			aseq[n-1][i-1] = '-';
		}
	}
	*/
}

void alignrefine::shift_gappy_left(int gb, int ge, int n) {

	int i, j, k;
	assert(gb <= ge);

	deal_with_gappy_individual(gb, ge, n, 1);
	/*
	for(i=gb;i<=ge;i++) {
		if(aseq[n-1][i-1]=='-') continue;
		for(j=i-1;j>=1;j--) {
			if(aseq[n-1][j-1]!='-') break;
		}
		if(j+1==i) continue;
		aseq[n-1][j]= aseq[n-1][i-1];
		aseq[n-1][i-1] = '-';
	}
	*/
}

void alignrefine::deal_with_gappy(int filter_length) {

	int i, j, k;
	int gb, ge;

	if(debug_here>1) cout << "=============" << endl;

	int original_alilen = alilen;

	// find the first gappy region bound by two core blocks
	for(i=1;i<=original_alilen;i++) {
		if(isgappy[i]==0) { break;}
	}
	for(j=i;j<=original_alilen;j++) {
		if(isgappy[j]) break;
	}
	for(k=j;k<=original_alilen;k++) {
		if(isgappy[k]==0) break;
	}
	
	gb = j;
 	ge = k-1;

	if(debug_here>-1) cout << "gb: " << gb << " ge: " << ge << endl;

	if(gb>ge) {
		if(debug_here>1) cout << "gb is larger than ge" << endl;
		return;
	}
	if(ge==alilen) return;

	deal_with_gappy(gb, ge, filter_length);


}

void alignrefine::deal_with_gappy(int gb, int ge, int filter_length) {

	int i, j, k;

  if(ge-gb+1< filter_length)  {

	for(i=1;i<=nal;i++) {
		cout << i << endl;
		deal_with_gappy_individual(gb, ge, i, 0);
	}

	cout << "==============" << endl;
	// find if there is any columns with entired gaps
	// if found, delete those positions
	int all_gap_begin=-1, all_gap_end=-1;
	for(i=gb;i<=ge;i++) {
		if(fabs(vaaw[i])<1e-6) {all_gap_begin=i; break;}
	}
	cout << "==============" << endl;
     if(all_gap_begin!=-1) {
	for(i=all_gap_begin+1;i<=ge;i++) {
		if(fabs(vaaw[i])>1e-6) {break;}
	}
	all_gap_end = i-1;
	if(debug_here>1) cout << "||||| " << all_gap_begin << "  " << all_gap_end << endl;
	delete_all_gap(all_gap_begin, all_gap_end);
	for(i=1;i<=nal;i++) {
		if(debug_here>1) cout << aseq[i-1] << endl;
	}
	alilen = alilen - all_gap_end + all_gap_begin  - 1;
	for(i=1;i<=alilen;i++) {
		if(debug_here>1) cout << i << "\t" << isgappy[i] << "\t" << vaaw[i] << endl;
	}
     }
	cout << "==============" << endl;
  }

	// go the the next gappy region
	// find the first gappy region bound by two core blocks
	for(i=gb;i<=alilen;i++) {
		if(isgappy[i]==0) { break;}
	}
	for(j=i;j<=alilen;j++) {
		if(isgappy[j]) break;
	}
	for(k=j;k<=alilen;k++) {
		if(isgappy[k]==0) break;
	}
	gb = j;
 	ge = k-1;

	if(gb>ge) return;

	if(debug_here>1) cout << "gb: " << gb << " ge: " << ge << " alilen: " << alilen << endl;
	if(ge == alilen) return;
	deal_with_gappy(gb, ge, filter_length);


}

// from b to e, are positions with all gaps
void alignrefine::delete_all_gap(int b, int e) {

	int i, j, k;

	string nstr, mstr, cstr, tstr;
	string allgap = "";
	for(i=b;i<=e;i++) allgap += "-";
	for(i=1;i<=nal;i++) {
		tstr = aseq[i-1];
		nstr = tstr.substr(0, b-1);
		mstr = tstr.substr(b-1, e-b+1);
		cstr = tstr.substr(e);
		if(debug_here>1) cout << nstr << endl;
		if(debug_here>1) cout << mstr << endl;
		if(debug_here>1) cout << cstr << endl;
		assert(mstr == allgap);
		tstr = nstr + cstr;
		strcpy(aseq[i-1], tstr.c_str());
		if(debug_here>1) cout << i << "\t" << aseq[i-1] << endl;
	}
	vaaw.erase(vaaw.begin()+b, vaaw.begin()+e+1);
	isgappy.erase(isgappy.begin()+b, isgappy.begin()+e+1);
	w_aa.erase(w_aa.begin()+b, w_aa.begin()+e+1);

	//for(i=0;i<a->nal;i++) delete [] a->aseq[i];
	//delete [] a->aseq;
	a->alilen = strlen(aseq[0]);
	for(i=0;i<nal;i++) {
		delete [] a->aseq[i];
		a->aseq[i] = new char [strlen(aseq[0])+1];
		strcpy(a->aseq[i], aseq[i]);
		delete [] a->alignment[i+1];
		a->alignment[i+1] = new int[alilen+1];
           	for(j=1;j<=alilen;j++) { a->alignment[i+1][j] = am2num(aseq[i][j-1]); }
	}
	//a->aseq = aseq;
	if(debug_here>1) a->printali(60);
	alilen = a->alilen;
}

// delete positions with all gaps in N and C terminals
void alignrefine::delete_NC_terminal_gaps() {
	
	int i, j;

        int allgapbegin = 1;
        int allgapend = 0;
        for(i=1;i<=alilen;i++) {
                if(vaaw[i]<1e-6) {
                        allgapend++;
                }
                else { break;   }
        }
        if(allgapend) delete_all_gap(allgapbegin, allgapend);

        allgapbegin = alilen+1;
        allgapend = alilen;
        for(i=alilen;i>=1;i--) {
                if(vaaw[i]<1e-6) { allgapbegin--; }
                else break;
        }
        if(debug_here>1) cout << "=============" << endl;
        if(debug_here>1) cout << allgapbegin << " " << allgapend << endl;
        if(allgapbegin<alilen+1) delete_all_gap(allgapbegin, allgapend);
}

// gb: gap beginning
// ge: gap ending
// n: sequence number
// lmr: left=-1; middle = 0; right = 1; shifting mode
void alignrefine::deal_with_gappy_individual(int gb, int ge, int n, int lmr) {

	int i, j, k;

	if(debug_here>-1) cout << "n: " << n << " lmr: " << lmr << endl;
	if(debug_here>-1) cout << aseq[n-1] << endl;

	// nstr: n_terminal fixed string
	// cstr: c terminal fixed string
	// mstr: middle string subject to change
	// tstr: total string = nstr+mstr+cstr
	// gappy_aa: amino acid string in the gappy region
	// naa, caa: amino acid string to the n or c terminal of the gappy region that 
	//      can be changed
	string nstr, cstr, gappy_aa, naa, caa, mstr, tstr;
	char *tmpstr;

	// lbm: number of letters in the gappy region
	// lmn: number of letters in the left movable region
	// rmn: number of letters in the right movable region
	// left_bound: start of the left movable region
	// right_bound: end of the right movable region
	int lbm, lmn, rmn, left_bound, right_bound;

	// left_core_pos: positions in the left movable region
	// right_core_pos: positions in the right movable region
	vector<int> left_core_pos, right_core_pos;

	// total string
	tstr = aseq[n-1];
	//cout << tstr << endl;
	
	// 1. get the amino acid string in the gappy region
	gappy_aa.clear();
	for(i=gb;i<=ge;i++) {
		if(aseq[n-1][i-1]!='-') gappy_aa += aseq[n-1][i-1];
	}
	if(gappy_aa.size()==0) {
		return;
	}
	lbm = gappy_aa.size();
	GetLower(gappy_aa);
	if(debug_here>-1) cout << "n: " << n << " " << lbm << endl;

	// 2. determine the left boundary
	lmn = 0;
	int left_core_count = 0;
	left_bound = gb;
	for(i=gb-1;i>0;i--) {
		if(debug_here>1) cout << tstr[i-1] << endl;
		if(tstr[i-1]=='-') {
			left_bound = i;
			lmn++;
			if(lmn==lbm) break;
		}
		else if(!isgappy[i]) {	
			left_core_count++;
			//if(left_core_count<=3) left_core_pos.push_back(i);
			if(left_core_count>CORE_POS_NUMBER) break;
		}
	}
	naa = "";
	for(i=left_bound;i<gb;i++) {
		if(tstr[i-1]!='-') naa += tstr[i-1];
		if(!isgappy[i]) left_core_pos.push_back(i);
	}
	//GetLower(naa);

	if(debug_here>1) cout << "left core position: " << endl;
	for(i=0;i<left_core_pos.size();i++) {
		if(debug_here>1) cout << left_core_pos[i] << endl;
	}
	
	// 3. determine the right boundary
	rmn = 0;
	int right_core_count = 0;
	right_bound = ge;
	for(i=ge+1;i<=alilen;i++) {
		if(tstr[i-1]=='-') {
			right_bound = i;
			rmn++;
			if(rmn==lbm) break;
		}
		else if(!isgappy[i]) {
			right_core_count++;
			//if(right_core_count<=3) right_core_pos.push_back(i);
			if(right_core_count>CORE_POS_NUMBER) break;
		}
	}
	caa = "";
	for(i=ge+1;i<=right_bound;i++) {
		if(tstr[i-1]!='-') caa += tstr[i-1];
		if(!isgappy[i]) right_core_pos.push_back(i);
	}
	//GetLower(caa);

	if(debug_here>-1) cout << "right core position: " << endl;
	for(i=0;i<right_core_pos.size();i++) {
		if(debug_here>1) cout << right_core_pos[i] << endl;
	}

	if(debug_here>-1) cout << "n: " << n << endl;
	if(debug_here>-1) cout << "left bound: " << left_bound << " " << nstr << " " << lmn << endl;
	if(debug_here>-1) cout << "right bound: " << right_bound << " " << cstr << " " << rmn << endl;
	if(debug_here>-1) cout << gappy_aa << " " << lbm << endl;

	nstr = tstr.substr(0, left_bound-1);
	cstr = tstr.substr(right_bound);
	mstr = tstr.substr(left_bound-1, right_bound-left_bound+1);
	//if(debug_here>1) cout << tstr << endl;
	if(debug_here>-1) cout << nstr << endl;
	if(debug_here>-1) cout << mstr << endl;
	if(debug_here>-1) cout << cstr << endl;
		
	// 4. if gappy_aa cannot be moved, push the amino acids aside depending on the position of the gaps
	if(lmn+rmn==0) {
		if(lmr==-1) mstr = add_left_gaps(gappy_aa, right_bound-left_bound+1);
		if(lmr==0) mstr = add_middle_gaps(gappy_aa, right_bound-left_bound+1);
		if(lmr==1) mstr = add_right_gaps(gappy_aa, right_bound-left_bound+1);
		tstr = nstr + mstr + cstr;
		if(debug_here>-1) cout << "tstr: " << tstr << endl;
		string oringal_str = aseq[n-1];
		cout << "**********"<<endl;
		strcpy(aseq[n-1], tstr.c_str());
		cout << "**********"<<endl;
		for(i=left_bound;i<=right_bound;i++) {
		cout << "**********"<<endl;
			cout << n << " " << i-1 << " " << oringal_str[i-1] << " " << sw[n] << endl;
			cout << am2num(oringal_str[i-1]) << endl;
			cout << nal << endl;
			cout <<  w_aa[i][am2num(oringal_str[i-1])] << endl;
			w_aa[i][am2num(oringal_str[i-1])] -= sw[n];
		cout << "**********"<<endl;
			if(am2num(oringal_str[i-1])) {aaw[i] -= sw[n]; vaaw[i] -= sw[n];}
		cout << "**********"<<endl;
			w_aa[i][am2num(aseq[n-1][i-1])] += sw[n];
		cout << "**********"<<endl;
			if(am2num(aseq[n-1][i-1])) {aaw[i] += sw[n]; vaaw[i] += sw[n];}
		cout << "**********"<<endl;
			if(vaaw[i]>=gappy_threshold) isgappy[i] = 0;
			else isgappy[i] = 1;
		}
		cout << "**********"<<endl;
		return;
	}

	// 5. if there is only one moving mode: lmn+rmn<=gappy_aa.size()
	string before_move, after_move;
	if(lmn+rmn<=gappy_aa.size()) {
		string gappy_aa_N = gappy_aa.substr(0, lmn);
		string gappy_aa_C = gappy_aa.substr(lbm-rmn, rmn);
		string gappy_aa_M = gappy_aa.substr(lmn, lbm-rmn-lmn);
		if(debug_here>1) cout << gappy_aa_N << "|" << gappy_aa_C << "|" << gappy_aa_M << endl;

		// judge if moving left is fine
		double oc_left, oc_right;
	     if(lmn>0) {
		before_move = tstr;
		after_move = nstr + naa + gappy_aa_N;
		oc_left = compare_score(before_move, after_move, left_core_pos, n);
	     }
	     if(rmn>0) {
		before_move = tstr;
		after_move = tstr.substr(0, ge) + gappy_aa_C + caa;
		oc_right = compare_score(before_move, after_move, right_core_pos, n);
	     }	
		mstr = gappy_aa_M;
		if(oc_left>=0) {
			nstr = nstr + naa + gappy_aa_N;
		}
		else {
			nstr = tstr.substr(0, gb-1);
			mstr = gappy_aa_N + mstr;
		}
		if(oc_right>=0) {
			cstr = gappy_aa_C + caa + cstr;
		}
		else {
			cstr = tstr.substr(ge);
			mstr = mstr + gappy_aa_C;
		}

		if(lmr==0) mstr = add_middle_gaps(mstr, ge-gb+1);
		if(lmr==-1) mstr = add_left_gaps(mstr, ge-gb+1);
		if(lmr==1) mstr = add_right_gaps(mstr, ge-gb+1);
		if(debug_here>1) cout << "tstr before: " << endl;
		if(debug_here>1) cout << tstr << endl;
		tstr = nstr + mstr + cstr;
		if(debug_here>1) cout << "tstr after: " << endl;
		if(debug_here>1) cout << tstr << endl;
		string oringal_str = aseq[n-1];
		strcpy(aseq[n-1], tstr.c_str());
		for(i=left_bound;i<=right_bound;i++) {
			w_aa[i][am2num(oringal_str[i-1])] -= sw[n];
			if(am2num(oringal_str[i-1])) {aaw[i] -= sw[n]; vaaw[i] -= sw[n];}
			w_aa[i][am2num(aseq[n-1][i-1])] += sw[n];
			if(am2num(aseq[n-1][i-1])) {aaw[i] += sw[n]; vaaw[i] += sw[n];}
			if(vaaw[i]>=gappy_threshold) isgappy[i] = 0;
			else isgappy[i] = 1;
		}

		return;
	}

	// 6. if there are multiple moving modes: lmn + rmn > lbm; choose the best one
	if(debug_here>1) cout << "=====++++++======" << endl;
	double oc_both, max_oc_both=-100;
	int move_index;
	string move_str;
	vector<int> both_core_pos;
	for(i=0;i<left_core_pos.size();i++) both_core_pos.push_back(left_core_pos[i]);
	for(i=0;i<right_core_pos.size();i++) both_core_pos.push_back(right_core_pos[i]);
	for(i=0;i<=lmn;i++) {
		if(debug_here>1) cout << "i: " << i << endl;
		if(lbm-i>rmn) continue;
		before_move = tstr;
		after_move = nstr;
		for(j=1;j<=gb-1-nstr.size()-i-naa.size();j++) after_move += '-';
		after_move = after_move + naa + gappy_aa.substr(0, i);
		for(j=1;j<=ge-gb+1;j++) after_move += '-';
		after_move += gappy_aa.substr(i);
		after_move += caa;
		for(j=1;j<=(nstr+mstr).size()-ge-caa.size()-(gappy_aa.substr(i)).size();j++) after_move += '-';
		after_move += cstr;
		if(debug_here>1) cout << before_move << endl;
		if(debug_here>1) cout << after_move<< endl;
		if(debug_here>1) cout << "=======================" << endl;
		if(debug_here>1) cout << endl;
		
		oc_both = compare_score(before_move, after_move, both_core_pos, n);
		if(debug_here>1) cout << "oc_both: " << oc_both << " " << max_oc_both << endl;
		if(oc_both>max_oc_both) {
			max_oc_both = oc_both;
			move_index = i;
			move_str = after_move;
			if(debug_here>1) cout << "here: " << endl;
		}
	}
	if(move_str==after_move) if(debug_here>1) cout << "+++++" << endl;
	if(debug_here>1) cout << after_move << endl;
	if(max_oc_both <0) return;
	if(debug_here>-1) cout << "tstr before: " << endl;
	if(debug_here>-1) cout << tstr << endl;
	if(debug_here>-1) cout << "tstr after: " << endl;
	if(debug_here>-1) cout << move_str << endl;
	string oringal_str = aseq[n-1];
	strcpy(aseq[n-1], move_str.c_str());
	for(i=left_bound;i<=right_bound;i++) {
		w_aa[i][am2num(oringal_str[i-1])] -= sw[n];
		if(am2num(oringal_str[i-1])) {aaw[i] -= sw[n]; vaaw[i] -= sw[n];}
		w_aa[i][am2num(aseq[n-1][i-1])] += sw[n];
		if(am2num(aseq[n-1][i-1])) {aaw[i] += sw[n]; vaaw[i] += sw[n];}
		if(vaaw[i]>=gappy_threshold) isgappy[i] = 0;
		else isgappy[i] = 1;
	}
	return;
}

// a1: sequence before move
// a2: sequence after move
// b: vector of core positions
// n: the n-th sequence
// return value: occupancy change
//    if resulting in significant frequency change (>0.5), return -1
double alignrefine::compare_score(string a1, string a2, vector<int> b, int n) {
	int i,j, k;

	vector <double> before_freq;
	vector <double> after_freq;
	double before_occupancy = 0;
	for(i=0;i<b.size();i++) {
		if(debug_here>1) cout << b[i] << " " << a1[b[i]-1]<<endl;
		if(a1[b[i]-1]!='-') {
			before_occupancy += vaaw[b[i]];
			before_freq.push_back(w_aa[b[i]][am2num(a1[b[i]-1])]);
		}
		else before_freq.push_back(0);
	}
	double after_occupancy = 0;
	double before_cummul_freq = 0;
	double after_cummul_freq = 0;
	double change_cummul_freq = 0;
	double change_freq = 0;
	for(i=0;i<b.size();i++) {
		if(a2[b[i]-1]!='-') {
			after_occupancy += vaaw[b[i]];
			after_freq.push_back(w_aa[b[i]][am2num(a2[b[i]-1])]);
		}
		else after_freq.push_back(0);
	}

	for(i=0;i<before_freq.size();i++) {
		if(debug_here>1) cout << "before-after freq: ";
		if(debug_here>1) cout << before_freq[i]-after_freq[i] << endl;
	}
	for(i=0;i<before_freq.size();i++) {
		if(before_freq[i]-after_freq[i]>0.5) return -1;
		change_freq += before_freq[i]-after_freq[i];
	}
	if(change_freq>0.5) return -1;

	if(debug_here>1) cout << before_occupancy << " " << after_occupancy << endl;

	return after_occupancy-before_occupancy;
}

string alignrefine::add_middle_gaps(string mstr, int length) {

	int i, j;
	
	string tmpstring = mstr.substr(0, mstr.size()/2);
	for(i=0;i<length-mstr.size();i++) {
		tmpstring += "-";
	}
	tmpstring += mstr.substr(mstr.size()/2);
	
	GetLower(mstr);

	return tmpstring;
}

string alignrefine::add_left_gaps(string mstr, int length) {

	int i, j;
	
	string tmpstring;
	for(i=0;i<length-mstr.size();i++) {
		tmpstring += "-";
	}
	tmpstring += mstr;

	GetLower(mstr);
	return tmpstring;
}


string alignrefine::add_right_gaps(string mstr, int length) {

	int i, j;
	
	string tmpstring;
	tmpstring = mstr;
	for(i=0;i<length-mstr.size();i++) {
		tmpstring += "-";
	}

	GetLower(mstr);
	return tmpstring;
}

void GetLower(string& str)
{
   //for(int i = 0; i < str.length(); i++) { str[i] = tolower(str[i]); }
   return;
}
