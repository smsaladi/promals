#include "all.h"
#include "hmm_profpair.h"

static int debug = 0;
hmm_profpair::hmm_profpair(subalign *x1, subalign *y1) {

	x = x1;
	y = y1;
	lenx = x->alilen;
	leny = y->alilen;

	gapopen_end = 0;

	delta = 0.019931; 
	epsilon = 0.79433; 
	tau = 0.19598;
}
	
	

void hmm_profpair::viterbi() {

	int i,j,k;
	int xi, yj, a;

	float mm, mxy, xym, xy, E;
	int maxalnlength = lenx + leny +1;
	int *reverse_path, *reverse_path1, *reverse_path2;

	mm = log(1-2*delta-tau);
	mxy = log(1-epsilon-tau);
	xym = log(delta);
	xy = log(epsilon);
	E = log(tau);

	gapopen_end = 1;

	fprintf(stdout, "tau: %f delta: %f epsilon: %f\n", tau, delta, epsilon);

	// initialize the matrices
	Vm = gmatrix<float>(lenx, leny);
	Vx = gmatrix<float>(lenx, leny);
	Vy = gmatrix<float>(lenx, leny);
	Tm = imatrix(lenx, leny);
	Tx = imatrix(lenx, leny);
	Ty = imatrix(lenx, leny);
	path = ivector(maxalnlength);
	path1 = ivector(maxalnlength);
	path2 = ivector(maxalnlength);
	reverse_path = ivector(maxalnlength);
	reverse_path1 = ivector(maxalnlength);
	reverse_path2 = ivector(maxalnlength);

	for(i=0;i<=lenx;i++) Vm[i][0] = Vx[i][0] = Vy[i][0] = -1000000;
	for(i=0;i<=leny;i++) Vm[0][i] = Vx[0][i] = Vy[0][i] = -1000000;
	Vm[0][0] = 0; //Vx[0][0] = Vy[0][0] = 0;
	for(i=1;i<=lenx;i++) {
		if(!gapopen_end) Vx[i][0] = max2(xy+Vm[i-1][0], xy+Vx[i-1][0], Tx[i][0]);
		// for end gap, waive the gap opening penalty: change
		else Vx[i][0] = max2(mxy+Vm[i-1][0], xy+Vx[i-1][0], Tx[i][0]);
		if(debug>1) cout << "Vx i 0: " << i << " " << Vx[i][0] << endl;
	}
	for(j=1;j<=leny;j++) {
		if(!gapopen_end) Vy[0][j] = max2(xy+Vm[0][j-1], xy+Vy[0][j-1], Ty[0][j]);
		else Vy[0][j] = max2(mxy+Vm[0][j-1], xy+Vy[0][j-1], Ty[0][j]);
		if(debug>1) cout << "Vy 0 j: " << j << " " << Vy[0][j] << endl;
	} 
	
	for(i=1;i<=lenx;i++) {
	   for(j=1;j<=leny;j++) {
		// testing
		// Vm[i][j] = log_odds_score(i,j) + max3(mm+Vm[i-1][j-1], mxy+Vx[i-1][j-1], mxy+Vy[i-1][j-1], Tm[i][j]);
		// Vx[i][j] = log_odds_x_d(i) + max2(mxy+Vm[i-1][j], xy+Vx[i-1][j], Tx[i][j]);
		// Vy[i][j] = log_odds_y_d(j) + max2(mxy+Vm[i][j-1], xy+Vy[i][j-1], Ty[i][j]);

		Vm[i][j] = log_odds_score(i,j) + max3(mm+Vm[i-1][j-1], mxy+Vx[i-1][j-1], mxy+Vy[i-1][j-1], Tm[i][j]);
		// for end gap, waive the gap opening penalty: change
		if(j==leny) {
		     if(!gapopen_end) Vx[i][j] = max2(xy+Vm[i-1][j], xy+Vx[i-1][j], Tx[i][j]);
		     else Vx[i][j] = max2(mxy+Vm[i-1][j], xy+Vx[i-1][j], Tx[i][j]);
		}	
		else Vx[i][j] = max2(mxy+Vm[i-1][j], xy+Vx[i-1][j], Tx[i][j]);
		// for end gap, waive the gap opening penalty: change
		if(i==lenx) {
		     if(!gapopen_end) Vy[i][j] = max2(xy+Vm[i][j-1], xy+Vy[i][j-1], Ty[i][j]);
		     else Vy[i][j] = max2(mxy+Vm[i][j-1], xy+Vy[i][j-1], Ty[i][j]);
		}
		else Vy[i][j] = max2(mxy+Vm[i][j-1], xy+Vy[i][j-1], Ty[i][j]);

		if(debug>1) fprintf(stdout, "%d %d %f %f %f %d %d %d %f\n", i, j, Vm[i][j], Vx[i][j], Vy[i][j], Tm[i][j], Tx[i][j], Ty[i][j], log_odds_score(i,j));
	   }
	}	

	VE = E + max3(Vm[lenx][leny], Vx[lenx][leny], Vy[lenx][leny], TE);

	// trace back
	xi = lenx;
	yj = leny;
	int mark_mxy;
	if(TE==1) {
		reverse_path[0] = 0;
		reverse_path1[0] = 0;
		reverse_path2[0] = 0;
		mark_mxy = 1;
		//xi--; yj--;
	}
	else if(TE==2) {
		reverse_path[0] = 1;
		reverse_path1[0] = 0;
		reverse_path2[0] = 1;
		mark_mxy = 2;
		//xi--;
	}
	else {
		reverse_path[0] = -1;
		reverse_path1[0] = 1;
		reverse_path2[0] = 0;
		mark_mxy = 3;
		//yj--;
	}
	i = 1;
	//fprintf(stdout, "=================\n");
	while( (xi!=0) || (yj!=0) ) {
		if(debug>1) fprintf(stdout, "%d %d    xi: %d yj: %d Tm %d Tx %d Ty %d\n", i, mark_mxy, xi, yj, Tm[xi][yj], Tx[xi][yj], Ty[xi][yj]);
		if(mark_mxy==1) {
			a = Tm[xi][yj];
			if(a==1) {
				//reverse_path[i+1] = 0;
				//xi--; yj--;
				mark_mxy = 1;
			}
			else if(a==2) {
				//reverse_path[i+1] = 1;
				//xi--;
				mark_mxy = 2;
			} 
			else if(a==3) {
				//reverse_path[i+1] = -1;
				//yj--;
				mark_mxy = 3;
			}
			reverse_path[i] = 0;
			reverse_path1[i] = 0; // 0 means amino acid
			reverse_path2[i] = 0;
			xi--; yj--;
			i++;
			continue;
		}
		if(mark_mxy==2) {
			a = Tx[xi][yj];
			if(a==1) {
				//reverse_path[i+1] = 0;
				// xi--; yj--;
				mark_mxy = 1;
			}
			else if(a==2) {
				//reverse_path[i+1] = 1;
				// xi--;
				mark_mxy = 2;
			}
			reverse_path[i] = 1;
			reverse_path1[i] = 0;
			reverse_path2[i] = 1; // 1 means gap
			xi--;
			i++;
			continue;
		}
		if(mark_mxy==3) {
			a = Ty[xi][yj];
			if(a==1) {
				//reverse_path[i+1] = 0;
				// xi--; yj--;
				mark_mxy = 1;
			}
			else if(a==2) {
				//reverse_path[i+1] = -1;
				// yj--;
				mark_mxy = 3;
			}
			reverse_path[i] = -1;
			reverse_path1[i] = 1; // 1 means gap
			reverse_path2[i] = 0;
			yj--;
			i++;
			continue;
		}
	}
	//fprintf(stdout, "=================\n");

	v_alilen = i-1;
	fprintf(stdout, "v_alilen: %d\n", v_alilen);
	for(i=1;i<=v_alilen;i++) {
		path[i] = reverse_path[v_alilen-i+1];
		path1[i] = reverse_path1[v_alilen-i+1];
		path2[i] = reverse_path2[v_alilen-i+1];
		if(debug>1) fprintf(stdout, "%d %d %d\n", i, reverse_path1[i], reverse_path2[i]);
	}
	if(debug>1) fprintf(stdout, "\n");
	
}
		
// this is correponds to sum-of-pairs of the log-odds RATIO scores of weighted amino acid pairs
float hmm_profpair::log_odds_score(int xi, int yj) {
	
	int i,j;

	float score = 0;
	
	for(i=1;i<=20;i++) {
	   for(j=1;j<=20;j++) {
		////score += (x->pseudoCnt[xi][i] * y->pseudoCnt[yj][j] * log(q_blosum62[i][j]/robinson_freq[i]/robinson_freq[j]) );
		score += (x->pseudoCnt[xi][i] * y->pseudoCnt[yj][j] * log_q_blosum62_ratio[i][j]);
		// testing
		// score += (x->pseudoCnt[xi][i] * y->pseudoCnt[yj][j] * log(q_blosum62[i][j]) );
	   }
	}

	// if(xi==lenx) if(yj==leny) { for(i=1;i<=20;i++) {
	//		fprintf(stdout, "*  %d %f  %f\n", i, x->pseudoCnt[xi][i], y->pseudoCnt[yj][i]); 
	//	}
	//}

	return score;
}
	
float hmm_profpair::max3(float a,float b, float c, int &t) {

	float max;
	
	if(a>=b) {
		if(a>=c) {
			max = a;
			t = 1;
		}
		else {
			max = c;
			t = 3;
		}
	}
	else {
		if(b>c) {
			max = b;
			t = 2;
		}
		else {
			max = c;
			t = 3;
		}
	}
	return max;
} 

	
float hmm_profpair::max2(float a, float b, int &t) {
	
	if(a>=b) {
		t= 1;
		return a;
	}
	else {
		t = 2;
		return b;
	}
}

void hmm_profpair::printHmmAlign(int blocksize) {

	int i,j,k;

	int nblocks;

	int max_name_len;

	if(x->mnamelen > y->mnamelen) max_name_len = x->mnamelen;
	else max_name_len =  y->mnamelen;

	fprintf(stdout, "%d %d\n", x->mnamelen, y->mnamelen);
	fprintf(stdout, "\n");
	fprintf(stdout, "The final alignment: \n\n");
	
	// number of blocks
	nblocks = (v_alilen-1)/blocksize + 1;

	int a = -1;
	int b = -1;
	int a_prev, b_prev, a_after, b_after;
	for(i=1; i<=nblocks; i++) {
		//fprintf(stdout, "=================\n\n");
		a_prev = a;
		for(j=1; j<=x->nal; j++) {
			cout << setiosflags(ios::left) << setw(max_name_len+3) << x->aname[j-1];
			for(k=(i-1)*blocksize+1; k<=i*blocksize; k++) {
				if(k>v_alilen) break;
				if(path[k]!=-1) {
					a++;
					cout << x->aseq[j-1][a];
				}
				else cout << "*";
			}
			fprintf(stdout, "\n");
			a_after = a;
			a = a_prev;
		}
		a = a_after;
		fprintf(stdout, "\n");

		b_prev = b;
		for(j=1; j<=y->nal; j++) {
			cout << setiosflags(ios::left) << setw(max_name_len+3) << y->aname[j-1];
			for(k=(i-1)*blocksize+1; k<=i*blocksize; k++) {
				if(k>v_alilen) break;
				if(path[k]!=1) {
					b++;
					cout << y->aseq[j-1][b];
				}
				else cout << "*";
			}
			fprintf(stdout, "\n");
			b_after = b;
			b = b_prev;
		}
		b = b_after;
		fprintf(stdout, "\n\n");
	}
}

// generate an alignment of the two alignments
subalign *hmm_profpair::productViterbiAlign() { 

	int i,j,k;

	//subalign *newAlign = new subalign();
	subalign *a = new subalign(); // newAlign;

	a->nal = x->nal + y->nal;
	a->alilen = v_alilen;
	if(x->mnamelen > y->mnamelen) a->mnamelen = x->mnamelen;
	else a->mnamelen = y->mnamelen;
	if(debug>1) cout << "a->mnamelen: " << a->mnamelen <<endl;
	if(debug>1) cout << "a->alilen: " << a->alilen <<endl;
	if(debug>1) cout << "a->nal: " << a->nal <<endl;
	a->aname = cmatrix(a->nal, a->mnamelen+1);
	a->aseq = cmatrix(a->nal, a->alilen+1);
	a->alignment = imatrix(a->nal, a->alilen);

	// copy the names
	for(i=0;i<x->nal;i++) {
		strcpy(a->aname[i], x->aname[i]);
	}
	for(i=x->nal;i<a->nal;i++) {
		strcpy(a->aname[i], y->aname[i - x->nal]);
	}

	if(debug>1) for(i=1;i<=v_alilen;i++) { cout << i << "\t" << path[i] << endl; }

	// derive the sequences from alignment path
	int xi=0, yi=0;
	for(i=0;i<a->alilen;i++) {
		for(j=0;j<x->nal;j++) {
			if(path[i+1]>=0) {
				a->aseq[j][i] = x->aseq[j][xi];
			}
			else {
				a->aseq[j][i] = '-';
			}
		}
		if(path[i+1]>=0) xi++;
		for(j=x->nal;j<a->nal;j++) {
			if(path[i+1]<=0) {
				a->aseq[j][i] = y->aseq[j - x->nal][yi];
			}
			else {
				a->aseq[j][i] = '-';
			}
		}
		if(path[i+1]<=0) yi++;
		if(debug>1) cout << "i: " << i<< endl;
	}
	for(i=0;i<a->nal;i++) {
		a->aseq[i][a->alilen]='\0';
	}

	// convert the letters to numbers
	a->convertAseq2Alignment();

	return a;
}

void hmm_profpair::forward() {

	int i,j,k;

	ScoreType mm, mxy, xym, xy, E;

	mm = LOG( 1-2*delta-tau );
	mxy = LOG( 1-epsilon-tau );
	xym = LOG( delta );
	xy = LOG( epsilon );
	E = LOG( tau );

	Fm = new ScoreType * [lenx+1];
	Fx = new ScoreType * [lenx+1];
	Fy = new ScoreType * [lenx+1];

	cout << "Forward algorithm:" << endl;

	for(i=0;i<=lenx;i++) {
		Fm[i] = new ScoreType [leny+1];
		Fx[i] = new ScoreType [leny+1];
		Fy[i] = new ScoreType [leny+1];
	}
	
	if(debug>1) fprintf(stdout, "===================\n");	
	Fm[0][0] = LOG_ONE; Fx[0][0] = LOG_ZERO; Fy[0][0] = LOG_ZERO;
	if(debug>1) fprintf(stdout, "===================\n");	
	for(i=0;i<=lenx;i++) {
		for(j=0;j<=leny;j++) {
			if( (i==0)&&(j==0) ) {; }
			//else if( (i==0)&&(j==1) ) {
			else if( (i==0)&&(j!=0) ) {
				Fm[i][j] = LOG_ZERO;
				Fx[i][j] = LOG_ZERO;
				//Fy[i][j] = xym * Fm[i][j-1] + xy * Fy[i][j-1];
				// waive gap open penalty for ends
				////Fy[i][j] = xy * Fm[i][j-1] + xy * Fy[i][j-1];
				////Fy[i][j] = Fy[i][j] * log_odds_y(j);
				Fy[i][j] = LOG_ADD(xy+Fm[i][j-1], xy+Fy[i][j-1])+log_odds_y(j);
			}
			//else if( (i==1)&&(j==0) ) {
			else if( (i!=0)&&(j==0) ) {
				Fm[i][j] = LOG_ZERO;
				//Fx[i][j] = xym * Fm[i-1][j] + xy * Fx[i-1][j];
				// waive gap open penalty for ends
				//Fx[i][j] = xy * Fm[i-1][j] + xy * Fx[i-1][j];
				//Fx[i][j] = Fx[i][j] * log_odds_x(i);
				Fx[i][j] = LOG_ADD(xy+Fm[i-1][j], xy+Fx[i-1][j])+log_odds_x(i);
				Fy[i][j] = LOG_ZERO;
			}	
			else {
				////Fm[i][j] = mm * Fm[i-1][j-1] + mxy * Fx[i-1][j-1] + mxy * Fy[i-1][j-1];
				////Fm[i][j] = Fm[i][j] * log_odds(i,j);
				Fm[i][j] = LOG_ADD(mm+Fm[i-1][j-1],mxy+Fx[i-1][j-1], mxy+Fy[i-1][j-1])+log_odds(i,j);
				// waive gap open penalty for ends
				////if(j==leny)Fx[i][j] = xy * Fm[i-1][j] + xy * Fx[i-1][j];
				////else Fx[i][j] = xym * Fm[i-1][j] + xy * Fx[i-1][j];
				//Fx[i][j] = xym * Fm[i-1][j] + xy * Fx[i-1][j];
				////Fx[i][j] = Fx[i][j] * log_odds_x(i);

				if(j==leny)Fx[i][j] = LOG_ADD(xy+Fm[i-1][j], xy+Fx[i-1][j]);
				else Fx[i][j] = LOG_ADD(xym+Fm[i-1][j], xy+Fx[i-1][j]);
				Fx[i][j] = Fx[i][j] + log_odds_x(i);

				// waive gap open penalty for ends
				////if(i==lenx)Fy[i][j] = xy * Fm[i][j-1] + xy * Fy[i][j-1];
				////else Fy[i][j] = xym * Fm[i][j-1] + xy * Fy[i][j-1];
				//Fy[i][j] = xym * Fm[i][j-1] + xy * Fy[i][j-1];
				////Fy[i][j] = Fy[i][j] * log_odds_y(j);
				
				if(i==lenx)Fy[i][j] = LOG_ADD(xy+Fm[i][j-1],xy+Fy[i][j-1]);
				else Fy[i][j] = LOG_ADD(xym+Fm[i][j-1],xy+Fy[i][j-1]);
				Fy[i][j] = Fy[i][j] + log_odds_y(j);
			}
			////if(debug>1) fprintf(stdout, "%d %d: %d %f   %d %f   %d %f\n", i,j,Fm[i][j].INT,Fm[i][j].DEC,Fx[i][j].INT,Fx[i][j].DEC,Fy[i][j].INT,Fy[i][j].DEC);
			if(debug>1) fprintf(stdout, "%d %d: %f   %f   %f\n", i,j,Fm[i][j],Fx[i][j],Fy[i][j]);
			fflush(stdout);
		}
	}

	////FE = E * (Fm[lenx][leny] + Fx[lenx][leny] + Fy[lenx][leny]);
	FE = E + LOG_ADD(Fm[lenx][leny],Fx[lenx][leny],Fy[lenx][leny]);
	fprintf(stdout, "FE: %f\n", FE);

}

// this is correponds to sum-of-pairs of the log-odds scores of weighted amino acid pairs
ScoreType hmm_profpair::log_odds(int xi, int yj) {
	
	int i,j;

	////double score = 0;
	float score = 0;

	////if(xi>lenx) return ScoreType(0);
	////if(yj>leny) return ScoreType(0);
	if(xi>lenx) return LOG_ZERO;
	if(yj>leny) return LOG_ZERO;
	
	for(i=1;i<=20;i++) {
	   for(j=1;j<=20;j++) {
		////score += (x->pseudoCnt[xi][i] * y->pseudoCnt[yj][j] * log(q_blosum62[i][j]) );
		score += (x->pseudoCnt[xi][i] * y->pseudoCnt[yj][j] * log_q_blosum62[i][j] );
	   }
	}

	// if(xi==lenx) if(yj==leny) { for(i=1;i<=20;i++) {
	//		fprintf(stdout, "*  %d %f  %f\n", i, x->pseudoCnt[xi][i], y->pseudoCnt[yj][i]); 
	//	}
	//}

	////return ScoreType(score,1);
	return score;
}
	
// this is correponds to log-odds scores of weighted amino acid frequencies
ScoreType hmm_profpair::log_odds_x(int xi) {
	
	int i,j;

	////double score = 0;
	float score = 0;
	
	////if(xi>lenx) return ScoreType(0);
	if(xi>lenx) return LOG_ZERO;
	
	for(i=1;i<=20;i++) {
		////score += (x->pseudoCnt[xi][i] * log(robinson_freq[i]) );
		score += (x->pseudoCnt[xi][i] * log_robinson_freq[i] );
	}

	// if(xi==lenx) if(yj==leny) { for(i=1;i<=20;i++) {
	//		fprintf(stdout, "*  %d %f  %f\n", i, x->pseudoCnt[xi][i], y->pseudoCnt[yj][i]); 
	//	}
	//}

	////return ScoreType(score,1);
	return score;
}
	
// this is correponds to log-odds scores of weighted amino acid frequencies
ScoreType hmm_profpair::log_odds_y(int yj) {
	
	int i,j;

	////double score = 0;
	float score = 0;
	
	////if(yj>leny) return ScoreType(0);
	if(yj>leny) return LOG_ZERO;
	
	for(i=1;i<=20;i++) {
		////score += (y->pseudoCnt[yj][i] * log(robinson_freq[i]) );
		score += (y->pseudoCnt[yj][i] * log_robinson_freq[i] );
	}

	// if(xi==lenx) if(yj==leny) { for(i=1;i<=20;i++) {
	//		fprintf(stdout, "*  %d %f  %f\n", i, x->pseudoCnt[xi][i], y->pseudoCnt[yj][i]); 
	//	}
	//}

	////return ScoreType(score,1);
	return score;
}
	
// this is correponds to log-odds scores of weighted amino acid frequencies
float hmm_profpair::log_odds_x_d(int xi) {
	
	int i,j;

	////double score = 0;
	float score = 0;
	
	//if(xi>lenx) return ScoreType(0);
	
	for(i=1;i<=20;i++) {
		////score += (x->pseudoCnt[xi][i] * log(robinson_freq[i]) );
		score += (x->pseudoCnt[xi][i] * log_robinson_freq[i]);
	}

	// if(xi==lenx) if(yj==leny) { for(i=1;i<=20;i++) {
	//		fprintf(stdout, "*  %d %f  %f\n", i, x->pseudoCnt[xi][i], y->pseudoCnt[yj][i]); 
	//	}
	//}

	//return ScoreType(score,1);
	return score;
}
	
// this is correponds to log-odds scores of weighted amino acid frequencies
float hmm_profpair::log_odds_y_d(int yj) {
	
	int i,j;

	float score = 0;
	
	//if(yj>leny) return ScoreType(0);
	
	for(i=1;i<=20;i++) {
		////score += (y->pseudoCnt[yj][i] * log(robinson_freq[i]) );
		score += (y->pseudoCnt[yj][i] * log_robinson_freq[i] );
	}

	// if(xi==lenx) if(yj==leny) { for(i=1;i<=20;i++) {
	//		fprintf(stdout, "*  %d %f  %f\n", i, x->pseudoCnt[xi][i], y->pseudoCnt[yj][i]); 
	//	}
	//}

	//return ScoreType(score,1);
	return score;
}
	
void hmm_profpair::backward() {

	int i,j,k;

	ScoreType mm, mxy, xym, xy, E;

	mm = LOG( 1-2*delta-tau );
	mxy = LOG( 1-epsilon-tau );
	xym = LOG( delta );
	xy = LOG( epsilon );
	E = LOG( tau );

	Bm = new ScoreType * [lenx+2];
	Bx = new ScoreType * [lenx+2];
	By = new ScoreType * [lenx+2];

	for(i=0;i<=lenx+1;i++) {
		Bm[i] = new ScoreType [leny+2];
		Bx[i] = new ScoreType [leny+2];
		By[i] = new ScoreType [leny+2];
	}
	
	if(debug>1) fprintf(stdout, "===================\n");	
	Bm[lenx][leny] = Bx[lenx][leny] = By[lenx][leny] = E;
	for(i=0;i<=lenx+1;i++) Bm[i][leny+1] = Bx[i][leny+1] = By[i][leny+1] = LOG_ZERO;
	for(i=0;i<=leny+1;i++) Bm[lenx+1][i] = Bx[lenx+1][i] = By[lenx+1][i] = LOG_ZERO;
	if(debug>1) fprintf(stdout, "===================\n");	
	for(i=lenx;i>=0;i--) {
		for(j=leny;j>=0;j--) {
			if( (i==lenx) && (j==leny) ) continue;
			if( i== lenx) {
				////Bm[i][j] = xy * log_odds_y(j+1) * By[i][j+1];
				Bm[i][j] = xy + log_odds_y(j+1) + By[i][j+1];
			}
			else if(j==leny){
				////Bm[i][j] = xy *log_odds_x(i+1) * Bx[i+1][j];
				Bm[i][j] = xy + log_odds_x(i+1) + Bx[i+1][j];
			}
			else 
			////Bm[i][j] = mm * log_odds(i+1,j+1) * Bm[i+1][j+1] + xym * (log_odds_x(i+1) * Bx[i+1][j] + log_odds_y(j+1) * By[i][j+1]);
			Bm[i][j] = LOG_ADD(mm + log_odds(i+1,j+1) + Bm[i+1][j+1], xym + LOG_ADD(log_odds_x(i+1) + Bx[i+1][j] , log_odds_y(j+1) + By[i][j+1]));

			if( (j==0)&&(i==0) ){
				////Bx[i][j] = xy * log_odds(i+1,j+1) * Bm[i+1][j+1] + xy * log_odds_x(i+1) * Bx[i+1][j];
				Bx[i][j] = LOG_ADD( xy + log_odds(i+1,j+1) + Bm[i+1][j+1] , xy + log_odds_x(i+1) + Bx[i+1][j]);
			}
			else 
			////Bx[i][j] = mxy * log_odds(i+1,j+1) * Bm[i+1][j+1] + xy * log_odds_x(i+1) * Bx[i+1][j];
			Bx[i][j] = LOG_ADD( mxy + log_odds(i+1,j+1) + Bm[i+1][j+1] , xy + log_odds_x(i+1) + Bx[i+1][j]);
			if( (j==0) &&(i==0) ){
				////By[i][j] = xy * log_odds(i+1,j+1) * Bm[i+1][j+1] + xy * log_odds_y(j+1) * By[i][j+1];	
				By[i][j] = LOG_ADD( xy + log_odds(i+1,j+1) + Bm[i+1][j+1] , xy + log_odds_y(j+1) + By[i][j+1]);	
			}
			else 
			////By[i][j] = mxy * log_odds(i+1,j+1) * Bm[i+1][j+1] + xy * log_odds_y(j+1) * By[i][j+1];	
			By[i][j] = LOG_ADD(mxy + log_odds(i+1,j+1) + Bm[i+1][j+1] , xy + log_odds_y(j+1) + By[i][j+1]);	
			////if(debug>1) fprintf(stdout, "%d %d: %d %f   %d %f   %d %f\n", i,j,Bm[i][j].INT,Bm[i][j].DEC,Bx[i][j].INT,Bx[i][j].DEC,By[i][j].INT,By[i][j].DEC);
			if(debug>1) fprintf(stdout, "%d %d: %f   %f   %f\n", i,j,Bm[i][j],Bx[i][j],By[i][j]);
		}
	}

	////BE = mm * Bm[1][1] * log_odds(1,1) + xy * Bx[1][0] * log_odds_x(1) + xy * By[0][1] * log_odds_y(1);
	////fprintf(stdout, "BE: %d %f Bm 0 0: %d %f\n", BE.INT, BE.DEC, Bm[0][0].INT, Bm[0][0].DEC);
	BE = LOG_ADD(mm + Bm[1][1] + log_odds(1,1) , xy + Bx[1][0] + log_odds_x(1) , xy + By[0][1] + log_odds_y(1) );
	fprintf(stdout, "BE: %f Bm 0 0: %f\n", BE, Bm[0][0]);

	//fprintf(stdout, "================\n");

	// not the full probabities
	////ScoreType full_prob;
	////for(i=1;i<=lenx;i++) { for(j=1;j<=leny;j++) { full_prob = Fm[i][j] * Bm[i][j] + Fx[i][j] * Bx[i][j] + Fy[i][j] * By[i][j]; //fprintf(stdout, "%d %d: %d %f\n", i, j, full_prob.INT, full_prob.DEC); } }

	//fprintf(stdout, "================\n");

	//ScoreType full_prob;
	cout << "Obtain posterior probabilities " << endl;
	ScoreType **aligned_pair_prob = gmatrix<ScoreType>(lenx, leny);
	probMat = gmatrix<float>(lenx, leny);
	for(i=1;i<=lenx;i++) {
		for(j=1;j<=leny;j++) {
			////aligned_pair_prob[i][j] = Fm[i][j] * Bm[i][j] / BE;
			aligned_pair_prob[i][j] = Fm[i][j] + Bm[i][j] - BE;
			////probMat[i][j] = aligned_pair_prob[i][j].real();
			////probMat[i][j] = EXP(aligned_pair_prob[i][j]);
			////why EXP gives negative values?
			probMat[i][j] = exp(aligned_pair_prob[i][j]);
			////if(debug>2) fprintf(stdout, "%d %d: %d %f %f\n", i, j, aligned_pair_prob[i][j].INT, aligned_pair_prob[i][j].DEC, aligned_pair_prob[i][j].real());
			if(debug>2) fprintf(stdout, "%d %d: %f %f\n", i, j, aligned_pair_prob[i][j], probMat[i][j]);
		}
	}

	free_gmatrix<ScoreType>(aligned_pair_prob, lenx, leny);

	//fprintf(stdout, "================\n");

}


	
void hmm_profpair::getPosterior() {
	
	; // nothing here now

}

	
int * hmm_profpair::posteriorAlignment() {

	int i,j,k;
	int *path;

	int **scoreMat = imatrix(lenx, leny);

	for(i=1;i<=lenx;i++) {
		for(j=1;j<=leny;j++) {
			scoreMat[i][j] = (int) (probMat[i][j]*1000);
		}
	}

	MM galign;
	galign.setM(lenx);
	galign.setN(leny);
	galign.set_g(0); // gap open penalty
	galign.set_h(0); // gap extension penalty

	cout << "Compute pairwise consistency alignment:" << endl;
	galign.dp(scoreMat);
	//cout << "Now here" << endl;
	path = new int [galign.print_ptr];
	for(i=1;i<galign.print_ptr;i++) {
		path[i] = galign.displ[i];
		//cout << "i: " << i << " " << path[i] << endl;
	}
	//len = galign.print_ptr - 1;
	//cout << "len: " << len << endl;
	//cout << "Pairwise consistency alignment ends here" << endl;	

	free_imatrix(scoreMat, lenx, leny);

	char *a, *b;
	a = new char[1000];
	b = new char[1000];

	int i1,i2;
	i1 = i2 = 0;
	int x1 = 0;
	for(i=1;i<galign.print_ptr;i++) {
		if(path[i]<0) {
			for(k=1;k<=0-path[i];k++) {
				a[x1] = x->aseq[0][i1];
				b[x1] = '-';
				i1++;
				x1++;
			}
			continue;
		}
		if(path[i]>0) {
			for(k=1;k<=path[i];k++) {
				a[x1] = '-';
				b[x1] = y->aseq[0][i2];
				i2++;
				x1++;
			}
			continue;
		}
		if(path[i]==0) {
			a[x1] = x->aseq[0][i1];
			b[x1] = y->aseq[0][i2];
			i1++;
			i2++;
			x1++;
		}
	}
	a[x1] = '\0';
	b[x1] = '\0';
	
	cout << left << setw(20) << x->aname[0] << a << endl;
	cout << left << setw(20) << y->aname[0] << b << endl;

	return path;
}

