#include "all.h"
#include "hmm_local.h"

extern const float LOG_ZERO;
static int debug = 0;
hmm_local::hmm_local(subalign *x1, subalign *y1) {

	x = x1;
	y = y1;
	lenx = x->alilen;
	leny = y->alilen;

	gapopen_end = 0;

	delta = 0.019931; 
	epsilon = 0.79433; 
	tau = 0.19598;
	//tau = 0;
	eta = 0.3;
	theta = 0.19598;
/*
New eta: 0.207316
New delta: 0.0241507
New epsilon: 0.863795
New tau: 0.00897139
*/

	eta = 0.207316;
	delta = 0.0241507;
	epsilon = 0.863795;
	tau = 0.00897139;

/*
New eta: 0.251283
New delta: 0.049276
New epsilon: 0.769577
New tau: 0.00985969
*/

	eta = 0.251283;
	delta = 0.049276;
	epsilon = 0.769577;
	tau = 0.00985969;
/*
NEW eta: 0.243145
NEW delta: 0.0490766
NEW epsilon: 0.769577
NEW tau: 0.00985969
*/

/*
NEW eta: 0.246077
NEW delta: 0.0668055
NEW epsilon: 0.708576
NEW tau: 0.00965374
*/
	eta = 0.246077;
	delta = 0.0668055;
	epsilon = 0.708576;
	tau = 0.00965374;

	Vm=Vx=Vy=0; Tm=Tx=Ty=0;
	path=path1=path2=0;
	Fm=Fx=Fy=Bm=Bx=By=0;
	Frx1 = Fn1 = 0;
	Fry1 = Fn2 = Fn3 = Frx2 = Fn4 = Fry2 = Fe = 0;
	Brx1 = Bn1 = Bry1 = Bn2 = Bn3 = Brx2 = 0;
	Bn4 = Bry2 = 0;
	probMat = 0;
}

void hmm_local::deleteSubalign() {

	delete x;
	delete y;
}
	
hmm_local::~hmm_local() {

	free_gmatrix<float>(Vm, lenx, leny);
	free_gmatrix<float>(Vx, lenx, leny);
	free_gmatrix<float>(Vy, lenx, leny);
	free_gmatrix<int>(Tm, lenx, leny);
	free_gmatrix<int>(Tx, lenx, leny);
	free_gmatrix<int>(Ty, lenx, leny);
	free_gmatrix<float>(Fm, lenx, leny);
	free_gmatrix<float>(Fx, lenx, leny);
	free_gmatrix<float>(Fy, lenx, leny);
	free_gmatrix<float>(Bm, lenx+1, leny);
	free_gmatrix<float>(Bx, lenx+1, leny);
	free_gmatrix<float>(By, lenx+1, leny);

	free_gvector<float>(Frx1);
	free_gvector<float>(Fn1);
	free_gmatrix<float>(Fry1, lenx,leny);
	free_gmatrix<float>(Fn2, lenx,leny);
	free_gmatrix<float>(Fn3, lenx,leny);
	free_gmatrix<float>(Frx2, lenx,leny);
	free_gmatrix<float>(Fry2, lenx,leny);
	free_gmatrix<float>(Fn4, lenx,leny);
	free_gmatrix<float>(Fe, lenx,leny);
	free_gvector<float>(Bry2);
	free_gvector<float>(Bn4);
	free_gmatrix<float>(Bry1, lenx+1,leny);
	free_gmatrix<float>(Bn2, lenx+1,leny);
	free_gmatrix<float>(Bn3, lenx+1,leny);
	free_gmatrix<float>(Brx2, lenx+1,leny);
	free_gmatrix<float>(Brx1, lenx+1,leny);
	free_gmatrix<float>(Bn1, lenx+1,leny);
	
	

	free_gmatrix<float>(probMat, lenx, leny);

	free_gvector<int>(path);
	free_gvector<int>(path1);
	free_gvector<int>(path2);
	

}	

void hmm_local::viterbi() {

	int i,j,k;
	int xi, yj, a;

	float mm, mxy, xym, xy, E;
	int maxalnlength = lenx + leny +1;
	int *reverse_path, *reverse_path1, *reverse_path2;

	mm = log(1-2*delta-tau);
	mxy = log(1-epsilon-tau);
	xym = log(delta);
	xy = log(epsilon);
	//E = log(tau);
	E = 0;

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
	//Vm[0][0] = 0; //Vx[0][0] = Vy[0][0] = 0;
	Vm[0][0] = log(theta) - xy; //Vx[0][0] = Vy[0][0] = 0;
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
	
	Vm[0][0] = log(1-2*theta);
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

	VE = E + max3(Vm[lenx][leny]+LOG(1-2*theta), Vx[lenx][leny]+LOG(theta), Vy[lenx][leny]+LOG(theta), TE);

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
	fprintf(stdout, "=================\n");
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
	fprintf(stdout, "=================\n");

	v_alilen = i-1;
	fprintf(stdout, "v_alilen: %d\n", v_alilen);
	for(i=1;i<=v_alilen;i++) {
		path[i] = reverse_path[v_alilen-i+1];
		path1[i] = reverse_path1[v_alilen-i+1];
		path2[i] = reverse_path2[v_alilen-i+1];
		if(debug>1) fprintf(stdout, "%d %d %d\n", i, reverse_path1[i], reverse_path2[i]);
	}
	if(debug>1) fprintf(stdout, "\n");
	delete [] reverse_path;
	delete [] reverse_path1;
	delete [] reverse_path2;
	
}
		
// this is correponds to sum-of-pairs of the log-odds RATIO scores of weighted amino acid pairs
float hmm_local::log_odds_score(int xi, int yj) {
	
	int i,j;

	float score = 0;
	
	for(i=1;i<=20;i++) {
	   if(!x->pseudoCnt[xi][i]) continue;
	   for(j=1;j<=20;j++) {
		////score += (x->pseudoCnt[xi][i] * y->pseudoCnt[yj][j] * log(q_blosum62[i][j]/robinson_freq[i]/robinson_freq[j]) );
		if(!y->pseudoCnt[yj][j]) continue;
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
	
float hmm_local::max3(float a,float b, float c, int &t) {

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

	
float hmm_local::max2(float a, float b, int &t) {
	
	if(a>=b) {
		t= 1;
		return a;
	}
	else {
		t = 2;
		return b;
	}
}

void hmm_local::printHmmAlign(int blocksize) {

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
		fprintf(stdout, "=================\n\n");
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
subalign *hmm_local::productViterbiAlign() { 

	int i,j,k;

	//subalign *newAlign = new subalign();
	subalign *a = new subalign(); // newAlign;

	a->nal = x->nal + y->nal;
	a->alilen = v_alilen;
	if(x->mnamelen > y->mnamelen) a->mnamelen = x->mnamelen;
	else a->mnamelen = y->mnamelen;
	if(debug>-1) cout << "a->mnamelen: " << a->mnamelen <<endl;
	if(debug>-1) cout << "a->alilen: " << a->alilen <<endl;
	if(debug>-1) cout << "a->nal: " << a->nal <<endl;
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

	if(debug>-1) for(i=1;i<=v_alilen;i++) { cout << i << "\t" << path[i] << endl; }

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

void hmm_local::getTransitions() {

	e = L_eta = log(eta);
	e1 = L1_eta = log(1-eta);
	d = L_delta = log(delta);
	d2t1 = L1_delta2_tau = log(1-2*delta-tau);
	t = L_tau = log(tau);
	ep = L_epsilon = log(epsilon);
	ept1 = L1_epsilon_tau = log(1-epsilon-tau);

}

void hmm_local::forward() {

	int i,j,k;

	ScoreType mm, mxy, xym, xy, E;

	mm = LOG( 1-2*delta-tau );
	mxy = LOG( 1-epsilon-tau );
	xym = LOG( delta );
	xy = LOG( epsilon );
	//E = LOG( tau );
	E = LOG_ONE;

	getTransitions();

	cout << "e:"<<e << "\te1:" <<  e1 << "\td:" <<  d << "\td2t1:" <<  d2t1 << "\tt:" <<  t << "\tep:" <<  ep << "\tept1:" <<  ept1 << endl;

	Fm = new ScoreType * [lenx+1];
	Fx = new ScoreType * [lenx+1];
	Fy = new ScoreType * [lenx+1];
	Frx1 = new ScoreType [lenx+1];
	Frx2 = new ScoreType * [lenx+1];
	Fry1 = new ScoreType * [lenx+1];
	Fry2 = new ScoreType * [lenx+1];
	Fn1 = new ScoreType [lenx+1];
	Fn2 = new ScoreType * [lenx+1];
	Fn3 = new ScoreType * [lenx+1];
	Fn4 = new ScoreType * [lenx+1];
	Fe = new ScoreType * [lenx+1];

	for(i=0;i<=lenx;i++) {
		Fm[i] = new ScoreType [leny+1];
		Fx[i] = new ScoreType [leny+1];
		Fy[i] = new ScoreType [leny+1];
		Frx2[i] = new ScoreType [leny+1];
		Fry1[i] = new ScoreType [leny+1];
		Fry2[i] = new ScoreType [leny+1];
		Fn2[i] = new ScoreType [leny+1];
		Fn3[i] = new ScoreType [leny+1];
		Fn4[i] = new ScoreType [leny+1];
		Fe[i] = new ScoreType [leny+1];
	}

	//cout << "Here: " << endl;

	// Frx1 and Fn1
	//cout << "Here: Frx1 and Fn1" << endl;
	Frx1[0] = LOG_ONE;
	Fn1[0] = e;
	for(i=1;i<=lenx;i++) {
		Frx1[i] = Frx1[i-1] + e1 + log_odds_x(i);
		Fn1[i] = Frx1[i] + e;
		//cout << i << ": " << Frx1[i] << "\t" << Fn1[i] << endl;
	}

	// Fry1 and Fn2
	//cout << "Here: Fry1 and Fn2" << endl;
	for(i=0;i<=lenx;i++) {
		j=0;
		Fry1[i][0] = Fn1[i];
		Fn2[i][0] = Fn1[i] + e;
		//cout << i << "\t" << j << ": " << Fry1[i][j] << "\t" << Fn2[i][j] << endl;
		for(j=1;j<=leny;j++) {
			Fry1[i][j] = Fry1[i][j-1] + e1 + log_odds_y(j);	
			Fn2[i][j] = Fry1[i][j] + e;
			//cout << i << "\t" << j << ": " << Fry1[i][j] << "\t" << Fn2[i][j] << endl;
		}
	}

	// Fm, Fx and Fy
	// initial conditions for Fm, Fx, Fy
	//cout << "Here: Fm Fx and Fy" << endl;
	cout << "LOG_ZERO: " << LOG_ZERO << endl;
	for(i=0;i<=lenx;i++) {
		Fm[i][0] = Fy[i][0] = LOG_ZERO;
		if(i==0) Fx[i][0] = LOG_ZERO;
		else Fx[i][0] = LOG_ADD(Fx[i-1][0]+ep, Fn2[i-1][0]+d)+log_odds_x(i);
	}
	for(j=0;j<=leny;j++) {
		Fm[0][j] = Fx[0][j] = LOG_ZERO;
		if(j==0) Fy[0][j] = LOG_ZERO;
		else Fy[0][j] = LOG_ADD(Fy[0][j-1]+ep, Fn2[0][j-1]+d)+log_odds_y(j);
	}
	for(i=1;i<=lenx;i++) {
		for(j=1;j<=leny;j++) {
			Fm[i][j] = LOG_ADD(Fm[i-1][j-1]+d2t1,Fx[i-1][j-1]+ept1,Fy[i-1][j-1]+ept1,Fn2[i-1][j-1]+d2t1)+log_odds(i,j);
			Fx[i][j] = LOG_ADD(Fm[i-1][j]+d, Fx[i-1][j]+ep, Fn2[i-1][j]+d)+log_odds_x(i);
			Fy[i][j] = LOG_ADD(Fm[i][j-1]+d, Fy[i][j-1]+ep, Fn2[i][j-1]+d)+log_odds_y(j);
		}
	}
	//for(i=1;i<=lenx;i++) { cout << log_odds_x(i) << endl; } cout << endl;
	//for(i=1;i<=leny;i++) { cout << log_odds_y(i) << endl; }

	//cout << "Fx 1 0: " << e+e+d+log_odds_x(1)<<endl;
	//cout << "Fy 1 0: " << e+e+d+log_odds_y(1)<<endl;
	//cout << "Fx 1 1: " << e+e1+log_odds_x(1)+e+d+log_odds_y(1) << endl;
	//cout << "Fm 1 1: " << e+e+d2t1+log_odds(1,1) << endl;
	//for(i=0;i<=lenx;i++) { for(j=0;j<=leny;j++) { cout << i << "\t" << j << ": " << Fm[i][j] << "\t" << Fx[i][j] << "\t" << Fy[i][j] << "\t" << log_odds(i,j) << endl; } }

	// Fn3
	//cout << "Here: Fn3" << endl;
	//cout << "Fn3 1 0: " << LOG_ADD(e+e+d+log_odds_x(1)+t, e1+log_odds_x(1)+e+e+t)<<endl;
	for(i=0;i<=lenx;i++) {
		for(j=0;j<=leny;j++) {
			if( (i==0) && (j==0) ) { Fn3[i][j] = Fn2[0][0] + t; continue; }
			Fn3[i][j] = LOG_ADD(Fm[i][j]+t, Fx[i][j]+t, Fy[i][j]+t, Fn2[i][j]+t);
			//cout << i << "\t" << j << ": " << Fn3[i][j] << endl;
		}
	}

	// Frx2 and Fn4
	//cout << "Here: Frx2 and Fn4" << endl;
	//cout << "Frx2 1 0: " << e+e+t+e1+log_odds_x(1)<<endl;
	//cout << "Fn4 1 0: " << e+e+t+e1+log_odds_x(1)+e<<endl;
	//cout << "Frx2 1 1: " << e+e1+log_odds_x(1)+e+d+log_odds_y(1) << endl;
	//cout << "Fn4 1 1: " << e+e+d2t1+log_odds(1,1) << endl;
	for(j=0;j<=leny;j++) {
		for(i=0;i<=lenx;i++) {
			if(i==0) Frx2[i][j] = LOG_ZERO;
			else Frx2[i][j] = LOG_ADD(Frx2[i-1][j]+e1, Fn3[i-1][j]+e1) + log_odds_x(i);;
			Fn4[i][j] = LOG_ADD(Frx2[i][j], Fn3[i][j]) + e;
			//cout << i << "\t" << j << ": " << Frx2[i][j] << "\t" << Fn4[i][j] << endl;
		}
	}

	// Fry2 an Fe
	//cout << "Here: Fry2 and Fe" << endl;
	for(i=0;i<=lenx;i++) {
		for(j=0;j<=leny;j++) {
			if(j==0) Fry2[i][j] = LOG_ZERO;
			else Fry2[i][j] = LOG_ADD(Fn4[i][j-1]+e1, Fry2[i][j-1]+e1) + log_odds_y(j);
			Fe[i][j] = LOG_ADD(Fn4[i][j]+e, Fry2[i][j]+e);
			//cout << i << "\t" << j << ": " << Fry2[i][j] << "\t" << Fe[i][j] << endl;
		}
	}
	//cout << endl;
		
	
	fprintf(stdout, "FE: Fry2: %f %f\n", Fe[lenx][leny], Fry2[lenx][leny]);
	//cout << "FE: " << LOG_ADD(2*e1+log_odds_x(1)+log_odds_y(1)+log(4.0), e1+d+log_odds_x(1)+log_odds_y(1)+log(4.0), d2t1+log_odds(1,1) ) + 4*e + t << endl;

}

// this is correponds to sum-of-pairs of the log-odds scores of weighted amino acid pairs
ScoreType hmm_local::log_odds(int xi, int yj) {
	
	int i,j;

	////double score = 0;
	float score = 0;

	////if(xi>lenx) return ScoreType(0);
	////if(yj>leny) return ScoreType(0);
	if(xi>lenx) return LOG_ZERO;
	if(yj>leny) return LOG_ZERO;
	if( (xi<=0) || ( yj<=0) ) return LOG_ZERO;
	
	for(i=1;i<=20;i++) {
	   if(!x->pseudoCnt[xi][i]) continue;
	   for(j=1;j<=20;j++) {
		////score += (x->pseudoCnt[xi][i] * y->pseudoCnt[yj][j] * log(q_blosum62[i][j]) );
		if(!y->pseudoCnt[yj][j]) continue;
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
ScoreType hmm_local::log_odds_x(int xi) {
	
	int i,j;

	////double score = 0;
	float score = 0;
	
	////if(xi>lenx) return ScoreType(0);
	if(xi>lenx) return LOG_ZERO;
	if(xi<=0) return LOG_ZERO;
	
	for(i=1;i<=20;i++) {
	   	if(!x->pseudoCnt[xi][i]) continue;
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
ScoreType hmm_local::log_odds_y(int yj) {
	
	int i,j;

	////double score = 0;
	float score = 0;
	
	////if(yj>leny) return ScoreType(0);
	if(yj>leny) return LOG_ZERO;
	if(yj<=0) return LOG_ZERO;
	
	for(i=1;i<=20;i++) {
		if(!y->pseudoCnt[yj][i]) continue;
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
float hmm_local::log_odds_x_d(int xi) {
	
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
float hmm_local::log_odds_y_d(int yj) {
	
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
	
void hmm_local::backward() {

	int i,j,k;

	ScoreType mm, mxy, xym, xy, E;

	mm = LOG( 1-2*delta-tau );
	mxy = LOG( 1-epsilon-tau );
	xym = LOG( delta );
	xy = LOG( epsilon );
	//E = LOG( tau );
	E = LOG_ONE;

	getTransitions();

	//cout << e << "\t" <<  e1 << "\t" <<  d << "\t" <<  d2t1 << "\t" <<  t << "\t" <<  ep << "\t" <<  ept1 << endl;

	Bm = new ScoreType * [lenx+2];
	Bx = new ScoreType * [lenx+2];
	By = new ScoreType * [lenx+2];
	Brx1 = new ScoreType * [lenx+2];
	Brx2 = new ScoreType * [lenx+2];
	Bry1 = new ScoreType * [lenx+2];
	Bry2 = new ScoreType [leny+2];
	Bn1 = new ScoreType * [lenx+2];
	Bn2 = new ScoreType * [lenx+2];
	Bn3 = new ScoreType * [lenx+2];
	Bn4 = new ScoreType  [leny+2];


	for(i=0;i<=lenx+1;i++) {
		Bm[i] = new ScoreType [leny+2];
		Bx[i] = new ScoreType [leny+2];
		By[i] = new ScoreType [leny+2];
		Bn1[i] = new ScoreType [leny+2];
		Brx1[i] = new ScoreType [leny+2];
		Brx2[i] = new ScoreType [leny+2];
		Bry1[i] = new ScoreType [leny+2];
		Bn2[i] = new ScoreType [leny+2];
		Bn3[i] = new ScoreType [leny+2];
	}

	cout << "LOG_ZERO: " << LOG_ZERO << endl;

	//cout << "Be Here: " << endl;

	// Bry2
	Bry2[leny] = e;
	for(j=leny-1;j>=0;j--) {
		Bry2[j] = Bry2[j+1]+log_odds_y(j+1)+e1;
	}
	//cout << "Be Here: " << endl;

	// Bn4
	Bn4[leny] = e;
	for(j=leny-1;j>=0;j--) {
		Bn4[j] = Bry2[j+1]+log_odds_y(j+1)+e1;
	}
	//cout << "Be Here: " << endl;

	// Brx2
	for(j=leny;j>=0;j--) {
		Brx2[lenx][j] = Bn4[j] + e;
	}
	for(i=lenx-1;i>=0;i--) {
	    for(j=leny;j>=0;j--) {
		Brx2[i][j] = Brx2[i+1][j] + log_odds_x(i+1) + e1;
	    }
	}
	//cout << "Be Here: " << endl;

	// Bn3
	for(j=leny;j>=0;j--) Bn3[lenx][j] = Bn4[j]+e;
	for(i=lenx-1;i>=0;i--) {
            for(j=leny;j>=0;j--) {
		Bn3[i][j] = Brx2[i+1][j] + log_odds_x(i+1) + e1;
	    }
	}
	//cout << "Be Here: " << endl;

	// Bm, Bx, By
	for(i=lenx+1;i>=0;i--) Bm[i][leny+1] = Bx[i][leny+1] = By[i][leny+1] = LOG_ZERO;
	for(i=leny+1;i>=0;i--) Bm[lenx+1][i] = Bx[lenx+1][i] = By[lenx+1][i] = LOG_ZERO;
	for(i=lenx;i>=0;i--) {
	    for(j=leny;j>=0;j--) {
		Bm[i][j] = LOG_ADD(Bm[i+1][j+1]+log_odds(i+1,j+1)+d2t1, Bx[i+1][j]+log_odds_x(i+1)+d, By[i][j+1]+log_odds_y(j+1)+d, Bn3[i][j]+t);
		Bx[i][j] = LOG_ADD(Bm[i+1][j+1]+log_odds(i+1,j+1)+ept1, Bx[i+1][j]+log_odds_x(i+1)+ep, Bn3[i][j]+t);
		By[i][j] = LOG_ADD(Bm[i+1][j+1]+log_odds(i+1,j+1)+ept1, By[i][j+1]+log_odds_y(j+1)+ep, Bn3[i][j]+t);
	   }
	}
	//cout << "Be Here: " << endl;
	
	// Bn2
	for(i=lenx;i>=0;i--) {
	    for(j=leny;j>=0;j--) {
		Bn2[i][j] = LOG_ADD(Bm[i+1][j+1]+d2t1+log_odds(i+1,j+1), Bx[i+1][j]+d+log_odds_x(i+1), By[i][j+1]+d+log_odds_y(j+1), Bn3[i][j]+t);
	    }
	}
	cout << "Be Here: " << endl;
	
	// Bry1
	for(i=lenx;i>=0;i--) {
		Bry1[i][leny+1]=LOG_ZERO;
	}
	for(j=leny;j>=0;j--) {
	    for(i=lenx;i>=0;i--) {
		Bry1[i][j] = LOG_ADD(Bry1[i][j+1]+log_odds_y(j+1)+e1, Bn2[i][j]+e);
	    }
	}
	cout << "Be Here: " << endl;

	// Bn1
	for(j=leny;j>=0;j--) {
            for(i=lenx;i>=0;i--) {
		Bn1[i][j] = LOG_ADD(Bry1[i][j+1]+e1+log_odds_y(j+1), Bn2[i][j]+e);
	    }
	}
	cout << "Be Here: " << endl;
	
	// Brx1
	for(j=leny;j>=0;j--) {
	    Brx1[lenx+1][j] = LOG_ZERO;
            for(i=lenx;i>=0;i--) {
		Brx1[i][j] = LOG_ADD(Brx1[i+1][j]+log_odds_x(i+1)+e1, Bn1[i][j]+e);
	    }
	}
	//cout << "Be Here: " << endl;

	BE = LOG_ADD(Brx1[1][0]+e1+log_odds_x(1), Bn1[0][0]+e);
	//Be[0][0] = LOG_ADD(Brx1[0][0]+e1, Bn1[0][0]+e);

	//for(i=0;i<=lenx;i++) { for(j=0;j<=leny;j++) { cout << i << "\t" << j << ": " << LOG_ADD(Fry1[i][j]+Bry1[i][j], Fn2[i][j]+Bn2[i][j], Fm[i][j]+Bm[i][j], Fx[i][j]+Bx[i][j], Fy[i][j]+By[i][j], Fn3[i][j]+Bn3[i][j], Frx2[i][j]+Brx2[i][j]) << endl; } }

	fprintf(stdout, "BE: %f\n", BE);

	// adjust BE:
	BE = (Fe[lenx][leny] + BE)/2;
	fprintf(stdout, "BE: %f\n", BE);

	fprintf(stdout, "================\n");

	// not the full probabities
	////ScoreType full_prob;
	////for(i=1;i<=lenx;i++) { for(j=1;j<=leny;j++) { full_prob = Fm[i][j] * Bm[i][j] + Fx[i][j] * Bx[i][j] + Fy[i][j] * By[i][j]; //fprintf(stdout, "%d %d: %d %f\n", i, j, full_prob.INT, full_prob.DEC); } }

	fprintf(stdout, "================\n");

	//ScoreType full_prob;
	cout << "Obtain posterior probabilities " << endl;
	//ScoreType **aligned_pair_prob = gmatrix<ScoreType>(lenx, leny);
	probMat = gmatrix<float>(lenx, leny);
	for(i=1;i<=lenx;i++) {
		for(j=1;j<=leny;j++) {
			////aligned_pair_prob[i][j] = Fm[i][j] * Bm[i][j] / BE;
			//aligned_pair_prob[i][j] = Fm[i][j] + Bm[i][j] - BE;
			//cout << aligned_pair_prob[i][j] << endl;
			////probMat[i][j] = aligned_pair_prob[i][j].real();
			////probMat[i][j] = EXP(aligned_pair_prob[i][j]);
			////why EXP gives negative values?
			probMat[i][j] = exp(Fm[i][j] + Bm[i][j] - BE);
			////if(debug>2) fprintf(stdout, "%d %d: %d %f %f\n", i, j, aligned_pair_prob[i][j].INT, aligned_pair_prob[i][j].DEC, aligned_pair_prob[i][j].real());
			if(debug>2) fprintf(stdout, "%d %d: %f\n", i, j, probMat[i][j]);
		}
	}

	//free_gmatrix<ScoreType>(aligned_pair_prob, lenx, leny);

	fprintf(stdout, "================\n");

}

void  hmm_local::baumwelch() {
	
	int i,j;

	getTransitions();

	float rx1rx1, rx1n1, n1ry1, n1n2, ry1ry1, ry1n2, n2m, n2x, n2y, n2n3, mm, mx, my, mn3, xm, xx, xn3, yy, ym, yn3, n3rx2, n3n4, rx2rx2, rx2n4, n4ry2, ry2ry2, n4e, ry2e, LOG_ZERO;
	rx1rx1=rx1n1=n1ry1=n1n2=ry1ry1=ry1n2=n2m=n2x=n2y=n2n3=mm=mx=my=mn3=xm=xx=xn3=yy=ym=yn3=n3rx2=n3n4=rx2rx2=rx2n4=n4ry2=ry2ry2=n4e=ry2e=-2e20;

	for(i=1;i<lenx;i++) {
		if(i<0) {
			rx1rx1 =  Frx1[i]+e1+log_odds_x(i+1)+Brx1[i+1][0];
			rx1n1 = Frx1[i]+e+Bn1[i][0];
		}
		else {
			if(i<lenx) rx1rx1 = LOG_ADD(rx1rx1, Frx1[i]+e1+log_odds_x(i+1)+Brx1[i+1][0]);
			if(i>0)    rx1n1 = LOG_ADD(rx1n1, Frx1[i]+e+Bn1[i][0]);
		}
		//cout << i << ": " << rx1rx1 << "\t" << rx1n1 << "\t" << Frx1[i] << "\t" << log_odds_x(i+1) <<"\t" <<  Brx1[i+1][0] << endl;
	}
	rx1rx1 -= BE;
	rx1n1 -= BE;
	cout << rx1rx1 << endl;
	cout << rx1n1 << endl;

	for(i=1;i<lenx;i++) {
		n1ry1 = LOG_ADD(n1ry1, Fn1[i]+e1+log_odds_y(1)+Bry1[i][1]);
		n1n2 = LOG_ADD(n1n2, Fn1[i]+e+Bn2[i][0]);
	}
	n1ry1 -= BE;
	n1n2 -= BE;
	cout << n1ry1 << endl;
	cout << n1n2 << endl;

	for(i=1;i<lenx;i++) {
		for(j=1;j<leny;j++) {
		    if(j<leny) ry1ry1 = LOG_ADD(ry1ry1, Fry1[i][j]+e1+log_odds_y(j+1)+Bry1[i][j+1]);
		    if(j>0) ry1n2 = LOG_ADD(ry1n2, Fry1[i][j]+e+Bn2[i][j]);
		}
	}
	ry1ry1 -= BE;
	ry1n2 -= BE;
	cout << ry1n2 << endl;
	cout << ry1ry1 << endl;

	for(i=1;i<lenx;i++) {
		for(j=1;j<leny;j++) {
		    if( (i<lenx)&&(j<leny) ) {
			n2m = LOG_ADD(n2m, Fn2[i][j]+d2t1+log_odds(i+1,j+1)+Bm[i+1][j+1]);
			mm = LOG_ADD(mm, Fm[i][j]+d2t1+log_odds(i+1,j+1)+Bm[i+1][j+1]);
			xm = LOG_ADD(xm, Fx[i][j]+ept1+log_odds(i+1,j+1)+Bm[i+1][j+1]);
			ym = LOG_ADD(ym, Fy[i][j]+ept1+log_odds(i+1,j+1)+Bm[i+1][j+1]);
		    }
		    if(i<lenx) {
			n2x = LOG_ADD(n2x, Fn2[i][j]+d+log_odds_x(i+1)+Bx[i+1][j]);
			xx = LOG_ADD(xx, Fx[i][j]+ep+log_odds_x(i+1)+Bx[i+1][j]);
			mx = LOG_ADD(mx, Fm[i][j]+d+log_odds_x(i+1)+Bx[i+1][j]);
		    }
		    if(j<leny) {
			n2y = LOG_ADD(n2y, Fn2[i][j]+d+log_odds_y(j+1)+By[i][j+1]);
			yy = LOG_ADD(yy, Fy[i][j]+ep+log_odds_y(j+1)+By[i][j+1]);
			my = LOG_ADD(my, Fm[i][j]+d+log_odds_y(j+1)+By[i][j+1]);
		    }
		    n2n3 = LOG_ADD(n2n3, Fn2[i][j]+t+Bn3[i][j]);
		    mn3 = LOG_ADD(mn3, Fm[i][j]+t+Bn3[i][j]);
		    xn3 = LOG_ADD(mn3, Fx[i][j]+t+Bn3[i][j]);
		    yn3 = LOG_ADD(mn3, Fy[i][j]+t+Bn3[i][j]);
		}
	}
	n2m -= BE;
	mm-= BE;
	xm -= BE;
	ym -= BE;
	n2x -= BE;
	xx -= BE;
	mx -= BE;
	n2y -= BE;
	yy -= BE;
	my -= BE;
	n2n3 -= BE;
	mn3 -= BE;
	xn3 -= BE;
	yn3 -= BE;

	for(i=1;i<lenx;i++) {
		for(j=1;j<leny;j++) {
		     if(i<lenx) {
			n3rx2 = LOG_ADD(n3rx2, Fn3[i][j]+e1+log_odds_x(i+1)+Brx2[i+1][j]);	
			rx2rx2 = LOG_ADD(rx2rx2, Frx2[i][j]+e1+log_odds_x(i+1)+Brx2[i+1][j]);	
		     }
		     if(i==lenx) {
		     	n3n4 = LOG_ADD(n3n4, Fn3[i][j]+e+Bn4[j]);
		     	rx2n4 = LOG_ADD(rx2n4, Frx2[i][j]+e+Bn4[j]);
		     }
		}
	}
	n3rx2 -= BE;
	rx2rx2 -= BE;
	n3n4 -= BE;
	rx2n4 -= BE;

	for(j=1;j<leny;j++) {
		if(j==leny) n4e = LOG_ADD(n4e, Fn4[lenx][j]+e);
		if(j<leny) {
			ry2ry2 = LOG_ADD(ry2ry2, Fry2[lenx][j]+e1+log_odds_y(j+1)+Bry2[j+1]);
			n4ry2 = LOG_ADD(n4ry2, Fn4[lenx][j]+e1+log_odds_y(j+1)+Bry2[j+1]);
		}
	}
	ry2ry2 -= BE;
	n4e -= BE;
	n4ry2 -= BE;

	ry2e = Fry2[lenx][leny]+e;
	ry2e -= BE;

	int index=0;
	cout << " rx1rx1 " << "\t" << exp(rx1rx1) << endl;
	cout << " rx1n1 " << "\t" << exp(rx1n1) << endl;
	cout << " n1ry1 " << "\t" << exp(n1ry1) << endl;
	cout << " n1n2 " << "\t" << exp(n1n2) << endl;
	cout << " ry1ry1 " << "\t" << exp(ry1ry1) << endl;
	cout << " ry1n2 " << "\t" << exp(ry1n2) << endl;
	cout << " n2m " << "\t" << exp(n2m) << endl;
	cout << " n2x " << "\t" << exp(n2x) << endl;
	cout << " n2y " << "\t" << exp(n2y) << endl;
	cout << " n2n3 " << "\t" << exp(n2n3) << endl;
	cout << " mm " << "\t" << exp(mm) << endl;
	cout << " mx " << "\t" << exp(mx) << endl;
	cout << " mn3 " << "\t" << exp(mn3) << endl;
	cout << " xx " << "\t" << exp(xx) << endl;
	cout << " xm " << "\t" << exp(xm) << endl;
	cout << " xn3 " << "\t" << exp(xn3) << endl;
	cout << " yy " << "\t" << exp(yy) << endl;
	cout << " ym " << "\t" << exp(ym) << endl;
	cout << " yn3 " << "\t" << exp(yn3) << endl;
	cout << " n3rx2 " << "\t" << exp(n3rx2) << endl;
	cout << " n3n4 " << "\t" << exp(n3n4) << endl;
	cout << " rx2rx2 " << "\t" << exp(rx2rx2) << endl;
	cout << " rx2n4 " << "\t" << exp(rx2n4) << endl;
	cout << " n4e " << "\t" << exp(n4e) << endl;
	cout << " n4ry2 " << "\t" << exp(n4ry2) << endl;
	cout << " ry2ry2 " << "\t" << exp(ry2ry2) << endl;
	cout << " ry2e " << "\t" << exp(ry2e) << endl;


	countE = LOG_ADD(rx1n1, n1n2, ry1n2, n3n4, rx2n4);
	countE1 = LOG_ADD(rx1rx1, n1ry1, ry1ry1, n3rx2, rx2rx2);
	estimatedE= countE-LOG_ADD(countE, countE1);
	cout << "Parameter esitmation for e: " << exp(e)  << "\t" << exp(estimatedE) << endl;

	countD = LOG_ADD(n2x, n2y, mx, my);
	countD1 = LOG_ADD(n2m, n2n3, mm, mn3);
	estimatedD = countD-LOG_ADD(countD, countD1)-log(2.0);
	cout << "Parameter esitmation for d: " << exp(d) << "\t" << exp(estimatedD) << endl;

	countEP = LOG_ADD(xx, yy);
	countEP1 = LOG_ADD(xm, xn3, ym, yn3);
	estimatedEP = countEP-LOG_ADD(countEP, countEP1);
	cout << "Parameter esitmation for ep: " << exp(ep) << "\t" << exp(estimatedEP) << endl;

	countT = LOG_ADD(xn3, n2n3, mn3, yn3);
	countT1 = LOG_ADD(n2m, n2x, n2y, mm, mx, my);
	countT1 = LOG_ADD(countT1, xm, xx, ym, yy);
	estimatedT = countT-LOG_ADD(countT, countT1);
	cout << "Parameter esitmation for t: " << exp(t) << exp(estimatedT) << endl;

	cout << "CountE and CountE1: " << exp(countE) << "\t" << exp(countE1) << endl;
	cout << "CountD and CountD1: " << exp(countD) << "\t" << exp(countD1) << endl;
	cout << "CountEP and CountEP1: " << exp(countEP) << "\t" << exp(countEP1) << endl;
	cout << "CountT and CountT1: " << exp(countT) << "\t" << exp(countT1) << endl;

}

void hmm_local::getPosterior() {
	
	; // nothing here now

}

/*float ne= 0.207316;
float nep = 0.863795;
float nt = 0.00897139;
float nd = 0.0241507;
*/
/*
float ne = 0.230949;
float nd = 0.0384072;
float nep = 0.810831;
float nt = 0.0101069;
*/

	float ne = 0.80;
	float nd = 0.0668055;
	float nep = 0.708576;
	float nt = 0.00965374;

void trainingLocalHMM(char *fastaFileList) {

	int i,j,k,m;

	vector<string> fastaFileNames;
	char tmpChar[1000], tmpChar1[5000], tmpChar2[100];
	string tmpString;

	ifstream fastaList(fastaFileList, ios::in);
	while(fastaList.good() ) {
		fastaList.getline(tmpChar, 1000);
		tmpString = tmpChar;
		fastaFileNames.push_back(tmpString);
	}
	fastaList.close();

	float e, e1, d, d1, ep, ep1, t, t1;
	e=e1=d=d1=ep=ep1=t=t1=0;

	float BE = 0;

	cout << "fastaFileNames size: " << fastaFileNames.size() << endl;
	for(i=0;i<fastaFileNames.size()-1;i++) {
	    strcpy(tmpChar, fastaFileNames[i].c_str() );
	    cout << "Fasta file name: " << tmpChar << endl;
	    sequences seqs(tmpChar, 1);
	    cout << "seqnum: " << seqs.nseqs << endl;
	    vector<subalign *> alns;
	    for(j=1;j<=seqs.nseqs;j++) {
		strcpy(tmpChar1, seqs.seq[j].c_str() );
		strcpy(tmpChar2, seqs.name[j].c_str() );
		alns.push_back(oneSeq2subalign(tmpChar1, tmpChar2) );
		alns[j-1]->gap_threshold = 1;
		alns[j-1]->beta=0;
		cout << alns[j-1]->nal << "\t" << alns[j-1]->alilen << endl;
		alns[j-1]->profile();
	    }
	    cout << "alns size: " << alns.size() << endl;
	    for(k=0;k<alns.size();k++) {
	       	for(m=k+1;m<alns.size();m++) {
			cout << "HMM Training: " << i << " " << k << "  " << m << endl;
			hmm_local trainHmm(alns[k], alns[m]);
			trainHmm.eta = ne;
			trainHmm.epsilon = nep;
			trainHmm.delta = nd;
			trainHmm.tau = nt;
			trainHmm.forward();
			trainHmm.backward();
			trainHmm.baumwelch();
			BE += trainHmm.BE;
			e += exp(trainHmm.countE);
			e1 += exp(trainHmm.countE1);
			ep += exp(trainHmm.countEP);
			ep1 += exp(trainHmm.countEP1);
			d += exp(trainHmm.countD);
			d1 += exp(trainHmm.countD1);
			t += exp(trainHmm.countT);
			t1 += exp(trainHmm.countT1);

			cout << "New eta: " << e << endl;
			cout << "New eta1: " << e1 << endl;
			cout << "New delta: " << d << endl;
			cout << "New delta1: " << d1 << endl;
			cout << "New epsilon: " << ep << endl;
			cout << "New epsilon1: " << ep1 << endl;
			cout << "New tau: " << t << endl;
			cout << "New tau1: " << t1 << endl;

	       	}
	    }
	    for(j=1;j<=seqs.nseqs;j++) { delete alns[j-1]; }
	    alns.erase(alns.begin(), alns.end() );
	}		

	ne = e/(e+e1);
	nd = d/(d+d1)/2;
	nep = ep/(ep+ep1);
	nt = t/(t+t1);

	cout << "NEW eta: " << ne << endl;
	cout << "NEW delta: " << nd << endl;
	cout << "NEW epsilon: " << nep << endl;
	cout << "NEW tau: " << nt << endl;
	cout << "Cumulative BE:"<<BE <<endl;
}
