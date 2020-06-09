#include "hmm_multim.h"
#include "mm.h"

static int debug = -1;

void hmm_parameters::read_parameters(char *filename) {

	int i, j, k;
	char line[500];

	if (debug> 1 )cout << "Number of states: " << num_states << endl;
	if (debug> 1) cout << "Number of match states: " << num_match_states << endl;

	// allocations
	transition_matrix = gmatrix<float>(num_states, num_states);
	indel_emission_vector1 = gvector<float>(20);
	indel_emission_vector2 = gvector<float>(20);
	match_emission_matrix = new float **[num_match_states +1];
	background_frequencies = gvector<float>(20);
	for(i=1;i<=num_match_states;i++) {
		match_emission_matrix[i] = gmatrix<float>(20, 20);
	}

	t_log = gmatrix<float>(num_states, num_states);
	indel_lor1 = gvector<float>(20);
	indel_lor2 = gvector<float>(20);
	m_lor = new float **[num_match_states +1];
	indel_log1 = gvector<float>(20);
	indel_log2 = gvector<float>(20);
	m_log = new float **[num_match_states +1];
	background_frequencies = gvector<float>(20);
	for(i=1;i<=num_match_states;i++) {
		m_lor[i] = gmatrix<float>(20, 20);
		m_log[i] = gmatrix<float>(20, 20);
	}
	ifstream fp(filename, ios::in);
	if (!fp) {
		cout << "Error: input file - " << filename << " not readable"<< endl;
		exit(0);
	}

	int ss1, solv1, unaligned1;
	char tmpstr[50];

	while(fp.good()){
	    fp.getline(line, 500);
	    if(strncmp(line, "ss", 2)==0) {
		break;
	    }
	}
	sscanf(line, "%s %d", tmpstr, &ss1);
	if(num_ss != ss1) {
	    cout << "Error: numbers of secondary structure types do not match" << endl;
	    cout << "       in parameter file and input option (-ss)" << endl;
	    exit(0);
	}

	while(fp.good()){
	    fp.getline(line, 500);
	    if(strncmp(line, "solv", 4)==0) {
		break;
	    }
	}
	sscanf(line, "%s%d", tmpstr, &solv1);
	if(num_solv != solv1) {
	    cout << "Error: numbers of solvent accessibility categories do not match" << endl;
	    cout << "       in parameter file and input option (-solv)" << endl;
	    exit(0);
	}

	while(fp.good()){
	    fp.getline(line, 500);
	    if(strncmp(line, "unaligned", 9)==0) {
		break;
	    }
	}
	sscanf(line, "%s %d", tmpstr, &unaligned1);
	if(unaligned1 != unaligned) {
	    cout << "Error: numbers of unaligned match states do not match" << endl;
	    cout << "       in parameter file and input option (-unaligned)" << endl;
	    exit(0);
	}

	while(fp.good()) {
		fp.getline(line, 500);
		if(strncmp(line, "Transition", 10)==0) {
			break;
		}
	}
	for(i=0;i<num_states;i++) {
		for(j=0;j<num_states;j++) {
			fp >> transition_matrix[i][j];
			//cout << i << " " << j << " " <<  transition_matrix[i][j] << " ";
		}
		//cout << endl;
	}
	
	while(fp.good()) {
		fp.getline(line, 500);
		if(strncmp(line, "Emission vector for indel state 1:",34)==0){
			break;
		}
	}
	for(i=0;i<=20;i++) {
		fp >> indel_emission_vector1[i];
	}
	while(fp.good()) {
		fp.getline(line, 500);
		if(strncmp(line, "Emission vector for indel state 2:",34)==0){
			break;
		}
	}
	for(i=0;i<=20;i++) {
		fp >> indel_emission_vector2[i];
	}

	while(fp.good()) {
		fp.getline(line, 500);
		if(strncmp(line, "Background frequencies:",23)==0){
			break;
		}
	}
	for(i=0;i<=20;i++) {
		fp >> background_frequencies[i];
	}

	while(fp.good()) {
		fp.getline(line, 500);
		if(strncmp(line, "Emission matrices:", 18)==0) {
			break;
		}
	}
	for(i=1;i<=num_match_states;i++) {
		fp.getline(line, 500);
		//cout << "===== " << line << endl;
		for(j=0;j<=20;j++) {
			for(k=0;k<=20;k++) {
				fp >> match_emission_matrix[i][j][k];
				//cout << setw(7) << match_emission_matrix[i][j][k];
			}
			//cout << endl;
		}
		fp.getline(line, 500);
		fp.getline(line, 500);
		//cout << "===++ " << line << endl;
	}

	// transform to log scale
	for(i=1;i<=20;i++) {
		indel_lor1[i] = log(indel_emission_vector1[i]/background_frequencies[i]);
		indel_lor2[i] = log(indel_emission_vector2[i]/background_frequencies[i]);
		indel_log1[i] = log(indel_emission_vector1[i]);
		indel_log2[i] = log(indel_emission_vector2[i]);
	}
	
	// transform to log scale
	for(i=1;i<=num_match_states;i++) {
	   for(j=1;j<=20;j++) {
		for(k=1;k<=20;k++) {
		  m_lor[i][j][k] = log(match_emission_matrix[i][j][k]/background_frequencies[j]/background_frequencies[k]);
		  m_log[i][j][k] = log(match_emission_matrix[i][j][k]);
		}
	   }
	}

	// transform to log scale
	for(i=0;i<num_states;i++) {
		for(j=0;j<num_states;j++) {
			if(transition_matrix[i][j]==0) {
				t_log[i][j] = -2e20;
			}
			else {
				t_log[i][j] = log(transition_matrix[i][j]);
			}
			// MAKING the transition probability to the ending state the same
			//if(j==num_states-1) t_log[i][j] = log(1.0/num_states);
		}
		//if(i==0) { for(j=1;j<=num_match_states+2;j++) t_log[i][j] = log(1.0/num_states); }
	}
	
	fp.close();
}

void hmm_parameters::print_parameters() {

	int i,j,k;
	cout << fixed;
	cout << setprecision(3);
	cout << setw(7);
	cout << left;
	
	cout << endl;

	cout << "Transition matrix:" << endl;
	for(i=0;i<num_states;i++) {
		for(j=0;j<num_states;j++) {
			cout << setw(9) << transition_matrix[i][j];
		}
		cout << endl;
	}
	cout <<endl;
	cout << "Transition matrix in LOG:" << endl;
	for(i=0;i<num_states;i++) {
		for(j=0;j<num_states;j++) {
			cout << setw(9) << setprecision(2) << t_log[i][j];
		}
		cout << endl;
	}
	cout <<endl;
	cout << "Emission vector for indel state 1:" << endl;
	for(i=1;i<=20;i++) {
		cout << am[i] << " " << indel_emission_vector1[i]<<" "<< indel_lor1[i]<<endl;
	}
	cout << endl;
	cout << "Emission vector for indel state 2:" << endl;
	for(i=1;i<=20;i++) {
		cout << am[i] << " " << indel_emission_vector2[i]<<" "<<indel_lor2[i]<<endl;
	}
	cout << endl;
	cout << setprecision(4);
	cout << "Emission matrices:" << endl;
	for(i=1;i<=num_match_states;i++) {
		cout << "Matrix: " << i << endl;
		for(j=1;j<=20;j++) {
			for(k=1;k<=20;k++) {
				cout << setw(7) << match_emission_matrix[i][j][k];
			}
			cout << endl;
		}
		cout << endl;
	}

	cout << endl;
	cout << "Emission matrices in LOR:" << endl;
	for(i=1;i<=num_match_states;i++) {
		cout << "Matrix: " << i << endl;
		for(j=1;j<=20;j++) {
			for(k=1;k<=20;k++) {
				cout << setw(7) << setprecision(3) << m_lor[i][j][k];
			}
			cout << endl;
		}
		cout << endl;
	}


}

hmm_multim::~hmm_multim() {

	int i;

	if (F) {
		for(i=1;i<=p->num_match_states +2;i++) {
			free_gmatrix<ScoreType>(F[i], lenx, leny);
		}
		delete [] F;
	}
	if (B) {
		for(i=1;i<=p->num_match_states +2;i++) {
			free_gmatrix<ScoreType>(B[i], lenx, leny);
		}
		delete [] B;
	}
	if (probMat) {
		free_gmatrix<float>(probMat, lenx, leny);
	}
		
}

void hmm_multim::viterbi(int lor) {

	int i,j,k,m;
	int *ax = x->alignment[1];
	int *ay = y->alignment[1];
	vector<float> vv;
	vector<float>::iterator it;
	
	// allocation
	V = new float **[p->num_match_states +3];
	for(i=1;i<=p->num_match_states +2;i++) {
		V[i] = gmatrix<float>(lenx, leny);
	}
	T = new int **[p->num_match_states +3];
	for(i=1;i<=p->num_match_states +2;i++) {
		T[i] = gmatrix<int>(lenx, leny);
	}

	//cout<< "IN viterbi here =============" << endl;
	cout << "Viterbi algorithm:" << endl;

	for(i=0;i<=lenx;i++) {
	   for(j=0;j<=leny;j++) {
	      	// match states
		//cout << "======== " << i << "  " << j << endl;
		for(k=1;k<=num_match_states;k++) {
		    if( (i==0) || (j==0) ) {
			V[k][i][j] = LOG_ZERO;
			T[k][i][j] = -1000;
			if(debug>1) fprintf(stdout, "%d %d %d:  %7.3f    %d", k, i, j, V[k][i][j], T[k][i][j]);
			if(debug>1) cout << endl;
			continue;
		    }
		    if( (i==1)&&(j==1) ) {
			V[k][i][j] = p->t_log[0][k];
			T[k][i][j] = 0;
			if(lor==1) V[k][i][j] += p->m_lor[k][ax[i]][ay[j]];
			else V[k][i][j] += p->m_log[k][ax[i]][ay[j]];
			if(debug>1) fprintf(stdout, "%d %d %d:  %7.3f    %d", k, i, j, V[k][i][j], T[k][i][j]);
			if(debug>1) cout << endl;
			continue;
		    }
		    vv.clear();
		    for(m=1;m<=num_match_states+2;m++) {
			vv.push_back(p->t_log[m][k] + V[m][i-1][j-1]);
		    }
		    it = max_element(vv.begin(), vv.end());
		    if(lor==1) V[k][i][j]=p->m_lor[k][ax[i]][ay[j]]+(*it);
		    else V[k][i][j]=p->m_log[k][ax[i]][ay[j]]+(*it);
		    vv.erase(++it, vv.end());
		    T[k][i][j]=vv.size();
		    if(debug>1) fprintf(stdout, "%d %d %d:  %7.3f    %d", k, i, j, V[k][i][j], T[k][i][j]);
		    if(debug>1) cout << endl;
	      }
		// insertion state 1
	     //cout << "+++++++++++" << endl;
 	     if(i==0) { V[p1][i][j] = LOG_ZERO; }
	     else {
		if( (j==0)&&(i==1) ) {
			V[p1][i][j] = p->t_log[0][p1];
			T[p1][i][j] = 0;
			if(lor==1) V[p1][i][j] += p->indel_lor1[ax[i]];
			else V[p1][i][j] += p->indel_log1[ax[i]];
			if(debug>1) fprintf(stdout, "%d %d %d:  %7.3f    %d", p1, i, j, V[p1][i][j], T[p1][i][j]);
			if(debug>1) cout << endl;
		}
		else {
			vv.clear();
			for(m=1;m<=num_match_states+2;m++) {
				vv.push_back(p->t_log[m][p1] + V[m][i-1][j]);
				//cout << p1 << " " << m << " " << i-1 << " " << j  << " " <<  p->t_log[m][p1] << " " << V[m][i-1][j] << endl;
			}
			it = max_element(vv.begin(), vv.end());
			if(lor==1) V[p1][i][j]=p->indel_lor1[ax[i]] + (*it);
			else V[p1][i][j]=p->indel_log1[ax[i]] + (*it);
			vv.erase(++it, vv.end());
			T[p1][i][j]=vv.size();
			//fprintf(stdout, "%d %d %d:  %7.3f    %d", p1, i, j, V[p1][i][j], T[p1][i][j]);
			//cout << endl;
		}
	     }
	     //cout << "+++++++++++" << endl;
		
		// insertion state 2
	     if(j==0) { V[p2][i][j] = LOG_ZERO; }
	     else {
		if( (i==0)&&(j==1) ) {
			V[p2][i][j] = p->t_log[0][p2];
			if(lor==1) V[p2][i][j] += p->indel_lor2[ay[j]];
			else V[p2][i][j] += p->indel_log2[ay[j]];
			T[p2][i][j] = -1000;
			if(debug>1) fprintf(stdout, "%d %d %d:  %7.3f    %d", p2, i, j, V[p2][i][j], T[p2][i][j]);
			if(debug>1) cout << endl;
		}
		else {
			vv.clear();
			for(m=1;m<=num_match_states+2;m++) {
				vv.push_back(p->t_log[m][p2]+V[m][i][j-1]);
				//cout << m << " " << p2 << " " << i << " " << j-1  << " " <<  p->t_log[m][p2] << " " << V[m][i][j-1] << endl;
			}
			//for(it=vv.begin();it<vv.end();it++) { cout << *it << " "; }
			//cout << endl;
			it=max_element(vv.begin(), vv.end());
			if(lor==1) V[p2][i][j]=p->indel_lor2[ay[j]] + (*it);
			else V[p2][i][j]=p->indel_log2[ay[j]] + (*it);
			vv.erase(++it, vv.end());
			T[p2][i][j]=vv.size();
			if(debug>1) fprintf(stdout, "%d %d %d:  %7.3f    %d", p2, i, j, V[p2][i][j], T[p2][i][j]);
			if(debug>1) cout << endl;
		}
	    }
	  }
	}
	if(debug>1) cout<< "IN viterbi here =============" << endl;

	// the ending state
	vv.clear();
	for(m=1;m<=num_match_states+2;m++) {
		vv.push_back(p->t_log[m][p_end]+V[m][lenx][leny]);
		if(debug>1) cout << m << " " << V[m][lenx][leny] << " " << p->t_log[m][p_end] << " " << V[m][lenx][leny]+p->t_log[m][p_end] << endl;
	}
	it=max_element(vv.begin(), vv.end());
	VE = (*it);
	vv.erase(++it, vv.end());
	TE = vv.size();
	vv.clear();

	if(debug>1) cout << "VE and TE: " << VE << "  " << TE << endl;

	// trace back
	int *trace = new int[lenx+leny+2];
	trace[0] = TE;
	int px = lenx; 
	int py = leny;
	//cout << "Tracing result: " << endl;
	for(i=1;i<=lenx+leny;i++) {
		trace[i] = T[trace[i-1]][px][py];
		if (trace[i]==-1000) { break;}
		if(trace[i-1]<=num_match_states) {
			px--;
			py--;
		}
		else if(trace[i-1]==num_match_states+1) {
			px--;
		}
		else if(trace[i-1]==num_match_states+2) {
			py--;
		}
		if(debug>1) cout << i << "  " << trace[i] << endl;
	}
	int trace_len = i;
	string sx, sy;
	sx = "";
	sy = "";
	k = 0; m = 0;
	for(i=0;i<trace_len;i++) {
		j = trace_len - i-1;
		if(trace[j]<=num_match_states) {
			sx.append(1, x->aseq[0][k]);
			sy.append(1, y->aseq[0][m]);
			k++;
			m++;
		}
		if(trace[j]==num_match_states+1) {
			sx.append(1, x->aseq[0][k]);
			sy.append(1, '-');
			k++;
		}
		if(trace[j]==num_match_states+2) {
			sx.append(1, '-');
			sy.append(1, y->aseq[0][m]);
			m++;
		}
	}
	cout << left << setw(20) << x->aname[0] << "  " << sx << endl;
	cout << left << setw(20) << y->aname[0] << "  " << sy << endl;
	cout << endl;

	// clear up
	for(i=1;i<=p->num_match_states +2;i++) {
		free_gmatrix<float>(V[i], lenx, leny);
		free_gmatrix<int>(T[i], lenx, leny);
	}
	delete [] V;
	delete [] T;
	delete [] trace;
}


void hmm_multim::forward() {

	int i,j,k,m;
	int *ax = x->alignment[1];
	int *ay = y->alignment[1];

	F = new ScoreType **[p->num_match_states+3];
	for(i=1;i<=p->num_match_states +2;i++) {
		F[i] = gmatrix<ScoreType>(lenx, leny);
	}

	if(debug>1) cout << "Forward algorithm:" << endl;

	for(i=0;i<=lenx;i++) {
	    for(j=0;j<=leny;j++) {
		// match states
		for(k=1;k<=p->num_match_states;k++) {
			if( (i==0)||(j==0) ) {
			    	F[k][i][j] = LOG_ZERO;
				continue;
			}
			if( (i==1) && (j==1) ) {
				F[k][i][j] = p->t_log[0][k];
				F[k][i][j] += p->m_log[k][ax[i]][ay[j]];
				if(debug>-1)fprintf(stdout, "%d %d %d:  %7.3f", k, i, j, F[k][i][j]);
				if(debug>-1)cout << endl;
				continue;
			}
			F[k][i][j] = LOG_ZERO;
			for(m=1;m<=p->num_match_states+2;m++) {
				F[k][i][j] = LOG_ADD(F[k][i][j], p->t_log[m][k]+F[m][i-1][j-1]);
			}
			F[k][i][j] += p->m_log[k][ax[i]][ay[j]];
			if(debug>-1)fprintf(stdout, "%d %d %d:  %7.3f", k, i, j, F[k][i][j]);
			if(debug>-1)cout << endl;
		}

		// insertion state 1
		if(i==0) { F[p1][i][j] = LOG_ZERO; }
		else {
		  if( (j==0)&&(i==1) ) {
			F[p1][i][j] = p->t_log[0][p1];
			F[p1][i][j] += p->indel_log1[ax[i]];
			if(debug>-1)fprintf(stdout, "%d %d %d:  %7.3f", p1, i, j, F[p1][i][j]);
			if(debug>-1)cout << endl;
		  }
		  else {
			F[p1][i][j] = LOG_ZERO;
			for(m=1;m<=p->num_match_states+2;m++) {
				F[p1][i][j] = LOG_ADD(F[p1][i][j], p->t_log[m][p1]+F[m][i-1][j]);
			}
			F[p1][i][j] += p->indel_log1[ax[i]];
			if(debug>-1)fprintf(stdout, "%d %d %d:  %7.3f", p1, i, j, F[p1][i][j]);
			if(debug>-1)cout << endl;
		   }
		}
			
		// insertion state 2
		if(j==0) { F[p2][i][j] = LOG_ZERO; }
		else {
		  if( (i==0)&&(j==1) ) {
			F[p2][i][j] = p->t_log[0][p2];
			F[p2][i][j] += p->indel_log2[ay[j]];
			if(debug>-1)fprintf(stdout, "%d %d %d:  %7.3f", p2, i, j, F[p2][i][j]);
			if(debug>-1)cout << endl;
		  }
		  else {
			F[p2][i][j] = LOG_ZERO;
			for(m=1;m<=p->num_match_states+2;m++) {
				F[p2][i][j] = LOG_ADD(F[p2][i][j], p->t_log[m][p2]+F[m][i][j-1]);
			}
			F[p2][i][j] += p->indel_log2[ay[j]];
			if(debug>-1)fprintf(stdout, "%d %d %d:  %7.3f", p2, i, j, F[p2][i][j]);
			if(debug>-1)cout << endl;
		   }
		}
	    }
	}
			
	// the ending state
	FE = LOG_ZERO;
	for(m=1;m<=num_match_states+2;m++) {
		if(debug>-1) fprintf(stdout, "%d %f %f\n", m, p->t_log[m][p_end], p->transition_matrix[m][p_end]);
		FE = LOG_ADD(FE, p->t_log[m][p_end]+F[m][lenx][leny]);
	}
	if(debug>1) cout << "FE: " << FE << endl;
}

void hmm_multim::backward() {

	int i,j,k,m;
	int *ax = x->alignment[1];
	int *ay = y->alignment[1];

	B = new ScoreType **[p->num_match_states+3];
	for(i=1;i<=p->num_match_states+2;i++) {
		B[i] = gmatrix<ScoreType>(lenx, leny);
	}

	for(i=lenx;i>=0;i--) {
	    for(j=leny;j>=0;j--) {
		//match states and insertion states together
		for(k=1;k<=p->num_match_states+2;k++) {
			if( (i==lenx) && (j==leny) ) {
				B[k][i][j] = p->t_log[k][p_end];
				if(debug>-1) cout << k << "\t" << i << "\t" << j << "\t" << B[k][i][j] << endl;
				continue;
			}
			B[k][i][j] = LOG_ZERO;
			if( i==lenx) {
			    B[k][i][j] = LOG_ADD(B[k][i][j], p->t_log[k][p2]+B[p2][i][j+1] + p->indel_log2[ay[j+1]]);
				if(debug>-1)cout << k << "\t" << i << "\t" << j << "\t" << B[k][i][j] << endl;
			    continue;
			}
			if(j==leny) {
			    B[k][i][j] = LOG_ADD(B[k][i][j], p->t_log[k][p1]+B[p1][i+1][j] + p->indel_log1[ax[i+1]]);
				if(debug>-1)cout << k << "\t" << i << "\t" << j << "\t" << B[k][i][j] << endl;
			    continue;
			}
			
			for(m=1;m<=num_match_states;m++) {
				//cout << m << "\t" << p->t_log[k][m] << "\t" << B[m][i+1][j+1] << "\t" << p->m_log[k][ax[i+1]][ay[j+1]] << endl;
				B[k][i][j] = LOG_ADD(B[k][i][j], p->t_log[k][m]+B[m][i+1][j+1]+p->m_log[m][ax[i+1]][ay[j+1]]);
			} 
			B[k][i][j] = LOG_ADD(B[k][i][j], p->t_log[k][p1]+B[p1][i+1][j] + p->indel_log1[ax[i+1]]);
			B[k][i][j] = LOG_ADD(B[k][i][j], p->t_log[k][p2]+B[p2][i][j+1] + p->indel_log2[ay[j+1]]);
				if(debug>-1)cout << k << "\t" << i << "\t" << j << "\t" << B[k][i][j] << endl;
		}
	    }
	}

	BE = LOG_ZERO;
	for(m=1;m<=p->num_match_states;m++) {
		BE = LOG_ADD(BE, p->t_log[0][m]+B[m][1][1]+p->m_log[m][ax[1]][ay[1]]);
	}
	//fprintf(stdout, "BE: %f", BE);
	BE = LOG_ADD(BE, p->t_log[0][p1]+B[p1][1][0]+p->indel_log1[ax[1]]);
	BE = LOG_ADD(BE, p->t_log[0][p2]+B[p2][0][1]+p->indel_log2[ay[1]]);

	if(debug>1) fprintf(stdout, "BE: %f", BE);
	if(debug>1) cout << endl;

	//ScoreType full_prob;
	if(debug>1) cout << "Obtain posterior probabilities (>0.01): " << x->aname[0] << "\t" << y->aname[0] << endl;
	ScoreType **aligned_pair_prob = gmatrix<ScoreType>(lenx, leny);
	probMat = gmatrix<float>(lenx, leny);
	int useful_num_match_states=num_match_states-unaligned;
	//int useful_num_match_states=num_match_states;
	for(i=1;i<=lenx;i++) {
		for(j=1;j<=leny;j++) {
			aligned_pair_prob[i][j] = LOG_ZERO;
			for(k=1;k<=useful_num_match_states;k++) {
				aligned_pair_prob[i][j]=LOG_ADD(aligned_pair_prob[i][j],F[k][i][j]+B[k][i][j]-FE);
			}
			probMat[i][j] = exp(aligned_pair_prob[i][j]);
			if(debug>-2) {
				if((probMat[i][j]>-0.01) && (probMat[i][j]<1.1))
					if(debug>-1) fprintf(stdout, "%d %d: %c %c %f %f %f\n", i, j, x->aseq[0][i-1], y->aseq[0][j-1], aligned_pair_prob[i][j], probMat[i][j], exp(F[num_match_states][i][j]+B[num_match_states][i][j]-FE) );
				}
		}
	}

	free_gmatrix<ScoreType>(aligned_pair_prob, lenx, leny);

	ScoreType tmp;

	if(debug>-1)
	for(i=1;i<=lenx;i++) {
		for(j=1;j<=leny;j++) {
			tmp = LOG_ZERO;
			for(k=1;k<=num_match_states+2;k++) {
				tmp = LOG_ADD(tmp, B[k][i][j]+F[k][i][j]);
			}
			cout << i << "\t" << j << "\t" << tmp << endl;
		}
	}

	//fprintf(stdout, "================\n");

}
	
int * hmm_multim::posteriorAlignment() {

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
	a = new char[3000];
	b = new char[3000];

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
	
	cout << left << setw(40) << x->aname[0] << a << endl;
	cout << left << setw(40) << y->aname[0] << b << endl;

	delete [] a;
	delete [] b;

	return path;
}

