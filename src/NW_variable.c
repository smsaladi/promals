#include "NW_variable.h"
#include "header_cpp.h"

NW_variable::NW_variable() {
	
	D = NULL;
	C = NULL;
	I = NULL;
	T = NULL;
	T_I = NULL;
	T_D = NULL;
	tr = NULL;
	r_tr = NULL;
	w = NULL;

	use_end_penalty = 0;
	debug = 0;

}

NW_variable::NW_variable(int m, int n, double *_go1, double *_ge1, double *_go2, double *_ge2) {

	int i;

	M = m;
	N = n;
	go1  = _go1;
	ge1  = _ge1;
	go2  = _go2;
	ge2  = _ge2;
	
	w = NULL;

	if(g>0) warning("gap opening penalty should be a negative integer");
	if(h>0) warning("gap extension penalty should be a negative integer");

	use_end_penalty = 0;
	debug = 0;

	D = new double * [M+1];
	C = new double * [M+1];
	I = new double * [M+1];
	T = new int * [M+1];
	T_I = new int * [M+1];
	T_D = new int * [M+1];
	for(i=0;i<=M;i++) {
	   D[i] = new double [N+1];
	   C[i] = new double [N+1];
	   I[i] = new double [N+1];
	   T[i] = new int [N+1];
	   T_I[i] = new int [N+1];
	   T_D[i] = new int [N+1];
	}
	tr = new int [M+N+2];
	r_tr = new int [M+N+2];
}

NW_variable::~NW_variable() {

	int i;

	if(D) {
	   for(i=0;i<=M;i++) delete [] D[i];
	   delete [] D;
	}
	if(C) {
           for(i=0;i<=M;i++) delete [] C[i];
           delete [] C;
        }
	if(I) {
           for(i=0;i<=M;i++) delete [] I[i];
           delete [] I;
        }
	if(T) {
	   for(i=0;i<=M;i++) delete [] T[i];
           delete [] T;
	}
       	if(T_I) {
           for(i=0;i<=M;i++) delete [] T_I[i];
           delete [] T_I;
        }
	if(T_D) {
	   for(i=0;i<=M;i++) delete [] T_D[i];
	   delete [] T_D;
	}
	if(tr) delete [] tr;
	if(r_tr) delete [] r_tr;
}

void NW_variable::setM(int x) {
	M = x;
}

void NW_variable::setN(int x) {
 	N = x;
}

void NW_variable::set_g(double x) {
	g = x;
	if(g>0) warning("gap opening penalty should be a negative integer");
}

void NW_variable::set_go2(double *x) {
	go2 = x;
}

void NW_variable::set_ge2(double *x) {
	ge2 = x;
}

void NW_variable::set_go1(double *x) {
	go1 = x;
}

void NW_variable::set_ge1(double *x) {
	ge1 = x;
}
void NW_variable::set_h(double x) {
	h = x;
	if(h>0) warning("gap extension penalty should be a negative integer");
}

void NW_variable::set_use_end_penalty(int x) {
	use_end_penalty = x;
}

void NW_variable::set_w(double **sm) {
	w = sm;
}

// x: I[-1]; y: C[0]
// Insertion and deletion in the second sequence
// I:    i -        D:  i-1 i
//     j-1 j              j -
double NW_variable::MAX_I(double x, double y) {

	if(x<y) {
		tracedirection_I = 0;
		return y;
	}
	else {
		tracedirection_I = 1;
		return x;
	}

}

// x: D[-1]; y: C[0]
double NW_variable::MAX_D(double x, double y) {

        if(x<y) {
		tracedirection_D = 0;
                return y;
        }
        else {
		tracedirection_D = -1;
                return x;
        }

}

// x: D [-1]; y: I [1]; z: C [0]
double NW_variable::MAX3(double x, double y, double z) {

	if(x>y) {
	   if(x>z) {
		tracedirection=-1;
		return x;
	   }
	   else {
		tracedirection = 0;
		//tracedirection_D=tracedirection_I=0;
		return z;
	   }
	}
	else {
	   if(y>z) {
		tracedirection = 1;
		return y;
	   }
	   else {
		tracedirection = 0;
		//tracedirection_D=tracedirection_I=0;
		return z;
	   }
	}
}

void NW_variable::warning(char *s) {

     fprintf(stderr, "%s\n", s);
}

void NW_variable::dp() {

	int i, j;
	double t;
	
	C[0][0] = 0;
	D[0][0] = 0;
	I[0][0] = 0;
	T[0][0] = 0;
	T_D[0][0] = 0;
	T_I[0][0] = 0;
	t = g = h = 0;

	tracedirection = tracedirection_I = tracedirection_D = 0;

	//for(i=0;i<=M;i++) {go1[i] = 5; ge1[i] = 1;}
	//for(i=0;i<=N;i++) {go2[i] = 5; ge2[i] = 1;}
	
	i = 0;
	for(j=1;j<=N;j++) {
		C[0][j] = C[0][j-1] - 10;
		I[0][j] = I[0][j-1] - 10;
		D[0][j] = D[0][j-1] - 10;
		T[0][j] = 1;
		T_D[0][j] = T_I[0][j] = 1;
		fprintf(stdout, "%d %d: %f %f %f %d %d %d\n", i, j, C[i][j], D[i][j],I[i][j],T[i][j],T_D[i][j],T_I[i][j]);
	}

	for(i=1;i<=M;i++) {
		j = 0;
		C[i][0] = C[i-1][0] -10;
		D[i][0] = D[i-1][0] -10;
		I[i][0] = I[i-1][0] -10;
		T[i][0] = T_D[i][0] = T_I[i][0] = -1;
		fprintf(stdout, "%d %d: %f %f %f %d %d %d\n", i, j, C[i][j], D[i][j],I[i][j],T[i][j],T_D[i][j],T_I[i][j]);
		
		for(j=1;j<=N;j++) {
		    	I[i][j] = MAX_I(I[i][j-1]-ge1[i], C[i][j-1]-go1[i]);
			D[i][j] = MAX_D(D[i-1][j]-ge2[j], C[i-1][j]-go2[j]);
			C[i][j] = MAX3(D[i-1][j-1], I[i-1][j-1], C[i-1][j-1])+w[i][j];
			T[i][j] = tracedirection;
			T_I[i][j] = tracedirection_I;
			T_D[i][j] = tracedirection_D;
			fprintf(stdout, "%d %d: %f %f %f %d %d %d %f %f %f %f\n", i, j, C[i][j], D[i][j],I[i][j],T[i][j],T_D[i][j],T_I[i][j], go1[i], ge1[i], go2[j], ge2[j]);
		}
	}
	
	/*
	tracedirection = tracedirection_I = tracedirection_D = 0;

	cout << "===============" << endl;

	for(j=1;j<=N;j++) {
	   if(use_end_penalty) {
	     C[0][j] = t = t+h;
	     D[0][j] = t+g;
	     T[0][j] = 1;
             T_I[0][j] = j;
             T_D[0][j] = 0;

	   }
	   else {
	     C[0][j] = 0;
	     D[0][j] = g;
             T[0][j] = 1;
             T_I[0][j] = j;
             T_D[0][j] = 0;
	   }
	}
	cout << "===============" << endl;
	
	t = g;

	for(i=1;i<=M;i++) {
		if(use_end_penalty) {
		   C[i][0] = t = t+h;
		   I[i][0] = t+g;
		}
		else {
		   C[i][0] = 0;
		   I[i][0] = g;
		}
		T[i][0] = -1;
		T_D[i][0] = i;
		T_I[i][0] = 0;
		//tracedirection_I = tracedirection_D = 0;
		for(j=1;j<=N;j++) {
		    	I[i][j] = MAX_I(I[i][j-1]-ge1[i], C[i][j-1]-go1[i]);
			D[i][j] = MAX_D(D[i-1][j]-ge2[j], C[i-1][j]-go2[j]);
			C[i][j] = MAX3(D[i][j], I[i][j], C[i-1][j-1]+w[i][j]);
			T[i][j] = tracedirection;
			if(tracedirection_I==0) T_I[i][j]=1;
			else T_I[i][j] = T_I[i][j-1]+tracedirection_I;
			if(tracedirection_D==0) T_D[i][j]=1;
			else T_D[i][j] = T_D[i-1][j]+tracedirection_D;
		}
	}
	*/

	cout << "===============" << endl;
	if(debug>-11) for(i=0;i<=M;i++) {
	   for(j=0;j<=N;j++) {
		cout << C[i][j] << "/" <<  T_I[i][j] << "/" << T_D[i][j] << "/" << T[i][j] << "\t";
	   }
	   cout << endl;
	}
	
}

void NW_variable::traceback() {
	
	int i, j, k, m;
	int px=-1, py=-1; // position maximum score resides in

	int end_direction;
	cout << C[M][N] << " "<< D[M][N]<<" " << I[M][N] << endl;
	if( (C[M][N]>=I[M][N]) && (C[M][N]>=D[M][N]) ) {
		maxscore = C[M][N];
		end_direction = 0;
	}
	else if( (I[M][N]>=C[M][N]) && (I[M][N]>=D[M][N]) ) {
		maxscore = I[M][N];
		end_direction = 1;
	}
	else if( (D[M][N]>=I[M][N]) && (D[M][N]>=C[M][N]) ) {
		maxscore = D[M][N];
		end_direction = -1;
	}

	cout << "end_direction " << end_direction << endl;

	i = M; j= N;
	pt = 0;
	while( (i>=0) || (j>=0) ) {
		
		if( (i==0) && (j==0) ) break;
		if(end_direction==0) {
			tr[pt] = 0;
			end_direction = T[i][j];
			cout << pt << "\t" << i << "\t"<<j<<"\t"<<tr[pt]<<endl;
			i--; j--;
		}
		else if(end_direction==1) {
			tr[pt] = 1;
			cout << pt << "\t" << i << "\t"<<j<<"\t"<<tr[pt]<<endl;
			end_direction = T_I[i][j];
			j--;
		}
		else if(end_direction==-1) {
			tr[pt] = -1;
			cout << pt << "\t" << i << "\t"<<j<<"\t"<<tr[pt]<<endl;
			end_direction = T_D[i][j];
			i--;
		}
		pt++;
	}
	// determine the reverse path
	for(i=pt-1;i>=0;i--) {
	    r_tr[pt-i] = tr[i];
	    if(debug>11) cout << tr[i] << " ";
	}

			
	/* 
			
	maxscore = C[M][N];

	pt = 0;
	// if no end gap penalty, determine the largest score ending with gaps
	if(!use_end_penalty) {
	    px = M;
	    for(i=M;i>=0;i--) {
		if(C[i][N]>maxscore) {
		    px = i;
		    maxscore = C[i][N];
		} 
	    }
	    for(j=N;j>=0;j--) {
		if(C[M][j]>maxscore) {
		   py = j;
		   maxscore = C[M][j];
		   px = -1;
		}
	    }
	
	    // deletion in the second dimension (x -)
	    if(px>=0) {
	        i = px; j = N;
		for(pt=0;pt<=M-px-1;pt++) {
		   tr[pt] = -1;
		}
	    }
	    // insertion in the second dimension (- y)
	    else if(py>=0) {
	        j = py; i = M;
		for(pt=0;pt<=N-py-1;pt++) {
		   tr[pt] = 1;
		}
	    }
	    else {
	        warning("ending problem");
	   	fprintf(stdout, "px %d py %d\n", px, py);
	    }
	}
  	else { i = M; j = N; }

	if(debug>11) cout << i << "\t" << j << endl;

	while(i>0 || j>0) {
	    if(debug>11) cout << i << "i\t" << j << "j\t" << T[i][j] << "\t" << T_D[i][j] << "\t" << T_I[i][j] << endl;
	    if(T[i][j]==0) {
		tr[pt] = 0;
		i--; j--;
		pt++;
	    }
	    else if(T[i][j] < 0) {
		m = T_D[i][j];
		for(k=1;k<=m;k++) {
			tr[pt] = -1;
			i--;
			pt++;
		}
	    }
	    else if(T[i][j] > 0) {
		m = T_I[i][j];
		for(k=1;k<=m;k++) {
		   	tr[pt] = 1;
			j--;
			pt++;
		}
	    }
	    else {
		warning("not a valid number of elements of T");
		cout << i << "\t" << j << "\t" << T[i][j] << endl;
		exit(0);
	    }
	}

	// determine the reverse path
	for(i=pt-1;i>=0;i--) {
	    r_tr[pt-i] = tr[i];
	    if(debug>11) cout << tr[i] << " ";
	}

	if(debug>11) cout << endl;
	*/

}

