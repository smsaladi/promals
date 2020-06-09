#include "all.h"

static int Debug = 1;

sparseMatrix::sparseMatrix() {

	nrows = ncols = 0;
	nelements = 0;
	// isCrs = 0;
	aux = 0;
}
	
sparseMatrix::sparseMatrix(const sparseMatrix & mat) {

	nrows = mat.nrows; ncols = mat.ncols;
	nelements = mat.nelements;
	rvalue = mat.rvalue;
	rindex = mat.rindex;
	rstart = mat.rstart;
	cvalue = mat.cvalue;
	cindex = mat.cindex;
	cstart = mat.cstart;
	
	aux = NULL;
}

sparseMatrix::sparseMatrix(float **regularMat, int m, int n) {

	aux = 0;
	// if(crs>0) isCrs = 1; else if(crs<0) crs = -1;
	// else { cout << "crs in sparseMatrix must be specified as 1 or -1" << endl; exit(0); }

	regular2Sparse(regularMat, m, n);
	// if(crs==1) regular2Crs(regularMat, m, n);
	// else if(crs==-1) regular2Ccs(regularMat, m, n);

}

sparseMatrix::~sparseMatrix() {

	if(aux) {
		assert(nrows); assert(ncols);
		free_gmatrix<float>(aux, nrows, ncols);
	}
	rvalue.clear();
	rindex.clear();
	rstart.clear();
	cvalue.clear();
	cindex.clear();
	cstart.clear();
}

void sparseMatrix::clear() {

	if(aux) {
		assert(nrows); assert(ncols);
		free_gmatrix<float>(aux, nrows, ncols);
		aux = 0;
	}
	rvalue.clear();
	rindex.clear();
	rstart.clear();
	cvalue.clear();
	cindex.clear();
	cstart.clear();

	nrows = ncols = 0; //isCrs = 0;
}

void sparseMatrix::regular2Sparse(float **regularMat, int m, int n) {

	regular2Crs(regularMat, m, n);
	regular2Ccs(regularMat, m, n);

}

void sparseMatrix::regular2Crs(float **regularMat, int m, int n) {

	int i,j;
	int counts=0;
	float **r = regularMat;

	nrows = m;
	ncols = n;
	nelements = 0;

	rstart.clear();
	rvalue.clear();
	rindex.clear();
	for(i=1;i<=m;i++) {
		rstart.push_back(counts);
		for(j=1;j<=n;j++) {
			if(r[i][j]) {
				rvalue.push_back(r[i][j]);
				rindex.push_back(j);
				nelements++;
				// cout << r[i][j] << "\t" << j << endl;
				counts++;
			}
		}
	}
	rstart.push_back(counts);
	// isCrs = 1; // is a compressed row storage
}


void sparseMatrix::regular2Ccs(float **regularMat, int m, int n) {

	int i,j;
	//int counts=0;
	float **r = regularMat;

	nrows = m;
	ncols = n;
	nelements = 0;

	cstart.clear();
	cindex.clear();
	cvalue.clear();

	for(i=1;i<=n;i++) {
		//cstart.push_back(counts);
		cstart.push_back(nelements);
		for(j=1;j<=m;j++) {
			if(r[j][i]) {
				cvalue.push_back(r[j][i]);
				cindex.push_back(j);
				nelements++;
				//counts++;
			}
		}
	}
	//cstart.push_back(counts);
	cstart.push_back(nelements);
	// isCrs = -1; // is a compressed column storage
}

float **sparseMatrix::sparseCrs2Regular() {

	int i,j;
	float **tmp = gmatrix<float>(nrows, ncols);
	for(i=1;i<=nrows;i++)for(j=1;j<=ncols;j++)tmp[i][j]=0;

	for(i=1;i<=nrows;i++) {
		for(j=rstart[i-1];j<rstart[i];j++) {
			tmp[i][rindex[j]] = rvalue[j];
		}
	}
	return tmp;
}

void sparseMatrix::Crs2Regular() {

	int i,j;

	if(!aux) {
		aux = gmatrix<float>(nrows, ncols);
	}
	else {
		free_gmatrix<float>(aux, nrows, ncols);	
		aux = gmatrix<float>(nrows, ncols);
	}

	for(i=1;i<=nrows;i++) { for(j=1;j<=ncols;j++) { aux[i][j] = 0; } }

	for(i=1;i<=nrows;i++) {
		for(j=rstart[i-1];j<rstart[i];j++) {
			aux[i][rindex[j]] = rvalue[j];
		}
	}
}


void sparseMatrix::Ccs2Regular() {

	int i,j;

	if(!aux) {
		aux = gmatrix<float>(nrows, ncols);
		//cout << "00000000000000000" << endl;
	}
	else {
		for(i=1;i<=nrows;i++)for(j=1;j<=ncols;j++)aux[i][j]=0;
	}

	for(i=1;i<=nrows;i++) { for(j=1;j<=ncols;j++) { aux[i][j] = 0; } }

	for(i=1;i<=ncols;i++) {
		for(j=cstart[i-1];j<cstart[i];j++) {
			aux[cindex[j]][i] = cvalue[j];
		}
	}
}

void sparseMatrix::Crs2Ccs() {

	int i,j;

	// if(isCrs==-1) return; // already in ccs form
	// else isCrs = -1;

	Crs2Regular();
	cvalue.clear();
	cindex.clear();
	cstart.clear();
	regular2Ccs(aux, nrows, ncols);
	free_gmatrix<float>(aux, nrows, ncols);
	aux = NULL;

}

void sparseMatrix::Ccs2Crs() {
	
	int i,j;

	// if(isCrs==1) return; // already in crs form
	// else isCrs = 1;

	 //cout << "++++++++++++" << endl;
	Ccs2Regular();
	 //cout << "++++++++++++" << endl;
	rvalue.clear();
	rindex.clear();
	rstart.clear();
	regular2Crs(aux, nrows, ncols);
	free_gmatrix<float>(aux, nrows, ncols);
	aux = 0;

}

void sparseMatrix::printSparseMatrix(int isCrs) {

	int i;
	if(isCrs==1) {
		printCrs();
	}
	else {
		printCcs();
	}
	
}

void sparseMatrix::printCrs() {

	int i,j;

	cout << "Number of elements in sparse matrix: " << nelements << endl;
	for(i=1;i<=nrows;i++) {
	    cout << "Row number " << i << ":" << endl;
	    for(j=rstart[i-1];j<rstart[i];j++) {
		cout << "\t" << rindex[j] << "\t" << rvalue[j] << endl;
	    }
	}

}

void sparseMatrix::printCcs() {

	int i,j;

	cout << "Number of elements in sparse matrix: " << nelements << endl;
	for(i=1;i<=ncols;i++) {
	    cout << "Column number " << i << ":" << endl;
	    for(j=cstart[i-1];j<cstart[i];j++) {
		cout << "\t" << cindex[j] << "\t" << cvalue[j] << endl;
	    }
	}
}

void sparseMatrix::printAuxMatrix() {
	
	int i,j;
	if(!aux) {
	 	cout << "Auxillary array does not exist" << endl;
		return;
	}
	cout << "Auxillary array: " << endl;
	for(i=1;i<=nrows;i++) {
	    for(j=1;j<=ncols;j++) {
		cout << setw(8) << aux[i][j] << " ";
	    }
	    cout << endl;
	}
}

sparseMatrix * sparseMatrix::transpose() {

	int i,j;

	sparseMatrix * smat = new sparseMatrix();
	
	smat->nrows = ncols;
	smat->ncols = nrows;
	smat->nelements = nelements;
	smat->aux = 0;

	smat->rstart = cstart;
	smat->rindex = cindex;
	smat->rvalue = cvalue;
	smat->cstart = rstart;
	smat->cindex = rindex;
	smat->cvalue = rvalue;

	return smat;

}

void sparseMatrix::multiplyConstant(float c) {

	int i;

	for(i=0;i<rvalue.size();i++) rvalue[i] *= c;
	for(i=0;i<cvalue.size();i++) cvalue[i] *= c;

}

// what this routine does: sum  = sum + a * b
void relaxTwoSparse(sparseMatrix *a, sparseMatrix *b, float **sum) {

	int i, j, k, m, n;
	int r = a->nrows, c = b->ncols;
	int d = a->ncols;

	//vector<int>::iterator p1, p2;
	//vector<float>::iterator v1, v2;
	int p1, p2;

	assert(a->ncols==b->nrows);

	if(Debug>1) cout << "Matrix relaxation r: " << r << "\t c: " << c << endl;
	for(i=1;i<=r;i++) {
 	    for(j=1;j<=c;j++) {
		p1 = a->rstart[i-1];
		p2 = b->cstart[j-1];
		while( (p1<a->rstart[i]) && (p2<b->cstart[j]) ) {
			if(a->rindex[p1]<b->cindex[p2]) {
				p1++;
			}
			else if(a->rindex[p1]>b->cindex[p2]) {
				p2++;
			}
			else {
				sum[i][j] += a->rvalue[p1] * b->cvalue[p2];
				p1++; p2++;
			}
		}
		//cout << setw(5) << sum[i][j] << " ";
			
	    }
	    //cout << endl;
	}
}

float sparseMatrix::getElement(int m, int n) {

	int i;

	assert(m<=nrows); assert(n<=ncols);

	for(i=rstart[m-1];i<rstart[m];i++) {
		if(rindex[i]==n) return rvalue[i];
	}

	return 0.0;
}
