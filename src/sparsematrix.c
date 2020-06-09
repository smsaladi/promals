#include "all.h"

static int Debug = 1;

sparseMatrix::sparseMatrix() {

	nrows = ncols = 0;
	nelements = 0;
	// isCrs = 0;
	aux = 0;
	rvalue = 0;
	rindex = 0;
	rstart = 0;
	cvalue = 0;
	cindex = 0;
	cstart = 0;
}
	
sparseMatrix::sparseMatrix(const sparseMatrix & mat) {

	int i;
	nrows = mat.nrows; ncols = mat.ncols;
	nelements = mat.nelements;

	if(mat.rvalue) {
		rvalue = gvector<float>(nelements);
		for(i=0;i<nelements;i++) rvalue[i] = mat.rvalue[i];
	}
	if(mat.rindex) {
		rindex = ivector(nelements);
		for(i=0;i<nelements;i++) rindex[i] = mat.rindex[i];
	}
	if(mat.rstart) {
		rstart = ivector(nrows);
		for(i=0;i<nrows;i++) rstart[i] = mat.rstart[i];
	}
	if(mat.cvalue) {
		cvalue = gvector<float>(nelements);
		for(i=0;i<nelements;i++) cvalue[i] = mat.cvalue[i];
	}
	if(mat.cindex) {
		cindex = ivector(nelements);
		for(i=0;i<nelements;i++) cindex[i] = mat.cindex[i];
	}
	if(mat.cstart) {
		cstart = ivector(ncols);
		for(i=0;i<ncols;i++) cstart[i] = mat.cstart[i];
	}
	/*
	rvalue = mat.rvalue;
	rindex = mat.rindex;
	rstart = mat.rstart;
	cvalue = mat.cvalue;
	cindex = mat.cindex;
	cstart = mat.cstart;
	*/
	
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
	delete [] rvalue;
	delete [] rindex;
	delete [] rstart;
	delete [] cvalue;
	delete [] cindex;
	delete [] cstart;
	/*
	rvalue.clear();
	rindex.clear();
	rstart.clear();
	cvalue.clear();
	cindex.clear();
	cstart.clear();
	*/
}

void sparseMatrix::clear() {

	if(aux) {
		assert(nrows); assert(ncols);
		free_gmatrix<float>(aux, nrows, ncols);
		aux = 0;
	}
	delete [] rvalue;
	delete [] rindex;
	delete [] rstart;
	delete [] cvalue;
	delete [] cindex;
	delete [] cstart;
	/*
	rvalue.clear();
	rindex.clear();
	rstart.clear();
	cvalue.clear();
	cindex.clear();
	cstart.clear();
	*/

	nrows = ncols = 0; //isCrs = 0;
}

void sparseMatrix::regular2Sparse(float **regularMat, int m, int n) {

	regular2Crs(regularMat, m, n);
	regular2Ccs(regularMat, m, n);

}

void sparseMatrix::regular2Crs(float **regularMat, int m, int n) {

	int i,j;
	//int counts=0;
	float **r = regularMat;

	nrows = m;
	ncols = n;
	nelements = 0;
	int tmp_n = 0;

	for(i=1;i<=m;i++) {
		for(j=1;j<=n;j++) {
			if(r[i][j]) nelements++;
		}
	}

	rstart = ivector(m);
	rvalue = gvector<float>(nelements);
	rindex = ivector(nelements);

	for(i=1;i<=m;i++) {
		rstart[i-1] = tmp_n;
		for(j=1;j<=n;j++) {
			if(r[i][j]) {
				rvalue[tmp_n] = r[i][j];
				rindex[tmp_n] = j;
				tmp_n++;
			}
		}
	}
	rstart[i-1] = tmp_n;

	/*
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
	*/
	// isCrs = 1; // is a compressed row storage
}


void sparseMatrix::regular2Ccs(float **regularMat, int m, int n) {

	int i,j;
	//int counts=0;
	float **r = regularMat;

	nrows = m;
	ncols = n;
	nelements = 0;
	int tmp_n = 0;

	for(i=1;i<=m;i++) {
		for(j=1;j<=n;j++) {
			if(r[i][j]) nelements++;
		}
	}

	cstart = ivector(n);
	cvalue = gvector<float>(nelements);
	cindex = ivector(nelements);

	for(i=1;i<=n;i++) {
		cstart[i-1] = tmp_n; 
		for(j=1;j<=m;j++) {
			if(r[j][i]) {
				cvalue[tmp_n] = r[j][i];
				cindex[tmp_n] = j;
				tmp_n++;
			}
		}
	}
	cstart[i-1] = tmp_n;

	/*
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
	*/
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
	delete [] cvalue;
	delete [] cindex;
	delete [] cstart;
	/*
	cvalue.clear();
	cindex.clear();
	cstart.clear();
	*/
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
	delete [] rvalue;
	delete [] rindex;
	delete [] rstart;
	/*
	rvalue.clear();
	rindex.clear();
	rstart.clear();
	*/
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

	if(rvalue) {
		smat->cvalue = gvector<float>(nelements);
		for(i=0;i<nelements;i++) smat->cvalue[i] = rvalue[i];
	}
	if(rindex) {
		smat->cindex = ivector(nelements);
		for(i=0;i<nelements;i++) smat->cindex[i] = rindex[i];
	}
	if(rstart) {
		smat->cstart = ivector(nrows);
		for(i=0;i<nrows;i++) smat->cstart[i] = rstart[i];
	}
	if(cvalue) {
		smat->rvalue = gvector<float>(nelements);
		for(i=0;i<nelements;i++) smat->rvalue[i] = cvalue[i];
	}
	if(cindex) {
		smat->rindex = ivector(nelements);
		for(i=0;i<nelements;i++) smat->rindex[i] = cindex[i];
	}
	if(cstart) {
		smat->rstart = ivector(ncols);
		for(i=0;i<ncols;i++) smat->rstart[i] = cstart[i];
	}
	/*
	smat->rstart = cstart;
	smat->rindex = cindex;
	smat->rvalue = cvalue;
	smat->cstart = rstart;
	smat->cindex = rindex;
	smat->cvalue = rvalue;
	*/

	return smat;

}

void sparseMatrix::multiplyConstant(float c) {

	int i;
	/*
	for(i=0;i<rvalue.size();i++) rvalue[i] *= c;
	for(i=0;i<cvalue.size();i++) cvalue[i] *= c;
	*/
	if(rvalue) for(i=0;i<nelements;i++) rvalue[i] *= c;
	if(cvalue) for(i=0;i<nelements;i++) cvalue[i] *= c;

}

// what this routine does: sum  = sum + a * b
void relaxTwoSparse(sparseMatrix *a, sparseMatrix *b, float **sum) {

	int i, j, k, m, n;
	int r = a->nrows, c = b->ncols;
	int d = a->ncols;
	float value_diff;

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
			/*
			value_diff = a->rindex[p1] - b->cindex[p2];
			if(value_diff<0) p1++;
			else if(value_diff>0) p2++;
			else {
				p1++; p2++;
				sum[i][j] += a->rvalue[p1] * b->cvalue[p2];
			}
			*/
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
