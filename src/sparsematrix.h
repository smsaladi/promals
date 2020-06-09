#ifndef _sparseMatrix__
#define _sparseMatrix__

// implement crs and ccs in the same class
class sparseMatrix {

    public:
	sparseMatrix();
	sparseMatrix(const sparseMatrix & mat);
	sparseMatrix(float **regularMat, int m, int n);
	~sparseMatrix();
	void clear();

	int nrows, ncols;
	int nelements;
	// int isCrs; // compressed row storage
	float *rvalue;
	int *rindex;
	int *rstart;
	float *cvalue;
	int *cindex;
	int *cstart;
	
	/*
	vector<float> rvalue; // compressed row storage
	vector<int> rindex;
	vector<int> rstart;
	vector<float> cvalue; // compressed column storage
	vector<int> cindex;
	vector<int> cstart;
	*/
	
	float **aux;

	void regular2Sparse(float **regularMat, int m, int n);
	void regular2Crs(float **regularMat, int m, int n);
	void regular2Ccs(float **regularMat, int m, int n);
	float **sparseCrs2Regular();
	void Crs2Regular();
	void Ccs2Regular();
	void Crs2Ccs();
	void Ccs2Crs();

	void printCrs();
	void printCcs();
	void printSparseMatrix(int isCrs);
	void printAuxMatrix();

	sparseMatrix multiply(sparseMatrix &m1, sparseMatrix &m2);
	sparseMatrix addition(sparseMatrix &m1, sparseMatrix &m2);
	sparseMatrix * transpose();
	void multiplyConstant(float c);

	float getElement(int m, int n);
};

// for probablistic consistency measure
extern void relaxTwoSparse(sparseMatrix *a, sparseMatrix *b, float **sum);


#endif


