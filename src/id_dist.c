#include "header_cpp.h"
#include "all.h"

int main(int argc, char **argv) {

	int i, j;

	sequences tmpseq(argv[1], 0);

	tmpseq.isAlign = 1;
	
	tmpseq.seqIdentity2DistMat();


	for(i=1;i<=tmpseq.nseqs;i++) {
		for(j=1;j<=tmpseq.nseqs;j++) {
			cout << tmpseq.name[i] << " " << tmpseq.name[j] << " " << tmpseq.distMat[i][j] << endl;
		}
	}

}
