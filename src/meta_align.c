#include "all.h"
#include "btree_template.h"
#include "progressiveAlignHMM.h"
#include "sequences.h"

//static int debug = 1;
vector<sequences> read_seq_from_alignments(char *name); 

int main(int argc, char **argv) {

	int i,j,k;
	int nseqs;

	getParameter(argc, argv, 0);
	//printParameters();

	get_log_robinson_freq();
	get_log_q_blosum62_ratio();
	get_log_q_blosum62();

	if(probconsBLOSUM) useProbconsFrequencies();

	// this function reads a set of alignments into sequences objects
	// the alignment names are in argv[1]
	vector<sequences> seq_vector = read_seq_from_alignments(argv[1]);

	// generate an additional "sequences" object that stores the sequences without any gaps
	sequences tmpseq(seq_vector[0]);
	for(i=1;i<=tmpseq.nseqs;i++) { tmpseq.zapLetters(tmpseq.seq[i]); }

	// calculate distance matrices for each sequences object in seq_vector
	for(i=0;i<seq_vector.size();i++) { seq_vector[i].seqIdentity2DistMat(); }

	// combine the distance matrices into one distance matrix
	nseqs = seq_vector[0].nseqs;
	double **distMat = dmatrix(nseqs, nseqs);
	for(i=1;i<=nseqs;i++) {
		for(j=1;j<=nseqs;j++) {
			for(k=0;k<seq_vector.size();k++) {
				distMat[i][j] += seq_vector[k].distMat[i][j]/seq_vector.size();
			}
			if(debug>100)fprintf(stdout, "%4.3f ", distMat[i][j]);
		}
		if(debug>100)fprintf(stdout, "\n");
	}

	// generate the tree
	btree<tnode> tree;
	tree.UPGMA(distMat, tmpseq.seq, tmpseq.name, tmpseq.nseqs);
	if(debug>1) tree.writeTree("tmptree.tre");

	// generate the consistency matrices
	tree.obtainPreAligned(tree.root);
	tree.alignment_consistency(seq_vector);
	
	// compute consistency alignment
	tree.computeConsistencyAlignment(tree.root);
	char outFileName[200];
	if(argc<=2) {strcpy(outFileName, argv[1]);
			strcat(outFileName, ".meta_align.aln");
	}
	else if(argv[2][0]=='-') {
	    strcpy(outFileName, argv[1]);
	    strcat(outFileName, ".meta_align.aln");
	}
	else {
	    strcpy(outFileName, argv[2]);
	}
	//if(!outFile.empty()) { strcpy(outFileName, outFile.c_str() ); }
	tree.printAlignmentFromAbs(tree.root, outFileName);

	cout << "Output file: " << outFileName << endl;
	cout << "meta_align finished" << endl;

}

// this function reads a set of alignments into sequences objects
// the alignment names are in argv[1]
vector<sequences> read_seq_from_alignments(char *name) {

	vector<sequences> sv;

	int i, j, k;
	char fileName[200];

	ifstream fp(name, ios::in);

	while(fp>>fileName) {
		if(strlen(fileName)==0) { break;}
		sequences sa;
		sa.readFasta(fileName, 0);
		sa.isAlign = 1;
		sv.push_back(sa);
	}

	return sv;
}

