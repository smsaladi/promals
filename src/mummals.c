#include "all.h"
#include "sequences.h"
#include "btree_template.h"
#include "progressiveAlignHMM.h"

//static int debug = 1;

int main(int argc, char **argv) {

	int i,j,k;

	getParameter(argc, argv, 1);
	//printParameters();

	get_log_robinson_freq();
	get_log_q_blosum62_ratio();
	get_log_q_blosum62();

	if(probconsBLOSUM) useProbconsFrequencies();

	multiple_alignment(argv[1]);

}

void multiple_alignment(char *input_file) {

	sequences tmpseq(input_file, 1);

	sequences tmpseq1(tmpseq);

	if(debug>1) tmpseq1.printSeqs();
	tmpseq1.toDayhoff6();
	//cout << endl << endl;
	if(debug>1) tmpseq1.printSeqs();

	tmpseq1.generateD6t(6);
	tmpseq1.diffCountD2t(tmpseq1.d6t[1], tmpseq1.d6t[2]);


if(debug>1) {
	cout << "Difference in Kmer counts:" << endl;
	for(i=1;i<=tmpseq1.nseqs;i++) {
		for(j=1;j<=tmpseq1.nseqs;j++) {
			fprintf(stdout, "%5d", tmpseq1.diffCountD2t(tmpseq1.d6t[i], tmpseq1.d6t[j]) );
		}
		cout << endl;
	}
	cout << endl;
	cout << "Common counts in Kmer counts:" << endl;
	for(i=1;i<=tmpseq1.nseqs;i++) {
		for(j=1;j<=tmpseq1.nseqs;j++) {
			fprintf(stdout, "%5d", tmpseq1.commonCountD2t(tmpseq1.d6t[i], tmpseq1.d6t[j]) );
		}
		cout << endl;
	}
	cout << endl;
	cout << "Estimated evolutionary distances using Kmer counts:" << endl;
	cout << "  " << tmpseq1.nseqs << endl;
}
	tmpseq1.d6t2DistMat(6);
if(debug>1) {
	for(i=1;i<=tmpseq1.nseqs;i++) {
		//cout << tmpseq1.name[i] << "  ";
		for(j=1;j<=tmpseq1.nseqs;j++) {
			cout << tmpseq1.name[i] << "  " << tmpseq1.name[j] << " " << tmpseq1.distMat[i][j] << endl;
			//fprintf(stdout, "%5f ", tmpseq1.distMat[i][j]) ;
		}
		//cout << endl;
	}
	//tmpseq1.printDistMat();
}


	// test the btree
	btree<tnode> tree;
	tree.UPGMA(tmpseq1.distMat, tmpseq.seq, tmpseq1.name, tmpseq1.nseqs);
	if(debug>1) tree.writeTree("tmptree.tre");

	// progressive alignment
	tree.progressiveAlignHMM(tree.root);

	// get a new sequence object
	sequences tmpseq2( *(tree.root->aln) );
	tmpseq2.get_map();
	if(debug>1) tmpseq2.printSeqs();
	// reestimate the distance matrix
	tmpseq2.seqIdentity2DistMat();
	if(debug>1) tmpseq2.printDistMat();
	// re-align using the new distance matrix
	btree<tnode> tree1;
	tree1.UPGMA(tmpseq2.distMat, tmpseq2.seq, tmpseq2.name, tmpseq2.nseqs);
	if(debug>1) tree1.writeTree("tmptree2.tre");
	if(debug>1) times(&tmsstart);
	//tree1.progressiveAlignHMM(tree1.root);
	
	// Fast stage of alignment
	tree1.progressiveAlignHMM_FastStage(tree1.root, 1-ave_grp_thr);
	if(debug>1) { times(&tmsend); timeDiff(); }
	tree1.obtainPreAligned(tree1.root);
	//cout << "Number of pre-aligned groups: " << tree1.preAligned.size() << endl << endl;
	if(debug>1) for(i=0;i<tree1.preAligned.size();i++) {
		cout <<  i << "\t" << tree1.preAligned[i]->p << "\t" << tree1.preAligned[i]->branchlen << endl;
	}
	if(debug>1) times(&tmsstart);
	tree1.getProfileForPreAligned();
	if(debug>1) {times(&tmsend); timeDiff();}
	if(debug>1) times(&tmsstart);
	switch (useLocal) {
		case 0: tree1.profileConsistency(); break;
		case 1: tree1.profileConsistency_local(); break;
		case 2: tree1.profileConsistency_glocal(weightG); break;
		case 3: {
			if(strlen(parameter_file)==0) {
				cout << "must input a parameter file for multim option" << endl;
				exit(0);
			}
			if(debug>1) cout << "Option 3: use multim" << endl;
			hmm_parameters params(solv,ss,unaligned);
			params.read_parameters(parameter_file);
			tree1.profileConsistency_multim(&params); 
			break;
			}
		case 4: {
			if( (strlen(parameter_file1)==0) || (strlen(parameter_file2)==0) ) {
				cout << "must input a parameter file for multim option" << endl;
				exit(0);
			}
			if(debug>1) cout << "Option 3: use multim" << endl;
			hmm_parameters params1(solv,ss,unaligned);
			params1.read_parameters(parameter_file1);
			hmm_parameters params2(solv,ss,unaligned);
			params2.read_parameters(parameter_file2);
			tree1.profileConsistency_multim2(&params1, &params2, &tmpseq2, 0.2); 
			break;
			}
		default: //tree1.profileConsistency(); 
			{
			if(strlen(parameter_file)==0) {
                                cout << "must input a parameter file for multim option" << endl;
                                exit(0);
                        }
                        if(debug>1) cout << "Option 3: use multim" << endl;
                        hmm_parameters params(solv,ss,unaligned);
                        params.read_parameters(parameter_file);
                        tree1.profileConsistency_multim(&params);
                        break;
			}
	}
	if(debug>1) { times(&tmsend); timeDiff();}
	//tree1.root->aln->printali("1csp.hmm2.aln", 50);

	if(debug>1) times(&tmsstart);
	tree1.computeConsistencyAlignment(tree1.root);
	//cout << "end here" << endl;
	if(debug>1) {times(&tmsend); timeDiff();}
	if(debug>1) times(&tmsstart);
	//cout << "============" << endl;
	char outFileName[200];
	if(argv[2][0]=='-') {
	    strcpy(outFileName, argv[1]);
	    int tmpLen = strlen(outFileName);
	    if( (outFileName[tmpLen-1]=='a')&&(outFileName[tmpLen-2]=='f')&&(outFileName[tmpLen-3]=='.')){
		outFileName[tmpLen-3] = '\0';
	    }
	    strcat(outFileName, ".mummals.aln");
	}
	else {
	    strcpy(outFileName, argv[2]);
	}
	if(!outFile.empty()) {
		strcpy(outFileName, outFile.c_str() );
	}
	cout << "  output file Name: " << outFileName << endl << endl;
	tree1.printAlignmentFromAbs(tree1.root, outFileName);
	cout << "  program finished"  << endl << endl;
	if(debug>1) {times(&tmsend); timeDiff();}
	
}
