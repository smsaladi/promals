#ifndef __progressiveAlignHmm_
#define __progressiveAlignHmm_

#include "util.h"
#include "amino.h"
#include "btree_template.h"
#include "hmm_profpair.h"
#include "hmm_profpair1.h"
#include "tnode.h"
#include "mm.h"
#include "sequences.h"
//#include "ScoreType.h"
#include "param.h"
#include "hmm_multim.h"
#include <vector>
#include "multiple.h"
#include "gap_refining.h"
#include "refine.h"
#include "profilehmm.h"
#include "hmm_local.h"
#include "refinegap.h"
#include "time.h"

#define MIN(x, y)       (((x) < (y)) ? (x) : (y)) 

static int Debug = 1;

char ssint2ss(int i);
extern FILE *logfp;

template <typename TNODE>
void btree<TNODE>::progressiveAlignHMM(TNODE *n) {
	
	int i,j;
	subalign *a, *b;
	
	if(n->aligned) return;
	if(!n->childL->aligned) {
		progressiveAlignHMM(n->childL);
	}
	if(!n->childR->aligned) {
		progressiveAlignHMM(n->childR);
	}

	//cout << "======" << endl;
	a = n->childL->aln; b = n->childR->aln;
	hmm_profpair1 profHMM(a, b);

	a->gap_threshold = b->gap_threshold = 1;
	a->beta=0; b->beta=0;  // make the pseudocount weak

	a->profile();b->profile();
	//cout << "++++++++++++++++" << endl;
	
	profHMM.viterbi();
	//profHMM.forward();
	//profHMM.backward();
	subalign *newalign = profHMM.productViterbiAlign();
	n->aln = newalign;
	n->aligned = 1;
	if(Debug>1) n->aln->printali(70);
	if(Debug>1) cout << "------------" << endl;
}
			
template <typename TNODE>
void btree<TNODE>::progressiveAlignHMM_FastStage(TNODE *n, float distCutoff) {
	
	int i,j;
	subalign *a, *b;
	
	if(n->aligned) return;
	if(!n->childL->aligned) {
		progressiveAlignHMM_FastStage(n->childL, distCutoff);
	}
	if(!n->childR->aligned) {
		progressiveAlignHMM_FastStage(n->childR, distCutoff);
	}

	//if( (!n->childL->aligned) || (!n->childR->aligned) ) { return; }

	// determine the distance from the TNODE to the leafs
	float dist=0;
	TNODE *k = n;
	while(k->childL!=0) {
		dist += k->childL->branchlen;
		k = k->childL;
	}
	//cout << "dist: " << dist << endl;
	// Debug here
	// NEW debug
	if(Debug>1) {
		if (dist>distCutoff) {
			if(n->childL->aligned) {
				cout << "Just aligned size: " << n->childL->aln->alilen << " " << n->childL->aln->nal << " " << n->childL->aln->alilen * n->childL->aln->nal << endl;
			}
			if(n->childR->aligned) {
				cout << "Just aligned size: " << n->childR->aln->alilen << " " << n->childR->aln->nal << " " << n->childR->aln->alilen * n->childR->aln->nal << endl;
			}
		}
	}
	if( (!n->childL->aligned) || (!n->childR->aligned) ) { return; }
	if(dist > distCutoff) return;

	//cout << "======" << endl;
	a = n->childL->aln; b = n->childR->aln;
	hmm_profpair1 profHMM(a, b);

	a->gap_threshold = b->gap_threshold = 1;
	a->beta=0; b->beta=0;  // make the pseudocount weak

	a->profile();b->profile();
	//cout << "++++++++++++++++" << endl;
	
	profHMM.viterbi();
	//profHMM.forward();
	//profHMM.backward();
	subalign *newalign = profHMM.productViterbiAlign();
	if(Debug>1) newalign->printali(70);
	if(Debug>1) cout << "------------" << endl;
	gap_refine *gr = new gap_refine();
	char **tmpseq = newalign->aseq;
	newalign->aseq = gr->batch_refine(newalign, 2);
	free_cmatrix(tmpseq, newalign->nal, newalign->alilen);
	delete_complete_gap_positions(newalign);
	newalign->alilen = strlen(newalign->aseq[0]);
	newalign->convertAseq2Alignment();
	n->aln = newalign;
	n->aligned = 1;
	delete gr;
	// delete the alignment in child nodes
	delete n->childL->aln;
	delete n->childR->aln;
	// NEW debug
	if(Debug>1) {
		cout << "The alignment size: " << n->aln->alilen << " " <<  n->aln->nal << " " << n->aln->alilen * n->aln->nal << endl;
	}
	if(Debug>1) n->aln->printali(70);
	if(Debug>1) cout << "------------" << endl;
	/*
	tnode *T = n->childL;
	double len = 0;
	while(T) { 
		len += T->branchlen;
		T = T->childL;
	}
	cout << "distance to the leaf node: " << len << endl;
	*/
}

// return an array of leaf tnode indexes, which will be used to extract sequences and their names
template <typename TNODE>
vector<int> btree<TNODE>::progressiveAlignHMM_FastStage_mafft(TNODE *n, float distCutoff) {
	
	int i,j;
	//subalign *a, *b;
        //cout << "This place" << endl;
        //
        //

        // 1. if distCutoff is negative
        if(distCutoff < 0) {
                if(n->childL==NULL) {
                        vector<int> index;
                        //cout << "n-n: " << n->n << endl;
                        index.push_back(n->n);
                        return index;
                }
                vector<int> index1 = progressiveAlignHMM_FastStage_mafft(n->childL, distCutoff);
                vector<int> index2 = progressiveAlignHMM_FastStage_mafft(n->childR, distCutoff);
                vector<int> index = index1;
                for(i=0;i<index2.size();i++) {
                        index.push_back(index2[i]);
                }
                return index;
        }

	// 2. determine the distance from the TNODE to the leafs
	float dist=0;
	TNODE *k = n;
	while(k->childL!=0) {
		dist += k->childL->branchlen;
		k = k->childL;
	}

        //cout << distCutoff << " " << dist << endl;

        // 3. if distance to the leafs is larger than distCutoff, go the the children
        if(dist>distCutoff) {
                progressiveAlignHMM_FastStage_mafft(n->childL, distCutoff);
                progressiveAlignHMM_FastStage_mafft(n->childR, distCutoff);
                vector<int> vnull;
                return vnull;
        }

        // 4. if distance to the leafs is smaller or equal to distCutoff, check the distance of its parent 
        // 4.1 a leaf node
        if(n->childL==NULL) {
                vector<int> index;
                //cout << "n-n: " << n->n << endl;
                index.push_back(n->n);
                return index;
        }
        // 4.2 an internal node
        float parent_dist;
        if(n->rootFlag) parent_dist = dist + 1.0;
        else parent_dist = dist + n->branchlen;
        vector<int> index1 = progressiveAlignHMM_FastStage_mafft(n->childL, distCutoff);
        vector<int> index2 = progressiveAlignHMM_FastStage_mafft(n->childR, distCutoff);
        vector<int> index = index1;
        for(i=0;i<index2.size();i++) {
                index.push_back(index2[i]);
        }
        // if the distance of its parent is smaller than distCutoff,  just return the vector of indexes
        // also set this node status to be 'aligned'
        if(parent_dist <= distCutoff) {
                n->aligned = 1;
                return index; 
        }

        // 5. distance is smaller or equal to distCutoff, while distance of parent is larger than distCutoff
        //    run mafft on these cases
        // 5.1. set up a directory
        char tmpdir[200];
        char command[200];
        long int a1 = time(NULL);
        srand(a1);
        int myrand = rand();
        sprintf(tmpdir, "/tmp/%d_%d", getpid(), myrand);
        sprintf(command, "mkdir %s", tmpdir);
        //cout << command << endl;
        system(command);

        // 5.2 write the sequences to a fasta file in the tmp directory
        char tmpfa[300];
        sprintf(tmpfa, "%s/tmp.fa", tmpdir);
        ofstream ofp(tmpfa, ios::out);
        for(i=0;i<index.size();i++) {
                //cout << i << " " << index[i] << endl;
                //cout << index[i] << endl; cout << v[index[i]]->aligned << endl;
                ofp << ">" << v[index[i]]->aln->aname[0] << endl;
                ofp << v[index[i]]->aln->aseq[0] << endl;
        }
        ofp.close();
        // debug
        //cout << "Debug of index" << endl;
        //cout << "index size: " << index.size() << endl;

        // 5.3 run mafft on the sequence
        //sprintf(command, "%s --maxiterate 1000 --localpair %s/tmp.fa 1>%s/tmp.aln.fa 2>%s/tmp.err", mafft, tmpdir, tmpdir, tmpdir);
        sprintf(command, "%s --auto %s/tmp.fa 1>%s/tmp.aln.fa 2>%s/tmp.err", mafft, tmpdir, tmpdir, tmpdir);
        system(command);
        char outalnfile[200];
        sprintf(outalnfile, "%s/tmp.aln.fa", tmpdir);
        //cout << "This place" << endl;
        subalign *a = new subalign(outalnfile, "fasta", 25);
        n->aln = a;
        n->aligned = 1;
        //cout << "Debug here" << endl;
        a->printali(80);
        //cout << "This place" << endl;
        //exit(0);
        // 5.4. clear up the tmp directory
        sprintf(command, "rm -rf %s", tmpdir);
        system(command);

        // 5.5. return a null vector
        vector<int> vnull;
        return vnull;

}
			
template <typename TNODE>
void  btree<TNODE>::obtainPreAligned(TNODE *n) {

	int i;
	if(n==0) return;

	/*if(n->aligned && (!n->parent->aligned) ) n->aln->printali(70);
	if(!n->rootFlag) cout << "aligned: " << n->aligned << " parent_aligned: " << n->parent->aligned << endl; 
	if(n->aligned) n->aln->printali(70);
	tnode *T = n->childL;
	double len = 0;
	while(T) { 
		len += T->branchlen;
		T = T->childL;
	}
	cout << "distance to the leaf node: " << len << endl;
	*/
	// if root is already aligned
	//if( (n->aligned) && (n->rootFlag) ) {return;}
	if( (n->aligned) && (n->rootFlag) ) {
		TNODE *a = n;
		//cout << "Root here" << endl;
		preAligned.push_back(a);
		n->p = preAligned.size()-1;
		// get the abstractSeq
		int *tmpList = new int [n->aln->alilen+1];
		for(i=1;i<=n->aln->alilen;i++) tmpList[i] = i;
		tmpList[0] = n->p;
		//cout << "tmpList 0: " <<  tmpList[0] << endl;
		n->abstractSeq.push_back(0);
		n->abstractSeq.push_back(tmpList);
		n->absSeqnum = 1;
		n->absSeqlength = n->aln->alilen;
		return;
	}
	//preAligned.clear();
	if(n->aligned && (!n->rootFlag) && (!n->parent->aligned) ) {
		//n->aln->printali(80);
		TNODE *a = n;
		preAligned.push_back(a);
		n->p = preAligned.size()-1;
		//cout << n->aln->nal << endl;
		// get the abstractSeq
		int *tmpList = new int [n->aln->alilen+1];
		for(i=1;i<=n->aln->alilen;i++) tmpList[i] = i;
		tmpList[0] = n->p;
		//cout << "tmpList 0: " <<  tmpList[0] << endl;
		n->abstractSeq.push_back(0);
		n->abstractSeq.push_back(tmpList);
		n->absSeqnum = 1;
		n->absSeqlength = n->aln->alilen;
		//cout <<"=======================" << endl;
		//n->aln->printali(80);
		//cout <<"+++++++++++++++++++++++" << endl;
		//if(!n->childL) {cout << "afasdfafd\n" << endl; }
	}
	obtainPreAligned(n->childL);
	obtainPreAligned(n->childR);
}

template <typename TNODE>
void btree<TNODE>::getProfileForPreAligned() {

	int i;

	for(i=0;i<preAligned.size();i++) {
		preAligned[i]->aln->gap_threshold=1; // this means all positions are included regardless of gap content
		preAligned[i]->aln->beta = 0; // this means do not mix with pseudocount frequencies; that is, only raw frequencies are used for the profile
		preAligned[i]->aln->profile(); // calculate the profile
	}
}

template <typename TNODE>
void btree<TNODE>::profileConsistency() {

	int i,j=0,k,l,m,n;
	float **realmat;
	int lenx, leny;
	int nonZeroCounts = 0;
	
	smat = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);
	smat1 = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);
	//smat = new sparseMatrix **[preAligned.size()];
	//for(i=0;i<preAligned.size();i++) {
	//	smat[i] = new sparseMatrix *[preAligned.size()];
	//	for(j=0;j<preAligned.size();j++) {
	//		smat[i][j] = new  sparseMatrix();
	//		cout << "i:" << i << "\t" << smat[i][j]->nrows << endl;
	//	}
	//}
	//exit(0);
	//sparseMatrix testSMAT;
	//cout << "PreAligned size: " << preAligned.size()<<endl;

	// pairwise profile alignments; generate probability matrices (sparse)
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
			//cout << "i: " << i << "\t" << "j: " << "\t" << j << endl;
			nonZeroCounts = 0;
			hmm_profpair1 *hmmProfPair = new hmm_profpair1(preAligned[i]->aln, preAligned[j]->aln);
		        hmmProfPair->viterbi();
			if(Debug>1) cout << "Pairwise viterbi alignment: " << i << "\t" << j << endl;
		        subalign *newalign = hmmProfPair->productViterbiAlign();
		        if(Debug>1) newalign->printali(70);
			delete newalign;
		        if(Debug>1) cout << "------------" << endl;

			hmmProfPair->forward();
			hmmProfPair->backward();
			realmat = hmmProfPair->probMat;
			lenx = hmmProfPair->lenx; leny = hmmProfPair->leny;
		        for(k=1;k<=lenx;k++) {
               		        for(l=1;l<=leny;l++) {
                       		    if(realmat[k][l]<minProb) realmat[k][l] = 0;
				    else nonZeroCounts++;
                       		    //cout << realmat[k][l] << " ";
                		} //cout << endl;
        		}
			if(debug>1) cout << "None-zero counts: " << i << "\t" << j << "\t" << nonZeroCounts << endl;
			//cout << "Here " << endl;
			smat[i][j] = new sparseMatrix(realmat, lenx, leny);
			if(Debug>1) smat[i][j]->printCrs();
			//smat[i][j] = new sparseMatrix;
			//smat[i][j]->nrows=0; smat[i][j]->ncols=0;
			//cout << smat[i][j]->nrows << "\t" << smat[i][j]->ncols<< endl;
			//smat[i][j]->regular2Sparse(realmat, lenx, leny);
			//testSMAT.regular2Sparse(realmat, lenx, leny);
			//cout << "Here " << endl;
			smat[j][i] = smat[i][j]->transpose();
			//cout << "Here " << endl;
			delete hmmProfPair;
		}
	}
	if(debug>1) {cout << "Time after the initial generation of consistency matrices: " << endl;
	times(&tmsend); timeDiff();
	times(&tmsstart);
	}
	
	for(i=0;i<relax_number;i++) {
	    relaxConsistMatrix();
	}
}

template <typename TNODE>
void btree<TNODE>::profileConsistency_multim(hmm_parameters *params) {


	int i,j=0,k,l,m,n;
	float **realmat, **realmat1;
	int lenx, leny, lenx1, leny1;
	int nonZeroCounts = 0;
	
	smat = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);
	smat1 = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);


	// pairwise profile alignments; generate probability matrices (sparse)
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
			//cout << "i: " << i << "\t" << "j: " << "\t" << j << endl;
			nonZeroCounts = 0;
			hmm_multim *hmmProfPair = new hmm_multim(preAligned[i]->aln, preAligned[j]->aln);
			hmmProfPair->set_parameters(params);
			//hmm_profpair1 *hmmProfPair = new hmm_profpair1(preAligned[i]->aln, preAligned[j]->aln);
		        //hmmProfPair->viterbi();
			//cout << "Pairwise viterbi alignment: " << i << "\t" << j << endl;
		        //subalign *newalign = hmmProfPair->productViterbiAlign();
		        //newalign->printali(70);
			//delete newalign;
		        //cout << "------------" << endl;

			hmmProfPair->forward();
			hmmProfPair->backward();
			realmat = hmmProfPair->probMat;
			lenx = hmmProfPair->lenx; leny = hmmProfPair->leny;
			if(reverse_align_order) {
				hmm_multim *hmmProfPair1 = new hmm_multim(preAligned[j]->aln, preAligned[i]->aln);
				hmmProfPair1->set_parameters(params);
				hmmProfPair1->forward();
				hmmProfPair1->backward();
				realmat1 = hmmProfPair1->probMat;
				lenx1 = hmmProfPair1->lenx; leny1 = hmmProfPair1->leny;
				//cout << lenx1 << "  " << leny1 << endl;
				//cout << lenx << "  " << leny << endl;
			        for(k=1;k<=lenx1;k++) {
       	        		        for(l=1;l<=leny1;l++) {
					    realmat[l][k] = (realmat1[k][l]+realmat[l][k])/2;
					}
                		}
				delete hmmProfPair1;
        		}
				
		        for(k=1;k<=lenx;k++) {
               		        for(l=1;l<=leny;l++) {
                       		    if(realmat[k][l]<minProb) realmat[k][l] = 0;
				    else nonZeroCounts++;
                       		    //cout << realmat[k][l] << " ";
                		} //cout << endl;
        		}
			if(debug>1) cout << "None-zero counts: " << i << "\t" << j << "\t" << nonZeroCounts << endl;
			//cout << "Here " << endl;
			smat[i][j] = new sparseMatrix(realmat, lenx, leny);
			if(Debug>1) smat[i][j]->printCrs();
			//smat[i][j] = new sparseMatrix;
			//smat[i][j]->nrows=0; smat[i][j]->ncols=0;
			//cout << smat[i][j]->nrows << "\t" << smat[i][j]->ncols<< endl;
			//smat[i][j]->regular2Sparse(realmat, lenx, leny);
			//testSMAT.regular2Sparse(realmat, lenx, leny);
			//cout << "Here " << endl;
			smat[j][i] = smat[i][j]->transpose();
			//cout << "Here " << endl;
			delete hmmProfPair;
			

		}
	}
	if(debug>1) {
	cout << "Time after the initial generation of consistency matrices: " << endl;
	times(&tmsend); timeDiff();
	times(&tmsstart);
	}
	
	//relaxConsistMatrix();
	//relaxConsistMatrix();
	for(i=0;i<relax_number;i++) {
	    relaxConsistMatrix();
	}
	/*
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
		    cout << "####### " << i << "   " << j << endl;
		    cout << preAligned[i]->aln->aname[0] << " " << preAligned[j]->aln->aname[0] << endl;
		    smat[i][j]->printSparseMatrix(1);
		}
	}
	*/
}

template <typename TNODE>
void btree<TNODE>::profileConsistency_multim(hmm_parameters *params, double **dist_matrix, double max_dist_cutoff) {


	int i,j=0,k,l,m,n;
	float **realmat, **realmat1;
	int lenx, leny, lenx1, leny1;
	int nonZeroCounts = 0;
	
	smat = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);
	smat1 = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);

	// pairwise profile alignments; generate probability matrices (sparse)
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
			//cout << "i: " << i << "\t" << "j: " << "\t" << j << endl;
			nonZeroCounts = 0;
			hmm_multim *hmmProfPair = new hmm_multim(preAligned[i]->aln, preAligned[j]->aln);
			hmmProfPair->set_parameters(params);
			//hmm_profpair1 *hmmProfPair = new hmm_profpair1(preAligned[i]->aln, preAligned[j]->aln);
		        //hmmProfPair->viterbi();
			//cout << "Pairwise viterbi alignment: " << i << "\t" << j << endl;
		        //subalign *newalign = hmmProfPair->productViterbiAlign();
		        //newalign->printali(60);
			//delete newalign;
		        //cout << "------------" << endl;

			hmmProfPair->forward();
			hmmProfPair->backward();
			realmat = hmmProfPair->probMat;
			lenx = hmmProfPair->lenx; leny = hmmProfPair->leny;
			if(reverse_align_order) {
				hmm_multim *hmmProfPair1 = new hmm_multim(preAligned[j]->aln, preAligned[i]->aln);
				hmmProfPair1->set_parameters(params);
				hmmProfPair1->forward();
				hmmProfPair1->backward();
				realmat1 = hmmProfPair1->probMat;
				lenx1 = hmmProfPair1->lenx; leny1 = hmmProfPair1->leny;
				//cout << lenx1 << "  " << leny1 << endl;
				//cout << lenx << "  " << leny << endl;
			        for(k=1;k<=lenx1;k++) {
       	        		        for(l=1;l<=leny1;l++) {
					    realmat[l][k] = (realmat1[k][l]+realmat[l][k])/2;
					}
                		}
				delete hmmProfPair1;
        		}
				
		        for(k=1;k<=lenx;k++) {
               		        for(l=1;l<=leny;l++) {
                       		    if(realmat[k][l]<minProb) realmat[k][l] = 0;
				    else nonZeroCounts++;
                       		    //cout << realmat[k][l] << " ";
                		} //cout << endl;
        		}
			if(debug>1) cout << "None-zero counts: " << i << "\t" << j << "\t" << nonZeroCounts << endl;
			//cout << "Here " << endl;
			smat[i][j] = new sparseMatrix(realmat, lenx, leny);
			//cout << "Here " << endl;
			if(Debug>1) smat[i][j]->printCrs();
			//smat[i][j] = new sparseMatrix;
			//smat[i][j]->nrows=0; smat[i][j]->ncols=0;
			//cout << smat[i][j]->nrows << "\t" << smat[i][j]->ncols<< endl;
			//smat[i][j]->regular2Sparse(realmat, lenx, leny);
			//testSMAT.regular2Sparse(realmat, lenx, leny);
			//cout << "Here " << endl;
			smat[j][i] = smat[i][j]->transpose();
			//cout << "Here " << endl;
			delete hmmProfPair;
			

		}
	}
	if(debug>1) {
	cout << "Time after the initial generation of consistency matrices: " << endl;
	times(&tmsend); timeDiff();
	times(&tmsstart);
	}
	
	//relaxConsistMatrix();
	//relaxConsistMatrix();
	for(i=0;i<relax_number;i++) {
	    //cout << " Before relax" << endl;
	    relaxConsistMatrix(dist_matrix, max_dist_cutoff);
	    //cout << " After relax" << endl;
	}
	/*
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
		    cout << "####### " << i << "   " << j << endl;
		    cout << preAligned[i]->aln->aname[0] << " " << preAligned[j]->aln->aname[0] << endl;
		    smat[i][j]->printSparseMatrix(1);
		}
	}
	*/
}

template <typename TNODE>
void btree<TNODE>::profileConsistency_psipred(hmm_psipred_parameters *params, double **dist_matrix, double max_dist_cutoff, int use_homologs) {

	int i,j=0,k,l,m,n;
	float **realmat, **realmat1;
	int lenx, leny, lenx1, leny1;
	int nonZeroCounts = 0;
	
	smat = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);
	smat1 = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);

	// pairwise profile alignments; generate probability matrices (sparse)
	for(i=0;i<(int)preAligned.size();i++) {
		fprintf(logfp, "*");
		fflush(logfp);
		for(j=i+1;j<(int)preAligned.size();j++) {
			//cout << "i: " << i << "\t" << "j: " << "\t" << j << endl;
			nonZeroCounts = 0;
			hmm_psipred *hmmProfPair = new hmm_psipred(preAligned[i]->aux_align, preAligned[j]->aux_align);
			//preAligned[i]->aln->printali(60);
			//preAligned[j]->aln->printali(60);
			hmmProfPair->set_parameters(params);
			//hmm_profpair1 *hmmProfPair = new hmm_profpair1(preAligned[i]->aln, preAligned[j]->aln);
		        //hmmProfPair->viterbi();
			//cout << "Pairwise viterbi alignment: " << i << "\t" << j << endl;
		        //subalign *newalign = hmmProfPair->productViterbiAlign();
		        //newalign->printali(60);
			//delete newalign;

			//hmmProfPair->get_scores(ss_w);
			if(adjust_weight) {
				int min_nal = hmmProfPair->x->nal;
				if(min_nal>hmmProfPair->y->nal) min_nal = hmmProfPair->y->nal;
				if(min_nal>100) {
					ss_w = 0.25; score_w = 0.8;
				}
				else if(min_nal>5) {
					ss_w = 0.2; score_w = 0.5;
				}
				else { ss_w = 0.35; score_w = 0.5; }
				hmmProfPair->get_scores(ss_w, score_w);
			}
			else {hmmProfPair->get_scores(ss_w, score_w); }
			
			//hmmProfPair->forward1();
			//hmmProfPair->backward1();
                        hmmProfPair->forward_no_end_penalty();
                        hmmProfPair->backward_no_end_penalty();
			realmat = hmmProfPair->probMat;
			lenx = hmmProfPair->lenx; leny = hmmProfPair->leny;
			if(reverse_align_order) {
				hmm_psipred *hmmProfPair1 = new hmm_psipred(preAligned[j]->aln, preAligned[i]->aln);
				hmmProfPair1->set_parameters(params);
				hmmProfPair1->forward1();
				hmmProfPair1->backward1();
				realmat1 = hmmProfPair1->probMat;
				lenx1 = hmmProfPair1->lenx; leny1 = hmmProfPair1->leny;
				//cout << lenx1 << "  " << leny1 << endl;
				//cout << lenx << "  " << leny << endl;
			        for(k=1;k<=lenx1;k++) {
       	        		        for(l=1;l<=leny1;l++) {
					    realmat[l][k] = (realmat1[k][l]+realmat[l][k])/2;
					}
                		}
				delete hmmProfPair1;
        		}
				
		        for(k=1;k<=lenx;k++) {
               		        for(l=1;l<=leny;l++) {
                       		    if(realmat[k][l]<minProb) realmat[k][l] = 0;
				    else nonZeroCounts++;
                       		    //cout << realmat[k][l] << " ";
                		} //cout << endl;
        		}
			//NEW debug
			if(debug>1) cout << "None-zero counts: " << i << "\t" << j << "\t" << nonZeroCounts << endl;
			//cout << "Here " << endl;
			smat[i][j] = new sparseMatrix(realmat, lenx, leny);
		        //cout << "####### " << i << "   " << j << endl;
		        //cout << preAligned[i]->aln->aname[0] << " " << preAligned[j]->aln->aname[0] << endl;
			//smat[i][j]->printCrs(preAligned[i]->aln->aseq[0], preAligned[j]->aln->aseq[0]);
		        //smat[i][j]->printSparseMatrix(1);
			//cout << "Here " << endl;
			if(Debug>1) smat[i][j]->printCrs();
			//smat[i][j] = new sparseMatrix;
			//smat[i][j]->nrows=0; smat[i][j]->ncols=0;
			//cout << smat[i][j]->nrows << "\t" << smat[i][j]->ncols<< endl;
			//smat[i][j]->regular2Sparse(realmat, lenx, leny);
			//testSMAT.regular2Sparse(realmat, lenx, leny);
			//cout << "Here " << endl;
			smat[j][i] = smat[i][j]->transpose();
			//cout << "Here " << endl;
			delete hmmProfPair;
			//cout << "Here " << endl;
		}
	}
	if(debug>1) {
	cout << "Time after the initial generation of consistency matrices: " << endl;
	times(&tmsend); timeDiff();
	times(&tmsstart);
	}
	
	//relaxConsistMatrix();
	//relaxConsistMatrix();
	fprintf(logfp, "\n");
	fflush(logfp);
	/*
	for(i=0;i<relax_number;i++) {
	    //cout << " Before relax" << endl;
	    relaxConsistMatrix(dist_matrix, max_dist_cutoff);
	    fprintf(logfp, "        finished relaxing consistency matrix - round %d of %d\n", i+1, relax_number);
	    fflush(logfp);
	    //cout << " After relax" << endl;
	}
	*/
	//relaxConsistMatrix(dist_matrix, max_dist_cutoff, minProb);
	//relaxConsistMatrix(dist_matrix, max_dist_cutoff, minProb*0.1);
	/*
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
		    cout << "####### " << i << "   " << j << endl;
		    cout << preAligned[i]->aln->aname[0] << " " << preAligned[j]->aln->aname[0] << endl;
		    smat[i][j]->printSparseMatrix(1);
		}
	}
	*/
}

template <typename TNODE>
void btree<TNODE>::profileConsistency_psipred_sum_of_pairs(hmm_psipred_parameters *params, double **dist_matrix, double max_dist_cutoff, int use_homologs) {

	int i,j=0,k,l,m,n;
	float **realmat, **realmat1;
	int lenx, leny, lenx1, leny1;
	int nonZeroCounts = 0;
	
	smat = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);
	smat1 = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);

	// pairwise profile alignments; generate probability matrices (sparse)
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
			//cout << "i: " << i << "\t" << "j: " << "\t" << j << endl;
			nonZeroCounts = 0;
			hmm_psipred *hmmProfPair = new hmm_psipred(preAligned[i]->aux_align, preAligned[j]->aux_align);
			//preAligned[i]->aln->printali(60);
			//preAligned[j]->aln->printali(60);
			hmmProfPair->set_parameters(params);
			//hmm_profpair1 *hmmProfPair = new hmm_profpair1(preAligned[i]->aln, preAligned[j]->aln);
		        //hmmProfPair->viterbi();
			//cout << "Pairwise viterbi alignment: " << i << "\t" << j << endl;
		        //subalign *newalign = hmmProfPair->productViterbiAlign();
		        //newalign->printali(60);
			//delete newalign;

			hmmProfPair->get_scores_sum_of_pairs(ss_w);
			//hmmProfPair->get_scores(ss_w);
			
			hmmProfPair->forward1();
			hmmProfPair->backward1();
			realmat = hmmProfPair->probMat;
			lenx = hmmProfPair->lenx; leny = hmmProfPair->leny;
			if(reverse_align_order) {
				hmm_psipred *hmmProfPair1 = new hmm_psipred(preAligned[j]->aln, preAligned[i]->aln);
				hmmProfPair1->set_parameters(params);
				hmmProfPair1->forward1();
				hmmProfPair1->backward1();
				realmat1 = hmmProfPair1->probMat;
				lenx1 = hmmProfPair1->lenx; leny1 = hmmProfPair1->leny;
				//cout << lenx1 << "  " << leny1 << endl;
				//cout << lenx << "  " << leny << endl;
			        for(k=1;k<=lenx1;k++) {
       	        		        for(l=1;l<=leny1;l++) {
					    realmat[l][k] = (realmat1[k][l]+realmat[l][k])/2;
					}
                		}
				delete hmmProfPair1;
        		}
				
		        for(k=1;k<=lenx;k++) {
               		        for(l=1;l<=leny;l++) {
                       		    if(realmat[k][l]<minProb) realmat[k][l] = 0;
				    else nonZeroCounts++;
                       		    //cout << realmat[k][l] << " ";
                		} //cout << endl;
        		}
			if(debug>1) cout << "None-zero counts: " << i << "\t" << j << "\t" << nonZeroCounts << endl;
			//cout << "Here " << endl;
			smat[i][j] = new sparseMatrix(realmat, lenx, leny);
			//cout << "Here " << endl;
			if(Debug>1) smat[i][j]->printCrs();
			//smat[i][j] = new sparseMatrix;
			//smat[i][j]->nrows=0; smat[i][j]->ncols=0;
			//cout << smat[i][j]->nrows << "\t" << smat[i][j]->ncols<< endl;
			//smat[i][j]->regular2Sparse(realmat, lenx, leny);
			//testSMAT.regular2Sparse(realmat, lenx, leny);
			//cout << "Here " << endl;
			smat[j][i] = smat[i][j]->transpose();
			//cout << "Here " << endl;
			delete hmmProfPair;
			//cout << "Here " << endl;
			

		}
	}
	if(debug>1) {
	cout << "Time after the initial generation of consistency matrices: " << endl;
	times(&tmsend); timeDiff();
	times(&tmsstart);
	}
	
	//relaxConsistMatrix();
	//relaxConsistMatrix();
	for(i=0;i<relax_number;i++) {
	    //cout << " Before relax" << endl;
	    relaxConsistMatrix(dist_matrix, max_dist_cutoff);
	    //cout << " After relax" << endl;
	}
	/*
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
		    cout << "####### " << i << "   " << j << endl;
		    cout << preAligned[i]->aln->aname[0] << " " << preAligned[j]->aln->aname[0] << endl;
		    smat[i][j]->printSparseMatrix(1);
		}
	}
	*/
}

// unfinished
template <typename TNODE>
void btree<TNODE>::profileConsistency_profilehmm(hmm_parameters *params, char *ss_dir_name, int use_ss, float ss_weight, double **dist_matrix, double max_dist_cutoff) {


	int i,j=0,k,l,m,n;
	float **realmat, **realmat1;
	int lenx, leny, lenx1, leny1;
	int nonZeroCounts = 0;
	
	smat = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);
	smat1 = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);
        char tmpstr[200];

	// pairwise profile alignments; generate probability matrices (sparse)
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
			//cout << "i: " << i << "\t" << "j: " << "\t" << j << endl;
			nonZeroCounts = 0;
			//hmm_multim *hmmProfPair = new hmm_multim(preAligned[i]->aln, preAligned[j]->aln);
			//hmmProfPair->set_parameters(params);
			//hmm_profpair1 *hmmProfPair = new hmm_profpair1(preAligned[i]->aln, preAligned[j]->aln);
		        //hmmProfPair->viterbi();
			//cout << "Pairwise viterbi alignment: " << i << "\t" << j << endl;
		        //subalign *newalign = hmmProfPair->productViterbiAlign();
		        //newalign->printali(60);
			//delete newalign;
		        //cout << "------------" << endl;

			profilehmm *hmmProfPair = new profilehmm(preAligned[i]->similarSet);
			hmmProfPair->set_align(preAligned[j]->similarSet);
			//hmmProfPair->set_parameters(params, ss_dir_name, use_ss);
			hmmProfPair->set_parameters(tmpstr, ss_dir_name, use_ss);
			//hmmProfPair->x->get_score_bg(use_ss);
			hmmProfPair->ss_weight = ss_weight;
			hmmProfPair->get_scores(use_ss);

			hmmProfPair->log_convert();
			hmmProfPair->forward();
			hmmProfPair->backward();
			hmmProfPair->get_match_prob();
			realmat = hmmProfPair->probMat;
			lenx = hmmProfPair->lenx; leny = hmmProfPair->leny;
			if(reverse_align_order) {
				profilehmm *hmmProfPair1 = new profilehmm(preAligned[j]->aln);
				hmmProfPair1->set_align(preAligned[i]->similarSet);
				//hmmProfPair1->set_parameters(params, ss_dir_name, use_ss);
				hmmProfPair1->set_parameters(tmpstr, ss_dir_name, use_ss);
				//hmmProfPair1->x->get_score_bg(use_ss);
				hmmProfPair1->ss_weight = ss_weight;
				hmmProfPair1->get_scores(use_ss);
	
				hmmProfPair1->log_convert();
				hmmProfPair1->forward();
				hmmProfPair1->backward();
				hmmProfPair1->get_match_prob();
				realmat1 = hmmProfPair1->probMat;
				lenx1 = hmmProfPair1->lenx; leny1 = hmmProfPair1->leny;

				//cout << lenx1 << "  " << leny1 << endl;
				//cout << lenx << "  " << leny << endl;
			        for(k=1;k<=lenx1;k++) {
       	        		        for(l=1;l<=leny1;l++) {
					    realmat[l][k] = (realmat1[k][l]+realmat[l][k])/2;
					}
                		}
				delete hmmProfPair1;
        		}
				
		        for(k=1;k<=lenx;k++) {
               		        for(l=1;l<=leny;l++) {
                       		    if(realmat[k][l]<minProb) realmat[k][l] = 0;
				    else nonZeroCounts++;
                       		    //cout << realmat[k][l] << " ";
                		} //cout << endl;
        		}
			if(debug>1) cout << "None-zero counts: " << i << "\t" << j << "\t" << nonZeroCounts << endl;
			//cout << "Here " << endl;
			smat[i][j] = new sparseMatrix(realmat, lenx, leny);
			//cout << "Here " << endl;
			if(Debug>1) smat[i][j]->printCrs();
			//smat[i][j] = new sparseMatrix;
			//smat[i][j]->nrows=0; smat[i][j]->ncols=0;
			//cout << smat[i][j]->nrows << "\t" << smat[i][j]->ncols<< endl;
			//smat[i][j]->regular2Sparse(realmat, lenx, leny);
			//testSMAT.regular2Sparse(realmat, lenx, leny);
			//cout << "Here " << endl;
			smat[j][i] = smat[i][j]->transpose();
			//cout << "Here " << endl;
			delete hmmProfPair;
			

		}
	}
	if(debug>1) {
	cout << "Time after the initial generation of consistency matrices: " << endl;
	times(&tmsend); timeDiff();
	times(&tmsstart);
	}
	
	//relaxConsistMatrix();
	//relaxConsistMatrix();
	for(i=0;i<relax_number;i++) {
	    //cout << " Before relax" << endl;
	    relaxConsistMatrix(dist_matrix, max_dist_cutoff);
	    //cout << " After relax" << endl;
	}
	/*
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
		    cout << "####### " << i << "   " << j << endl;
		    cout << preAligned[i]->aln->aname[0] << " " << preAligned[j]->aln->aname[0] << endl;
		    smat[i][j]->printSparseMatrix(1);
		}
	}
	*/
}


template <typename TNODE>
void btree<TNODE>::profileConsistency_multim2(hmm_parameters *params1, hmm_parameters *params2, sequences *tmpSeq2, double id_cutoff) {


	int i,j=0,k,l,m,n;
	float **realmat, **realmat1;
	int lenx, leny, lenx1, leny1;
	int nonZeroCounts = 0;
	double dist;
	
	smat = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);
	smat1 = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);

	// pairwise profile alignments; generate probability matrices (sparse)
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
			//cout << "i: " << i << "\t" << "j: " << "\t" << j << endl;
			nonZeroCounts = 0;
			dist = tmpSeq2->distMat[tmpSeq2->name2index.find(preAligned[i]->aln->aname[0])->second][tmpSeq2->name2index.find(preAligned[j]->aln->aname[0])->second];
			//cout << preAligned[i]->aln->aname[0] << " " << tmpSeq2->name2index.find(preAligned[i]->aln->aname[0])->second << " " << preAligned[j]->aln->aname[0] << " " << tmpSeq2->name2index.find(preAligned[j]->aln->aname[0])->second << " " << dist << endl;
			hmm_multim *hmmProfPair = new hmm_multim(preAligned[i]->aln, preAligned[j]->aln);
			if(dist>1-id_cutoff) hmmProfPair->set_parameters(params1);
			else hmmProfPair->set_parameters(params2);

			hmmProfPair->forward();
			hmmProfPair->backward();
			realmat = hmmProfPair->probMat;
			lenx = hmmProfPair->lenx; leny = hmmProfPair->leny;
			/*
			if(reverse_align_order) {
				hmm_multim *hmmProfPair1 = new hmm_multim(preAligned[j]->aln, preAligned[i]->aln);
				hmmProfPair1->set_parameters(params);
				hmmProfPair1->forward();
				hmmProfPair1->backward();
				realmat1 = hmmProfPair1->probMat;
				lenx1 = hmmProfPair1->lenx; leny1 = hmmProfPair1->leny;
				//cout << lenx1 << "  " << leny1 << endl;
				//cout << lenx << "  " << leny << endl;
			        for(k=1;k<=lenx1;k++) {
       	        		        for(l=1;l<=leny1;l++) {
					    realmat[l][k] = (realmat1[k][l]+realmat[l][k])/2;
					}
                		}
				delete hmmProfPair1;
        		}
			*/	
		        for(k=1;k<=lenx;k++) {
               		        for(l=1;l<=leny;l++) {
                       		    if(realmat[k][l]<minProb) realmat[k][l] = 0;
				    else nonZeroCounts++;
                       		    //cout << realmat[k][l] << " ";
                		} //cout << endl;
        		}
			if(debug>1) cout << "None-zero counts: " << i << "\t" << j << "\t" << nonZeroCounts << endl;
			//cout << "Here " << endl;
			smat[i][j] = new sparseMatrix(realmat, lenx, leny);
			if(Debug>1) smat[i][j]->printCrs();
			//smat[i][j] = new sparseMatrix;
			//smat[i][j]->nrows=0; smat[i][j]->ncols=0;
			//cout << smat[i][j]->nrows << "\t" << smat[i][j]->ncols<< endl;
			//smat[i][j]->regular2Sparse(realmat, lenx, leny);
			//testSMAT.regular2Sparse(realmat, lenx, leny);
			//cout << "Here " << endl;
			smat[j][i] = smat[i][j]->transpose();
			//cout << "Here " << endl;
			delete hmmProfPair;
			

		}
	}
	if(debug>1) {
	cout << "Time after the initial generation of consistency matrices: " << endl;
	times(&tmsend); timeDiff();
	times(&tmsstart);
	}
	
	//relaxConsistMatrix();
	//relaxConsistMatrix();
	for(i=0;i<relax_number;i++) {
	    relaxConsistMatrix();
	}
	/*
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
		    cout << "####### " << i << "   " << j << endl;
		    cout << preAligned[i]->aln->aname[0] << " " << preAligned[j]->aln->aname[0] << endl;
		    smat[i][j]->printSparseMatrix(1);
		}
	}
	*/
}


template <typename TNODE>
void btree<TNODE>::profileConsistency_local() {

	int i,j=0,k,l,m,n;
	float **realmat;
	int lenx, leny;
	int nonZeroCounts = 0;
	
	smat = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);
	smat1 = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);

	// pairwise profile alignments; generate probability matrices (sparse)
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
			nonZeroCounts = 0;
			hmm_local *hmmProfPair = new hmm_local(preAligned[i]->aln, preAligned[j]->aln);
			hmmProfPair->forward();
			hmmProfPair->backward();
			realmat = hmmProfPair->probMat;
			lenx = hmmProfPair->lenx; leny = hmmProfPair->leny;
		        for(k=1;k<=lenx;k++) {
               		        for(l=1;l<=leny;l++) {
                       		    if(realmat[k][l]<minProb) realmat[k][l] = 0;
				    else nonZeroCounts++;
                		}
        		}
			if(debug>1) cout << "None-zero counts: " << i << "\t" << j << "\t" << nonZeroCounts << endl;
			smat[i][j] = new sparseMatrix(realmat, lenx, leny);
			if(Debug>1) smat[i][j]->printCrs();
			smat[j][i] = smat[i][j]->transpose();
			delete hmmProfPair;
		}
	}
	cout << "Time after the initial generation of consistency matrices: " << endl;
	times(&tmsend); timeDiff();
	times(&tmsstart);
	
	//relaxConsistMatrix();
	//relaxConsistMatrix();
	for(i=0;i<relax_number;i++) {
	    relaxConsistMatrix();
	}
}

template <typename TNODE>
void btree<TNODE>::profileConsistency_glocal(float wg) {

	int i,j=0,k,l,m,n;
	float **realmat;
	int lenx, leny;
	
	smat = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);
	smat1 = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);

	cout << "weightG: " << wg << endl;

	// pairwise profile alignments; generate probability matrices (sparse)
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
			hmm_local *hmmProfPair = new hmm_local(preAligned[i]->aln, preAligned[j]->aln);
			hmmProfPair->forward();
			hmmProfPair->backward();
			
			hmm_profpair1 *hmmProfPair_g = new hmm_profpair1(preAligned[i]->aln, preAligned[j]->aln);
			hmmProfPair_g->forward();
			hmmProfPair_g->backward();
			realmat = hmmProfPair->probMat;
			lenx = hmmProfPair->lenx; leny = hmmProfPair->leny;
			for(k=1;k<=lenx;k++) {
				for(l=1;l<=leny;l++) {
				   realmat[k][l] = wg * hmmProfPair_g->probMat[k][l] + (1-wg) * realmat[k][l];
				}
			}
		        for(k=1;k<=lenx;k++) {
               		        for(l=1;l<=leny;l++) {
                       		    if(realmat[k][l]<minProb) realmat[k][l] = 0;
                		}
        		}
			smat[i][j] = new sparseMatrix(realmat, lenx, leny);
			if(Debug>1) smat[i][j]->printCrs();
			smat[j][i] = smat[i][j]->transpose();
			delete hmmProfPair;
			delete hmmProfPair_g;
		}
	}
	if(debug>1) {
	cout << "Time after the initial generation of consistency matrices: " << endl;
	times(&tmsend); timeDiff();
	times(&tmsstart);
	}
	
	//relaxConsistMatrix();
	//relaxConsistMatrix();
	for(i=0;i<relax_number;i++) {
	    relaxConsistMatrix();
	}
}

// consistency derived from a set of alignments
template <typename TNODE> 
void btree<TNODE>::alignment_consistency(vector<sequences> seq_vector) {


	int i,j=0,k,l,m,n;
	float **realmat;
	int lenx, leny;
	
	smat = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);
	smat1 = gmatrix<sparseMatrix *>(preAligned.size()-1, preAligned.size()-1);

	for(i=0;i<seq_vector.size();i++) {
		seq_vector[i].get_map();
	}

	// pairwise profile alignments; generate probability matrices (sparse)
	int count1, count2;
	for(i=0;i<(int)preAligned.size();i++) {
		for(j=i+1;j<(int)preAligned.size();j++) {
			lenx = preAligned[i]->aln->alilen; leny = preAligned[j]->aln->alilen;
			realmat = gmatrix<float>(lenx, leny);
			for(k=0;k<=lenx;k++) {for(l=0;l<=leny;l++) realmat[k][l] = 0;}

			string name_string1(preAligned[i]->aln->aname[0]);
			string name_string2(preAligned[j]->aln->aname[0]);
			for(k=0;k<seq_vector.size();k++) {
				seq_vector[k].name2seq;
				string seq_string1 = seq_vector[k].name2seq.find(name_string1)->second;
				string seq_string2 = seq_vector[k].name2seq.find(name_string2)->second;

				count1 = count2 = 0;
				for(l=0;l<seq_string1.size();l++) {
					if((seq_string1[l]>='A')&&(seq_string1[l]<='Z')) {
						count1+=1;
						if((seq_string2[l]>='A')&&(seq_string2[l]<='Z')) {
							count2+=1;
							realmat[count1][count2] += 1;
						}
					}
					else {
						if((seq_string2[l]>='A')&&(seq_string2[l]<='Z')) {
							count2+=1;
						}
					}
				}
			}

			if(Debug>1) cout << i << "\t" << j << endl;
			for(l=1;l<=lenx;l++) {
				for(m=1;m<=leny;m++) {
					realmat[l][m] /= seq_vector.size();
					//cout << l << " " << m << " " << realmat[l][m] << endl;
					if(Debug>1) fprintf(stdout, "%3.2f ", realmat[l][m]);
				}
				if(Debug>1) fprintf(stdout, "\n");
			}
			
		        for(k=1;k<=lenx;k++) {
               		        for(l=1;l<=leny;l++) {
                       		    if(realmat[k][l]<minProb) realmat[k][l] = 0;
                		}
        		}
			smat[i][j] = new sparseMatrix(realmat, lenx, leny);
			if(Debug>1) smat[i][j]->printCrs();
			smat[j][i] = smat[i][j]->transpose();
			free_gmatrix<float>(realmat, lenx, leny);
		}
	}
	if(debug>1) {
	cout << "Time after the initial generation of consistency matrices: " << endl;
	}
	
	for(i=0;i<relax_number;i++) {
	    relaxConsistMatrix();
	}
}


template <typename TNODE> 
void btree<TNODE>::relaxConsistMatrix() {

	int i, j, k,l;
	float **tmpMat;

	if(debug>1) {
	cout << "Number of elements in the matrices before consistency: " << endl;
	for(i=0;i<preAligned.size();i++) { // Update the matrices
		for(j=0;j<preAligned.size();j++) {
			if(i==j) { cout << "0\t"; if(i==preAligned.size()-1) cout << endl; continue; }
			cout << smat[i][j]->nelements<<"\t";
		}
		cout << endl;
	}
	}
	// consistency measurement among the probability matrices
	for(i=0;i<preAligned.size();i++) {
		for(j=i+1;j<preAligned.size();j++) {
			//cout << "I: " << i << "\tJ: " << j <<endl;
			//smat1[i][j] = new sparseMatrix(&smat[i][j]);
			tmpMat = smat[i][j]->sparseCrs2Regular(); 
			// x-x-y and x-y-y
			for(k=1;k<=preAligned[i]->aln->alilen;k++) {
				for(l=1;l<=preAligned[j]->aln->alilen;l++) {
					tmpMat[k][l] *= 2;
				}
			}
			// relaxation 
			for(k=0;k<preAligned.size();k++) {
				if(k==i) continue; if(k==j) continue;
				relaxTwoSparse(smat[i][k], smat[k][j], tmpMat);
			}
			//cout << "HERE: " << endl;
			// normalize
		        for(k=1;k<=preAligned[i]->aln->alilen;k++) {
               		        for(l=1;l<=preAligned[j]->aln->alilen;l++) {
				    tmpMat[k][l] /= preAligned.size();
				    //if(relax_number>1) if(tmpMat[k][l]<minProb) tmpMat[k][l] = 0;
				    if(tmpMat[k][l]<minProb) tmpMat[k][l] = 0;
                		}
        		}
			//cout << "HERE: " << endl;
			smat1[i][j] = new sparseMatrix(tmpMat,preAligned[i]->aln->alilen, preAligned[j]->aln->alilen); 
			//cout << "HERE: " << endl;
			// normalize
			//smat1[i][j]->multiplyConstant(1.0/preAligned.size() );
			smat1[j][i] = smat1[i][j]->transpose();
			//smat1[i][j]->printCrs();
			free_gmatrix<float>(tmpMat, preAligned[i]->aln->alilen, preAligned[j]->aln->alilen);
			tmpMat = NULL;
		}
	}
	for(i=0;i<preAligned.size();i++) { // Update the matrices
		for(j=i+1;j<preAligned.size();j++) {
			//cout << "II: "<< i << "\tJJ: " << j << endl;
			smat[i][j]->clear();
			smat[i][j] = smat1[i][j];
			if(Debug>1) smat[i][j]->printCrs();
			smat[j][i]->clear();
			smat[j][i] = smat1[j][i];
		}
	}
	if(debug>1) {
	cout << "Number of elements in the matrices after consistency: " << endl;
	for(i=0;i<preAligned.size();i++) { // Update the matrices
		for(j=0;j<preAligned.size();j++) {
			if(i==j) { cout << "0\t"; if(i==preAligned.size()-1) cout << endl; continue; }
			cout << smat[i][j]->nelements<<"\t";
		}
		cout << endl;
	}
	cout << "Time after the relaxation: " << endl;
	times(&tmsend); timeDiff();
	}

} 

static int relaxcount = 0;
template <typename TNODE> 
void btree<TNODE>::relaxConsistMatrix(double **dist_matrix, double max_dist_cutoff, double min_cutoff) {

	int i, j, k,l;
	float **tmpMat;

        relaxcount += 1;
	int totalElements = 0;
	if(debug>1) {
	cout << "Number of elements in the matrices before consistency: " << endl;
	for(i=0;i<preAligned.size();i++) { // Update the matrices
		for(j=0;j<preAligned.size();j++) {
			if(i==j) { cout << "0\t"; if(i==preAligned.size()-1) cout << endl; continue; }
			cout << smat[i][j]->nelements<<"\t";
		}
		cout << endl;
	}
	}
	// consistency measurement among the probability matrices
	for(i=0;i<preAligned.size();i++) {
		for(j=i+1;j<preAligned.size();j++) {
			//cout << "I: " << i << "\tJ: " << j <<endl;
			//smat1[i][j] = new sparseMatrix(&smat[i][j]);
			tmpMat = smat[i][j]->sparseCrs2Regular(); 
                        float sum_of_elements = 0;
			// x-x-y and x-y-y
			for(k=1;k<=preAligned[i]->aln->alilen;k++) {
				for(l=1;l<=preAligned[j]->aln->alilen;l++) {
                                        sum_of_elements += tmpMat[k][l];
					tmpMat[k][l] *= 2;
				}
			}
			// relaxation 
			int number_for_relax=2;
			for(k=0;k<preAligned.size();k++) {
				if(k==i) continue; if(k==j) continue;
				if( (dist_matrix[i][k]==0)&&(dist_matrix[j][k]==0) ) continue;
				if(dist_matrix[i][k]>max_dist_cutoff) continue;
				if(dist_matrix[j][k]>max_dist_cutoff) continue;
				relaxTwoSparse(smat[i][k], smat[k][j], tmpMat);
				//relaxTwoSparse(*smat[i][k], *smat[k][j], tmpMat);
				number_for_relax+=1;
			}
			if(Debug>1) cout << i << " " << j << " Number for relax: " << number_for_relax << endl;
			//cout << "HERE: " << endl;
			// normalize
                        float sum_of_elements_after=0;
		        for(k=1;k<=preAligned[i]->aln->alilen;k++) {
               		        for(l=1;l<=preAligned[j]->aln->alilen;l++) {
				    //tmpMat[k][l] /= preAligned.size();
				    tmpMat[k][l] /= number_for_relax;
				    //if(relax_number>1) if(tmpMat[k][l]<min_cutoff) tmpMat[k][l] = 0;
                                    //tmpMat[k][l] *= 1.5;
				    if(tmpMat[k][l]<min_cutoff) tmpMat[k][l] = 0;
                                    sum_of_elements_after += tmpMat[k][l];
                		}
        		}
                        cout << "sum_of_elements: " << sum_of_elements << " " << sum_of_elements_after << " " << relaxcount << endl;
			//cout << "HERE: " << endl;
			smat1[i][j] = new sparseMatrix(tmpMat,preAligned[i]->aln->alilen, preAligned[j]->aln->alilen); 
			//cout << "HERE: " << endl;
			// normalize
			//smat1[i][j]->multiplyConstant(1.0/preAligned.size() );
			smat1[j][i] = smat1[i][j]->transpose();
			//smat1[i][j]->printCrs();
			free_gmatrix<float>(tmpMat, preAligned[i]->aln->alilen, preAligned[j]->aln->alilen);
			tmpMat = NULL;
			// NEW debug
			if(debug>1) cout << "Elements after relax: " << smat1[i][j]->nelements << endl;
			totalElements += smat1[i][j]->nelements;
		}
	}
	// NEW debug
	if(debug>1) cout << "totalElements: " << totalElements << endl;
	for(i=0;i<preAligned.size();i++) { // Update the matrices
		for(j=i+1;j<preAligned.size();j++) {
			//cout << "II: "<< i << "\tJJ: " << j << endl;
			//smat[i][j]->clear();
			delete smat[i][j];
			smat[i][j] = smat1[i][j];
			if(Debug>1) smat[i][j]->printCrs();
			delete smat[j][i];
			//smat[j][i]->clear();
			smat[j][i] = smat1[j][i];
		}
	}
	if(debug>1) {
	cout << "Number of elements in the matrices after consistency: " << endl;
	for(i=0;i<preAligned.size();i++) { // Update the matrices
                if(i>10) continue;
		for(j=0;j<preAligned.size();j++) {
			if(i==j) { cout << "0\t"; if(i==preAligned.size()-1) cout << endl; continue; }
                        cout << "element " << i << " and element " << j << endl;
			cout << smat[i][j]->nelements<<"\n";
                        smat[i][j]->printCrs();
		}
		cout << endl;
	}
	cout << "Time after the relaxation: " << endl;
	times(&tmsend); timeDiff();
	}

} 

template <typename TNODE> 
void btree<TNODE>::relaxConsistMatrix(double **dist_matrix, double max_dist_cutoff) {

	int i, j, k,l;
	float **tmpMat;

	int totalElements = 0;
	if(debug>1) {
	cout << "Number of elements in the matrices before consistency: " << endl;
	for(i=0;i<preAligned.size();i++) { // Update the matrices
		for(j=0;j<preAligned.size();j++) {
			if(i==j) { cout << "0\t"; if(i==preAligned.size()-1) cout << endl; continue; }
			cout << smat[i][j]->nelements<<"\t";
		}
		cout << endl;
	}
	}
	// consistency measurement among the probability matrices
	for(i=0;i<preAligned.size();i++) {
		for(j=i+1;j<preAligned.size();j++) {
			//cout << "I: " << i << "\tJ: " << j <<endl;
			//smat1[i][j] = new sparseMatrix(&smat[i][j]);
			tmpMat = smat[i][j]->sparseCrs2Regular(); 
			// x-x-y and x-y-y
			for(k=1;k<=preAligned[i]->aln->alilen;k++) {
				for(l=1;l<=preAligned[j]->aln->alilen;l++) {
					tmpMat[k][l] *= 2;
				}
			}
			// relaxation 
			int number_for_relax=2;
			for(k=0;k<preAligned.size();k++) {
				if(k==i) continue; if(k==j) continue;
				if( (dist_matrix[i][k]==0)&&(dist_matrix[j][k]==0) ) continue;
				if(dist_matrix[i][k]>max_dist_cutoff) continue;
				if(dist_matrix[j][k]>max_dist_cutoff) continue;
				relaxTwoSparse(smat[i][k], smat[k][j], tmpMat);
				number_for_relax+=1;
			}
			if(Debug>1) cout << i << " " << j << " Number for relax: " << number_for_relax << endl;
			//cout << "HERE: " << endl;
			// normalize
		        for(k=1;k<=preAligned[i]->aln->alilen;k++) {
               		        for(l=1;l<=preAligned[j]->aln->alilen;l++) {
				    //tmpMat[k][l] /= preAligned.size();
				    tmpMat[k][l] /= number_for_relax;
				    //if(relax_number>1) if(tmpMat[k][l]<minProb) tmpMat[k][l] = 0;
				    if(tmpMat[k][l]<minProb) tmpMat[k][l] = 0;
                		}
        		}
			//cout << "HERE: " << endl;
			smat1[i][j] = new sparseMatrix(tmpMat,preAligned[i]->aln->alilen, preAligned[j]->aln->alilen); 
			//cout << "HERE: " << endl;
			// normalize
			//smat1[i][j]->multiplyConstant(1.0/preAligned.size() );
			smat1[j][i] = smat1[i][j]->transpose();
			//smat1[i][j]->printCrs();
			free_gmatrix<float>(tmpMat, preAligned[i]->aln->alilen, preAligned[j]->aln->alilen);
			tmpMat = NULL;
			// NEW debug
			if(debug>1) cout << "Elements after relax: " << smat1[i][j]->nelements << endl;
			totalElements += smat1[i][j]->nelements;
		}
	}
	// NEW debug
	if(debug>1) cout << "totalElements: " << totalElements << endl;
	for(i=0;i<preAligned.size();i++) { // Update the matrices
		for(j=i+1;j<preAligned.size();j++) {
			//cout << "II: "<< i << "\tJJ: " << j << endl;
			//smat[i][j]->clear();
			delete smat[i][j];
			smat[i][j] = smat1[i][j];
			if(Debug>1) smat[i][j]->printCrs();
			delete smat[j][i];
			//smat[j][i]->clear();
			smat[j][i] = smat1[j][i];
		}
	}
	if(debug>-1) {
	cout << "Number of elements in the matrices after consistency: " << endl;
	for(i=0;i<preAligned.size();i++) { // Update the matrices
		for(j=0;j<preAligned.size();j++) {
			if(i==j) { cout << "0\t"; if(i==preAligned.size()-1) cout << endl; continue; }
			cout << smat[i][j]->nelements<<"\n";
                        smat[i][j]->printCrs();
		}
		cout << endl;
	}
	cout << "Time after the relaxation: " << endl;
	times(&tmsend); timeDiff();
	}

}

template <typename TNODE> 
void btree<TNODE>::computeConsistencyAlignment(TNODE *a) {

	int i,j,k;
	int **consistScoringMatrix;
	int *path;
	int len;
	int alignmentLength=0;
	int oi, ai;

	if(a->aligned) return;

	if(!a->childL->aligned) computeConsistencyAlignment(a->childL);
	if(!a->childR->aligned) computeConsistencyAlignment(a->childR);

	consistScoringMatrix = computeConsistMatrix(a->childL, a->childR);
	path = computePairwiseAlignment(consistScoringMatrix, a->childL->absSeqlength, a->childR->absSeqlength, len);
	//free_imatrix(consistScoringMatrix, a->childL->absSeqlength, a->childR->absSeqlength);

	//get the abstractSeq for the node a
	if(debug>1) cout << "len: " << len << endl;
	for(i=1;i<=len;i++) {
		if(path[i]==0) alignmentLength++; 
		else {
			alignmentLength += abs(path[i]);
		}
	}
	if(debug>1) cout << "alignmentLength: " << alignmentLength << endl;
	a->abstractSeq.push_back(0);
	if(debug>1) cout << "Number of absSeq: " << a->childL->absSeqnum+a->childR->absSeqnum<<endl;
	for(i=1;i<=a->childL->absSeqnum;i++) {
		int *newSeq = ivector(alignmentLength);
		oi=0, ai=0;
		newSeq[0] = a->childL->abstractSeq[i][0];
		for(j=1;j<=len;j++) {
			if(path[j]==0) {
				oi++; ai++;
				newSeq[ai] = a->childL->abstractSeq[i][oi];
			}
			else if(path[j]<0) {
				for(k=1;k<=abs(path[j]);k++) {
					oi++; ai++;
					newSeq[ai] = a->childL->abstractSeq[i][oi];
				}
			}
			else {
				for(k=1;k<=abs(path[j]);k++) {
					ai++;
					newSeq[ai] = 0;
				}
			}
		}
		if(debug>1) cout << "sequence i: " << i << endl;
		//for(j=1;j<=alignmentLength;j++) { cout << "j: " << j << "\t" << newSeq[j]; << endl; }
		if(debug>1) {for(j=1;j<=alignmentLength;j++) {cout << newSeq[j] << " ";} cout << endl; }
		a->abstractSeq.push_back(newSeq);
	}
	for(i=a->childL->absSeqnum+1;i<=a->childL->absSeqnum+a->childR->absSeqnum;i++) {
		int *newSeq = ivector(alignmentLength);
		oi=0, ai=0;
		if(debug>1) cout << "===========" <<endl;
		newSeq[0] = a->childR->abstractSeq[i-a->childL->absSeqnum][0];
		if(debug>1) cout << "===========" <<endl;
		for(j=1;j<=len;j++) {
			//cout << "j: " << j << "\t" << path[j] << endl;
			if(path[j]==0) {	
				oi++; ai++;
				newSeq[ai] = a->childR->abstractSeq[i-a->childL->absSeqnum][oi];
			}
			else if(path[j]>0) {
				for(k=1;k<=abs(path[j]);k++) {
					oi++; ai++;
					newSeq[ai] = a->childR->abstractSeq[i-a->childL->absSeqnum][oi];
				}
			}
			else {
				for(k=1;k<=abs(path[j]);k++) {
					ai++;
					newSeq[ai] = 0;
				}
			}
		}
		//cout << endl;
		if(debug>1) cout << "sequence i: " << i << endl;
		//for(j=1;j<=alignmentLength;j++) { cout << "j: " << j << "\t" << newSeq[j] << endl; }
		if(debug>1) for(j=1;j<=alignmentLength;j++) {cout << newSeq[j] << " ";} if(debug>1) cout << endl;
		a->abstractSeq.push_back(newSeq);
	}
	a->absSeqnum = a->childL->absSeqnum + a->childR->absSeqnum;
	a->absSeqlength = alignmentLength;
	//cout << "Number of letters for node a: " << a->absSeqnum << " * " << a->absSeqlength << " + " << a->absSeqlength*a->absSeqnum<<endl;
	if(debug>1) cout << "Obtain the abstract alignment" << endl;
	//for(i=1;i<=a->absSeqnum;i++) { for(j=0;j<=a->absSeqlength;j++) { cout << "IJ: " << i << "\t" << j << "\t" << a->abstractSeq[i][j]<< endl; } }
        if(debug>1) {
                cout << "printAlignmentFromAbs" << endl;
                printAlignmentFromAbs(a->childL);
                printAlignmentFromAbs(a->childR);
	        printAlignmentFromAbs(a);
        }
	if(debug>1) cout << "----------" <<endl;

	free_imatrix(consistScoringMatrix, a->childL->absSeqlength, a->childR->absSeqlength); 
	delete [] path;

	// print the intermediate results
	//cout << "Here: " << endl;
	//printAlignmentFromAbs(a);
	return;

}
template <typename TNODE> 
void btree<TNODE>::computeConsistencyAlignment(TNODE *a, float divergent_cutoff) {

	int i,j;
	int **consistScoringMatrix;
	int *path;
	int len;
	int alignmentLength=0;
	int oi, ai;

	if(a->aligned) return;

	if(!a->childL->aligned) computeConsistencyAlignment(a->childL);
	if(!a->childR->aligned) computeConsistencyAlignment(a->childR);

	// calculate the distance from the leaf to the TNODE a
	float dist=0;
	TNODE *k = a;
	while(k->childL!=0) {
		dist += k->childL->branchlen;
		k = k->childL;
	}
	//cout << "dist: " << dist << endl;
	if(dist > divergent_cutoff) return;

	

	consistScoringMatrix = computeConsistMatrix(a->childL, a->childR);
	path = computePairwiseAlignment(consistScoringMatrix, a->childL->absSeqlength, a->childR->absSeqlength, len);

	//get the abstractSeq for the node a
	if(debug>1) cout << "len: " << len << endl;
	for(i=1;i<=len;i++) {
		if(path[i]==0) alignmentLength++; 
		else {
			alignmentLength += abs(path[i]);
		}
	}
	if(debug>1) cout << "alignmentLength: " << alignmentLength << endl;
	a->abstractSeq.push_back(0);
	if(debug>1) cout << "Number of absSeq: " << a->childL->absSeqnum+a->childR->absSeqnum<<endl;
	for(i=1;i<=a->childL->absSeqnum;i++) {
		int *newSeq = ivector(alignmentLength);
		oi=0, ai=0;
		newSeq[0] = a->childL->abstractSeq[i][0];
		for(j=1;j<=len;j++) {
			if(path[j]==0) {
				oi++; ai++;
				newSeq[ai] = a->childL->abstractSeq[i][oi];
			}
			else if(path[j]<0) {
				for(k=1;k<=abs(path[j]);k++) {
					oi++; ai++;
					newSeq[ai] = a->childL->abstractSeq[i][oi];
				}
			}
			else {
				for(k=1;k<=abs(path[j]);k++) {
					ai++;
					newSeq[ai] = 0;
				}
			}
		}
		if(debug>1) cout << "sequence i: " << i << endl;
		//for(j=1;j<=alignmentLength;j++) { cout << "j: " << j << "\t" << newSeq[j]; << endl; }
		if(debug>1) {for(j=1;j<=alignmentLength;j++) {cout << newSeq[j] << " ";} cout << endl; }
		a->abstractSeq.push_back(newSeq);
	}
	for(i=a->childL->absSeqnum+1;i<=a->childL->absSeqnum+a->childR->absSeqnum;i++) {
		int *newSeq = ivector(alignmentLength);
		oi=0, ai=0;
		if(debug>1) cout << "===========" <<endl;
		newSeq[0] = a->childR->abstractSeq[i-a->childL->absSeqnum][0];
		if(debug>1) cout << "===========" <<endl;
		for(j=1;j<=len;j++) {
			//cout << "j: " << j << "\t" << path[j] << endl;
			if(path[j]==0) {	
				oi++; ai++;
				newSeq[ai] = a->childR->abstractSeq[i-a->childL->absSeqnum][oi];
			}
			else if(path[j]>0) {
				for(k=1;k<=abs(path[j]);k++) {
					oi++; ai++;
					newSeq[ai] = a->childR->abstractSeq[i-a->childL->absSeqnum][oi];
				}
			}
			else {
				for(k=1;k<=abs(path[j]);k++) {
					ai++;
					newSeq[ai] = 0;
				}
			}
		}
		//cout << endl;
		if(debug>1) cout << "sequence i: " << i << endl;
		//for(j=1;j<=alignmentLength;j++) { cout << "j: " << j << "\t" << newSeq[j] << endl; }
		if(debug>1) for(j=1;j<=alignmentLength;j++) {cout << newSeq[j] << " ";} if(debug>1) cout << endl;
		a->abstractSeq.push_back(newSeq);
	}
	a->absSeqnum = a->childL->absSeqnum + a->childR->absSeqnum;
	a->absSeqlength = alignmentLength;
	if(debug>1) cout << "Obtain the abstract alignment" << endl;
	//for(i=1;i<=a->absSeqnum;i++) { for(j=0;j<=a->absSeqlength;j++) { cout << "IJ: " << i << "\t" << j << "\t" << a->abstractSeq[i][j]<< endl; } }
	//printAlignmentFromAbs(a);
	if(debug>1) cout << "----------" <<endl;

	free_imatrix(consistScoringMatrix, a->childL->absSeqlength, a->childR->absSeqlength); 
	delete [] path;
	

	return;

}

template <typename TNODE> 
void btree<TNODE>::iterativeRefinement(int maxround) {

	int i,j;
	int len;
	int alignmentLength=0;
	int oi, ai;

        if(maxround==0) return;

        // 
        if(debug>-1) cout << "before refinement: " << endl;
        if(debug>-1) printAlignmentFromAbs(root);

        // 1. sort the tnode array according to distance to the leaf
        // if the tree is made by UPGMA, then there is actually no need to sort, but anyway sort it
        int sorted_index[2*preAligned.size()];
        float dist_to_leaf[2*preAligned.size()];
        int count = 1;
        TNODE *tmptnode;
        //cout << "size: " <<  size << endl;
        for(i=1;i<=size;i++) {
                if(v[i]->abstractSeq.size()>0) {
                        if(v[i]->isRoot()) continue;
                        sorted_index[count] = i;
                        dist_to_leaf[count] = 0;
                        tmptnode = v[i];
                        while(tmptnode->childL!=NULL) {
                                dist_to_leaf[count] += tmptnode->childL->branchlen;
                                tmptnode = tmptnode->childL;
                        }
                        count++;
                }
        }
        count--;
        assert(count==2*preAligned.size()-2);
     if(debug>1){
        cout << "count: " << count << endl;
        cout << "2*preAligned.size()-1: " << 2*preAligned.size()-1 << endl;
        for(i=1;i<=count;i++) {
                cout << i << " " << sorted_index[i] << " " << dist_to_leaf[i] << endl;
        }
        sort2<float, int>(count, dist_to_leaf, sorted_index);
        cout << "after sorted" << endl;
        for(i=1;i<=count;i++) {
                cout << i << " " << sorted_index[i] << " " << dist_to_leaf[i] << endl;
        }
     }

        // 
        for(int round=1;round<=maxround;round++) {
                for(i=1;i<=count;i++) {
                        recomputeAlignment(v[sorted_index[i]]);
                }
        }

        // reorder the sequences in
        vector<int *> newabstractSeq;
        newabstractSeq.push_back(NULL);
        for(i=0;i<preAligned.size();i++) {
                for(j=1;j<=root->absSeqnum;j++) {
                }

        }
        return;

}

template <typename TNODE> 
void btree<TNODE>::recomputeAlignment(TNODE *a) {
        
        int i, j, k;

        // 1. get the partition of root abstractSeq into two parts
        // store partion in a array of 0's and 1's
        int *partition = ivector(root->absSeqnum);
        int nabsn1 = 0;
        for(i=1;i<=root->absSeqnum;i++) {
                partition[i] = 0;
                for(j=1;j<=a->absSeqnum;j++) {
                        if(root->abstractSeq[i][0] == a->abstractSeq[j][0]) {
                                partition[i] = 1;
                                nabsn1++;
                                break;
                        }
                }
        }

        //cout << "partition: " << endl;
        //for(i=1;i<=root->absSeqnum;i++) { cout << "i: " << i << " " << partition[i] << endl; }

        // 2. create two tnodes corresponding to the partition
        //    delete all gapped positions from the abstractSeq
        //    then get the scoring matrix
        TNODE *t0 = new TNODE();
        TNODE *t1 = new TNODE();
        t0->abstractSeq.push_back(NULL);
        t1->abstractSeq.push_back(NULL);
        t1->absSeqnum = nabsn1;
        t0->absSeqnum = root->absSeqnum-nabsn1;
        // 2.1 find all-gap positions, mark them in partition 1 and partition 0
        int allgap_mark1[root->absSeqlength+1];
        int allgap_mark0[root->absSeqlength+1];
        allgap_mark1[0] = 0; // 0 means is not all gapped region
        allgap_mark0[0] = 0;
        t1->absSeqlength = 0;
        for(i=1;i<=root->absSeqlength;i++) {
                allgap_mark1[i] = 1;
                for(j=1;j<=root->absSeqnum;j++) {
                        if(partition[j]==1) {
                                //cout << "i: " << i << " j: " << root->abstractSeq[j][i] <<" "<<t1->absSeqlength<< endl;
                                if(root->abstractSeq[j][i]) {allgap_mark1[i]=0; t1->absSeqlength++;break;}
                        }
                }
        }
        t0->absSeqlength = 0;
        for(i=1;i<=root->absSeqlength;i++) {
                allgap_mark0[i] = 1;
                for(j=1;j<=root->absSeqnum;j++) {
                        if(partition[j]==0) {if(root->abstractSeq[j][i]) {allgap_mark0[i]=0; t0->absSeqlength++;break;}}
                }
        }
        //cout <<"absSeqlength for root: " << root->absSeqlength << endl;
        //cout <<"absSeqlength for t1: " << t1->absSeqlength << endl;
        //cout <<"absSeqlength for t0: " << t0->absSeqlength << endl;
        // 2.2 assign abstractSeq for t1 and t0
        int count;
        for(i=1;i<=root->absSeqnum;i++) {
                count = 0;
                int *tmpabsSeq = ivector(root->absSeqlength);
                if(partition[i]==1) {
                        t1->abstractSeq.push_back(tmpabsSeq);
                        for(j=0;j<=root->absSeqlength;j++) {
                                if(!allgap_mark1[j]) {
                                        tmpabsSeq[count] = root->abstractSeq[i][j];
                                        count++;
                                }
                        }
                }
                if(partition[i]==0) {
                        t0->abstractSeq.push_back(tmpabsSeq);
                        for(j=0;j<=root->absSeqlength;j++) {
                                if(!allgap_mark0[j]) {
                                        tmpabsSeq[count] = root->abstractSeq[i][j];
                                        count++;
                                }
                        }
                }
        }
        if(debug>1) for(i=1;i<=t1->absSeqnum;i++) { for(j=0;j<=t1->absSeqlength;j++) { cout << "t1: " << j << " " << t1->abstractSeq[i][j] << endl; } }
        // 2.3 create score matrix
        int **consistScoringMatrix = computeConsistMatrix(t1, t0);
	int *path;
        int len;
	int alignmentLength=0;
	path = computePairwiseAlignment(consistScoringMatrix, t1->absSeqlength, t0->absSeqlength, len);

        // 3. make alignment
	//get the abstractSeq for the node a
	if(debug>1) {
                cout << "len: " << len << endl;
                for(i=1;i<=len;i++) {
                        cout << path[i] << endl;
                }
        }
	for(i=1;i<=len;i++) {
		if(path[i]==0) alignmentLength++; 
		else {
			alignmentLength += abs(path[i]);
		}
	}
	if(debug>1) cout << "alignmentLength: " << alignmentLength << endl;
        // now update the root abstractSeq
        // first delete the old one
        for(i=1;i<=root->absSeqnum;i++) {
                delete [] root->abstractSeq[i];
        }
        root->abstractSeq.clear();
	root->abstractSeq.push_back(0);
        int oi, ai;
	for(i=1;i<=t1->absSeqnum;i++) {
		int *newSeq = ivector(alignmentLength);
		oi=0, ai=0;
		newSeq[0] = t1->abstractSeq[i][0];
		for(j=1;j<=len;j++) {
			if(path[j]==0) {
				oi++; ai++;
				newSeq[ai] = t1->abstractSeq[i][oi];
			}
			else if(path[j]<0) {
				for(k=1;k<=abs(path[j]);k++) {
					oi++; ai++;
					newSeq[ai] = t1->abstractSeq[i][oi];
				}
			}
			else {
				for(k=1;k<=abs(path[j]);k++) {
					ai++;
					newSeq[ai] = 0;
				}
			}
		}
		if(debug>1) cout << "sequence i: " << i << endl;
		//for(j=1;j<=alignmentLength;j++) { cout << "j: " << j << "\t" << newSeq[j]; << endl; }
		if(debug>1) {for(j=1;j<=alignmentLength;j++) {cout << newSeq[j] << " ";} cout << endl; }
		root->abstractSeq.push_back(newSeq);
	}
	for(i=t1->absSeqnum+1;i<=t1->absSeqnum+t0->absSeqnum;i++) {
		int *newSeq = ivector(alignmentLength);
		oi=0, ai=0;
		if(debug>1) cout << "===========" <<endl;
		newSeq[0] = t0->abstractSeq[i-t1->absSeqnum][0];
		if(debug>1) cout << "===========" <<endl;
		for(j=1;j<=len;j++) {
			//cout << "j: " << j << "\t" << path[j] << endl;
			if(path[j]==0) {	
				oi++; ai++;
				newSeq[ai] = t0->abstractSeq[i-t1->absSeqnum][oi];
			}
			else if(path[j]>0) {
				for(k=1;k<=abs(path[j]);k++) {
					oi++; ai++;
					newSeq[ai] = t0->abstractSeq[i-t1->absSeqnum][oi];
				}
			}
			else {
				for(k=1;k<=abs(path[j]);k++) {
					ai++;
					newSeq[ai] = 0;
				}
			}
		}
		//cout << endl;
		if(debug>1) cout << "sequence i: " << i << endl;
		//for(j=1;j<=alignmentLength;j++) { cout << "j: " << j << "\t" << newSeq[j] << endl; }
		if(debug>1) for(j=1;j<=alignmentLength;j++) {cout << newSeq[j] << " ";} if(debug>1) cout << endl;
		root->abstractSeq.push_back(newSeq);
	}
	root->absSeqnum = t1->absSeqnum + t0->absSeqnum;
	root->absSeqlength = alignmentLength;
	if(debug>1)cout << "Number of letters for node root: " << root->absSeqnum << " * " << root->absSeqlength << " + " << root->absSeqlength*root->absSeqnum<<endl;
	if(debug>1) {cout << "Obtain the abstract alignment" << endl;
	//for(i=1;i<=a->absSeqnum;i++) { for(j=0;j<=a->absSeqlength;j++) { cout << "IJ: " << i << "\t" << j << "\t" << a->abstractSeq[i][j]<< endl; } }

	printAlignmentFromAbs(t1);
	printAlignmentFromAbs(t0);
	printAlignmentFromAbs(root);
	cout << "----------" <<endl;
        }

	free_imatrix(consistScoringMatrix, t1->absSeqlength, t0->absSeqlength); 
	delete [] path;

        delete t0;
        delete t1;

	// print the intermediate results
	//cout << "Here: " << endl;
	//printAlignmentFromAbs(a);
	return;

              
}
        
template <typename TNODE> 
int **btree<TNODE>::computeConsistMatrix(TNODE *a, TNODE *b) {

	int i,j,k,m;
	int lenx, leny;
	float **tmpMat;
	int **scoreMat;
	int index1, index2;
        float scaling_factor = 1000.;

        scaling_factor = scaling_factor/a->absSeqnum/b->absSeqnum;
	
	lenx = a->absSeqlength;
	leny = b->absSeqlength;

	if(debug>1) {
	cout << "Computing consistency matrix" <<endl;
	cout << "lenx: " << lenx << "\tleny: " << leny << endl;
	cout << "Seqnum: " << a->absSeqnum << "\tseqnum: " << b->absSeqnum<< endl;
	}
	
	tmpMat = gmatrix<float>(lenx, leny);
	for(i=1;i<=lenx;i++)for(j=1;j<=leny;j++)tmpMat[i][j]=0;
	scoreMat = imatrix(lenx, leny);

	for(i=1;i<=lenx;i++) {
	    if(Debug>1) cout << i<< ": ";
	    for(j=1;j<=leny;j++) {
		for(k=1;k<=a->absSeqnum;k++) {
			if(a->abstractSeq[k][i]==0) {
				//cout << "a->abstractSeq[k][i]==0" << endl;
				continue;
			}
			index1 = a->abstractSeq[k][0];
			for(m=1;m<=b->absSeqnum;m++) {
				if(b->abstractSeq[m][j]==0) continue;
				index2 = b->abstractSeq[m][0];
				//cout << "k: " << k << "\tm: " << m << "\t" << index1 << "\t" << index2 << "\t" << a->abstractSeq[k][i] << "\t" << b->abstractSeq[m][j] << endl;
				tmpMat[i][j] += smat[index1][index2]->getElement(a->abstractSeq[k][i], b->abstractSeq[m][j]);
				//cout << tmpMat[i][j] << endl;
			}
		}
		scoreMat[i][j] = (int) (scaling_factor * tmpMat[i][j]);
                if(scoreMat[i][j]<0) {cout << "less than zero: " << scoreMat[i][j] << endl; }
		if(Debug>1) cout <<  scoreMat[i][j] << " "; // << endl;
	    }
	    if(Debug>1) cout << endl;
	}
	if(debug>1) cout << "Matrix computation ends here" << endl;

	free_gmatrix<float>(tmpMat, lenx, leny);

	return scoreMat;

}

template <typename TNODE> 
int * btree<TNODE>::computePairwiseAlignment(int **scoreMat, int m, int n, int &len) {

	int i,j,k;
	int *path;

	MM galign;
	galign.setM(m);
	galign.setN(n);
	galign.set_g(1); // gap open penalty
	galign.set_h(1); // gap extension penalty



	if(debug>1) cout << "Computer pairwise consistency alignment" << endl;
	galign.dp(scoreMat);
	if(debug>1)cout << "Now here" << endl;
	path = new int [galign.print_ptr];
	for(i=1;i<galign.print_ptr;i++) {
		path[i] = galign.displ[i];
		if(debug>1)cout << "i: " << i << " " << path[i] << endl;
	}
	len = galign.print_ptr - 1;
	//cout << "len: " << len << endl;
	if(debug>1) cout << "Pairwise consistency alignment ends here" << endl;	
        int a1=0, a2=0;
        int maxprobs = 0;
        for(i=1;i<galign.print_ptr;i++) {
                if(path[i]>0) {a2+=path[i]; continue;}
                if(path[i]<0) {a1-=path[i]; continue;}
                a1++; a2++;
                maxprobs += scoreMat[a1][a2];
        }

        cout << "maximum probablity: " << maxprobs << " "<< galign.maxscore << " " << m << " " << n << " " << MIN(m, n) << endl;

	return path;
}

template <typename TNODE> 
void btree<TNODE>::printAlignmentFromAbs(TNODE *n, char *outfilename, vector<subalign *>similaraln, vector<char *>repnames) {

	int i, j,k,l,m;
	int seqNUM=0;
	subalign *tmpSeq = new subalign;
	char *seqStr;
	tmpSeq->nal = 0;
	tmpSeq->alilen = n->absSeqlength;
	tmpSeq->mnamelen = 0;
	subalign *tmpSeq1;
	//cout << "printSeq here" << endl;

	if(n->abstractSeq.empty() ) {
		cout << "No abstract sequences, print the subalign" << endl;
		//refinegap(n->aln, 0.5, 1, 1, 1);
		//delete_complete_gap_positions(n->aln);
		
		n->aln->printali(blocksize);
		n->aln->printali(outfilename, blocksize);
		cout << "------------------" << endl;
	 	return;
	}

	assert(n->absSeqnum >0);

	for(i=1;i<=n->absSeqnum;i++) { seqNUM += preAligned[n->abstractSeq[i][0]]->aln->nal; }
	//for(i=1;i<=n->absSeqnum;i++) { preAligned[n->abstractSeq[i][0]]->aln->printali(60); }
	//cout << "seqNUM: " << seqNUM << endl;
	tmpSeq->aseq = new char *[seqNUM+1];//cmatrix(0, seqNUM-1, 0, n->absSeqlength);
	for(i=0;i<=seqNUM;i++) tmpSeq->aseq[i] = new char [n->absSeqlength+1];
	tmpSeq->aname = new char *[seqNUM+1];
	for(i=0;i<=seqNUM;i++) tmpSeq->aname[i] = new char [100];
	//cout << "finished assign tmpSeq" << endl;

	if(debug>1) 
	for(i=1;i<=n->absSeqnum;i++) {
	    cout << "i: " << i << endl;
	    for(j=1;j<=n->absSeqlength;j++) {
	   	cout << n->abstractSeq[i][j] << " ";
	    }
	    cout << endl;
	}
	int S=0; for(i=1;i<=n->absSeqnum;i++) { S+=preAligned[n->abstractSeq[i][0]]->aln->nal; } tmpSeq->nal = S;
	if(debug>1) cout << "Nal: " << tmpSeq->nal << endl;

	//cout << n->absSeqnum << endl;

	S=0;
        /*
	for(i=1;i<=n->absSeqnum;i++) {
	   if(tmpSeq->mnamelen<preAligned[n->abstractSeq[i][0]]->aln->mnamelen) tmpSeq->mnamelen = preAligned[n->abstractSeq[i][0]]->aln->mnamelen;
	   for(j=1;j<=preAligned[n->abstractSeq[i][0]]->aln->nal;j++) {
		seqStr = tmpSeq->aseq[S]; //new char [n->absSeqlength+1];
		//seqStr = new char [n->absSeqlength+1];
		l = 0;
		for(k=1;k<=n->absSeqlength;k++) {
			if(n->abstractSeq[i][k]==0) {
				seqStr[k-1] = '-';
			}
			else {
				seqStr[k-1] = preAligned[n->abstractSeq[i][0]]->aln->aseq[j-1][n->abstractSeq[i][k]-1];
				l++;
			}
		}
		seqStr[k-1]='\0';
		assert(k-1==n->absSeqlength);
		strcpy(tmpSeq->aname[S], preAligned[n->abstractSeq[i][0]]->aln->aname[j-1]);
		//cout << left << setw(30) << tmpSeq->aname[S] << tmpSeq->aseq[S] << endl;
		
		//cout << setw(20) << seqStr << endl;
		//(tmpSeq->nal)++;
		S++;
	   }
	}
        */
        for(j=0;j<preAligned.size();j++) {
                for(i=1;i<=n->absSeqnum;i++) {
                        if(n->abstractSeq[i][0] == j) break;
                }
                if(tmpSeq->mnamelen<preAligned[j]->aln->mnamelen) tmpSeq->mnamelen = preAligned[j]->aln->mnamelen;
                for(m=1;m<=preAligned[j]->aln->nal;m++) {
                        seqStr = tmpSeq->aseq[S];
                        for(k=1;k<=n->absSeqlength;k++) {
                                if(n->abstractSeq[i][k] == 0) {seqStr[k-1] = '-';}
                                else { seqStr[k-1] = preAligned[j]->aln->aseq[m-1][n->abstractSeq[i][k]-1]; }
                        }
                        seqStr[k-1]='\0';
                        assert(k-1==n->absSeqlength);
                        strcpy(tmpSeq->aname[S], preAligned[j]->aln->aname[m-1]);
                        S++;
                }
        }
                                        
        /*
	for(i=1;i<=n->absSeqnum;i++) {
	   if(tmpSeq->mnamelen<preAligned[n->abstractSeq[i][0]]->aln->mnamelen) tmpSeq->mnamelen = preAligned[n->abstractSeq[i][0]]->aln->mnamelen;
	   for(j=1;j<=preAligned[n->abstractSeq[i][0]]->aln->nal;j++) {
		seqStr = tmpSeq->aseq[S]; //new char [n->absSeqlength+1];
		//seqStr = new char [n->absSeqlength+1];
		l = 0;
		for(k=1;k<=n->absSeqlength;k++) {
			if(n->abstractSeq[i][k]==0) {
				seqStr[k-1] = '-';
			}
			else {
				seqStr[k-1] = preAligned[n->abstractSeq[i][0]]->aln->aseq[j-1][n->abstractSeq[i][k]-1];
				l++;
			}
		}
		seqStr[k-1]='\0';
		assert(k-1==n->absSeqlength);
		strcpy(tmpSeq->aname[S], preAligned[n->abstractSeq[i][0]]->aln->aname[j-1]);
		//cout << left << setw(30) << tmpSeq->aname[S] << tmpSeq->aseq[S] << endl;
		
		//cout << setw(20) << seqStr << endl;
		//(tmpSeq->nal)++;
		S++;
	   }
	}
        */
        cout << "alignment of representatives is below" << endl;
        tmpSeq->printali(90);

	// add similar sequences
	int found = 0;
	//cout << found << endl;
        vector<subalign *> similarset_aln;
        vector<char *> similarset_repnames;
	for(i=0;i<preAligned.size();i++) {
		if(preAligned[i]->similarSet) {
			//cout << "similarSet " << i << endl;
                        //cout << preAligned[i]->similarSet->nal << endl;
			found = 0;
                        similarset_aln.push_back(preAligned[i]->similarSet);
			for(j=0;j<preAligned[i]->similarSet->nal;j++) {
				for(k=0;k<tmpSeq->nal;k++) {
					if(strcmp(preAligned[i]->similarSet->aname[j], tmpSeq->aname[k])==0) {
						found = 1;
						break;
					}
				}
				if(found==1) break;
			}
			//cout << tmpSeq->aname[k] << endl;
			if(found==0) {
				cout << "same name does not found" << endl;
				exit(0);
			}
                        similarset_repnames.push_back(tmpSeq->aname[k]);
			//tmpSeq1 = merge_align_by_one_sequence_insert(tmpSeq, preAligned[i]->similarSet, tmpSeq->aname[k]);	
			//delete tmpSeq;
			//***old *** for(k=0;k<=tmpSeq->nal;k++) {delete [] tmpSeq->aseq[k]; delete [] tmpSeq->aname[k];}
			//***old *** delete [] tmpSeq->aseq; delete [] tmpSeq->aname;
			//delete tmpSeq;
			//tmpSeq = tmpSeq1;
		}
	}
        if(similarset_aln.size()) {
                tmpSeq1 = merge_master_and_slaves(tmpSeq, similarset_aln, similarset_repnames);
                delete tmpSeq;
                tmpSeq = tmpSeq1;
                similarset_aln.clear();
                similarset_repnames.clear();
                print_time_diff("add_similarset");
        }

        // add similaraln
        /*
        if(!exclude_similar) {
                for(i=0;i<similaraln.size();i++) {
                        tmpSeq1 = merge_align_by_one_sequence_insert(tmpSeq, similaraln[i], repnames[i]);
                        delete tmpSeq;
                        tmpSeq = tmpSeq1;
                }
        }
        */

	tmpSeq->alignment = imatrix(tmpSeq->nal, tmpSeq->alilen);
	for(i=1;i<=tmpSeq->nal;i++) {
	    for(j=1;j<=tmpSeq->alilen;j++) {
		tmpSeq->alignment[i][j] = am2num(tmpSeq->aseq[i-1][j-1]);
	    }
	}

        if(Debug>1) {
                cout << "alignment after adding similarset" << endl;
	        tmpSeq->printali(80);
        }

	
	//delete_complete_gap_positions(tmpSeq);
	fprintf(logfp, "\nStart refining alignment ...\n");
	fprintf(logfp, "\t- step 1\n");
                
	fflush(logfp);
	if( (tmpSeq->nal <= 1000)&&(tmpSeq->alilen<=4000) )  {refine_align_new(tmpSeq);cout<<"refinenew"<<endl;}
	fprintf(logfp, "\t- step 2\n");
	fflush(logfp);
        if(Debug>1) { cout << "alignment after refine_align_new, step 1" <<endl;  tmpSeq->printali(70);}
        print_time_diff("refine_align_new, step 1");

	refinegap(tmpSeq, 0.8, 1, 1, 1);
	fprintf(logfp, "\t- step 3\n");
	fflush(logfp);
	//cout << "after refine1: " << endl;
	delete_complete_gap_positions(tmpSeq);
        if(Debug>1) { cout << "alignment after refine_align_new, step 2" <<endl;  tmpSeq->printali(70);}
        print_time_diff("refinegap, step 2");

	if( (tmpSeq->nal <= 1000)&&(tmpSeq->alilen<=4000) ) {
		alignrefine *y = new alignrefine(tmpSeq, 0);
		y->treat_single_and_doublet();
                cout <<"refinenew2" << endl;
	}
        if(Debug>1) { cout << "alignment after alignrefine, step 3" <<endl;  tmpSeq->printali(70);}
        print_time_diff("alignrefine, step 3");
	//treat_single_residues(tmpSeq, 0.5);

	//gap_refine *gr = new gap_refine();
	//tmpSeq->aseq = gr->batch_refine(tmpSeq, 2);
	//tmpSeq->alilen = strlen(tmpSeq->aseq[0]);
	//delete gr;

	// now re-order sequences and add predicted secondary structure information
	subalign *tmpSeq_ss = new subalign;
	tmpSeq_ss->alilen = tmpSeq->alilen;
	tmpSeq_ss->nal = tmpSeq->nal + preAligned.size() + 1;
	tmpSeq_ss->aseq = cmatrix(tmpSeq_ss->nal, tmpSeq_ss->alilen);
        cout << "tmpSeq->mnamelen: " << tmpSeq->mnamelen+2 << endl;
	tmpSeq_ss->mnamelen = tmpSeq->mnamelen;
        if(tmpSeq_ss->mnamelen < strlen("Consensus_ss: ") ) tmpSeq_ss->mnamelen = strlen("Consensus_ss: ");
	tmpSeq_ss->aname = cmatrix(tmpSeq_ss->nal, tmpSeq_ss->mnamelen+2);
	int tmpSeqnum = 0;
	int tmppos = 0;
	char *consensus_ss = new char [tmpSeq_ss->alilen+1];
	int **consensus_counts = imatrix(tmpSeq_ss->alilen+1, 3);
        // a re-implementation, finding the indexes of the representatives and use them as marks of boundaries
        int *repindex = ivector(preAligned.size());
        for(j=0;j<tmpSeq->nal;j++) {
                for(i=0;i<preAligned.size();i++) {
                        if(strcmp(preAligned[i]->aln->aname[0], tmpSeq->aname[j])==0) {
                                repindex[i] = j;
                        }
                }
        }
        repindex[preAligned.size()] = tmpSeq->nal;
        int myindex = 0;
	for(i=0;i<preAligned.size();i++) {
		// find the name
		for(j=0;j<tmpSeq->nal;j++) {
			if(strcmp(preAligned[i]->aln->aname[0], tmpSeq->aname[j])==0) {break;}
		}

		// get secondary structure
		tmppos = 0;
		for(k=0;k<tmpSeq->alilen;k++) {
			if(tmpSeq->aseq[j][k]=='-') {
				tmpSeq_ss->aseq[tmpSeqnum][k] = '-';
			}
			else {
				tmpSeq_ss->aseq[tmpSeqnum][k] = ssint2ss(preAligned[i]->aln->ss->sstype[tmppos+1]);
				consensus_counts[k][preAligned[i]->aln->ss->sstype[tmppos+1]]++;
				tmppos++;
			}
		}
		tmpSeq_ss->aseq[tmpSeqnum][tmpSeq->alilen] = '\0';
		strcpy(tmpSeq_ss->aname[tmpSeqnum], "ss:");
		tmpSeqnum++;
                /*
		// get the representative 
		strcpy(tmpSeq_ss->aname[tmpSeqnum], tmpSeq->aname[j]);
		strcpy(tmpSeq_ss->aseq[tmpSeqnum], tmpSeq->aseq[j]);
		tmpSeqnum++;

		// get sequences from similar set
		if(!preAligned[i]->similarSet) continue;
		for(k=0;k<preAligned[i]->similarSet->nal;k++) {
			if(strcmp(preAligned[i]->similarSet->aname[k], tmpSeq->aname[j])==0) continue;
			for(l=0;l<tmpSeq->nal;l++) {
				if(strcmp(tmpSeq->aname[l], preAligned[i]->similarSet->aname[k])==0) break;
			}
			strncpy(tmpSeq_ss->aseq[tmpSeqnum], tmpSeq->aseq[l], tmpSeq_ss->alilen);
			strcpy(tmpSeq_ss->aname[tmpSeqnum], tmpSeq->aname[l]);
			//cout << tmpSeq->aname[l] << endl;
			//cout << tmpSeq_ss->aname[tmpSeqnum] << endl;
			tmpSeqnum++;
		}
                */
                // get sequences based on repindex
                //cout << repindex[i] << endl; cout << repindex[i+1] << endl;
                for(k=repindex[i];k<repindex[i+1];k++) {
                        strncpy(tmpSeq_ss->aseq[tmpSeqnum], tmpSeq->aseq[k], tmpSeq_ss->alilen);
                        tmpSeq_ss->aseq[tmpSeqnum][tmpSeq_ss->alilen] = '\0';
                        strcpy(tmpSeq_ss->aname[tmpSeqnum], tmpSeq->aname[k]);
                        tmpSeqnum++;
                }
	}
	// get the consensus secondary structure sequence
	int totalcounts = 0;
	for(i=0;i<tmpSeq->alilen;i++) {
		totalcounts = consensus_counts[i][1]+consensus_counts[i][2]+consensus_counts[i][3];
		// modified; totalcounts equals prealigned size
		totalcounts = preAligned.size();
		if(totalcounts==0) {consensus_ss[i] = '.'; continue;}
		if(consensus_counts[i][1]*1.0/totalcounts>=0.5) consensus_ss[i] = 'h';
		else if(consensus_counts[i][2]*1.0/totalcounts>=0.5) consensus_ss[i] = 'e';
		else consensus_ss[i] = '.';
	}
	consensus_ss[i] = '\0';
        cout << "consensus_secondary structure:" << endl << consensus_ss << endl;
	strcpy(tmpSeq_ss->aname[tmpSeqnum], "Consensus_ss:");
	strcpy(tmpSeq_ss->aseq[tmpSeqnum], consensus_ss);

	if(tmpSeq_ss->mnamelen<12) tmpSeq_ss->mnamelen=12;
                
        // now re-modified the sequences to remove small letters (which are added to representatives)
        for(i=0;i<tmpSeq->nal;i++) {
                for(j=0;j<tmpSeq->alilen;j++) {
                        if( (tmpSeq->aseq[i][j]>='a') && (tmpSeq->aseq[i][j]<='z') ) {
                                tmpSeq->aseq[i][j] = '-';
                        }
                }
        }
        for(i=0;i<tmpSeq_ss->nal;i++) {
                if(strcmp("ss:", tmpSeq_ss->aname[i])==0) continue;
                if(strcmp("Consensus_ss:", tmpSeq_ss->aname[i])==0) continue;
                for(j=0;j<tmpSeq_ss->alilen;j++) {
                        if( (tmpSeq_ss->aseq[i][j]>='a') && (tmpSeq_ss->aseq[i][j]<='z') ) {
                                tmpSeq_ss->aseq[i][j] = '-';
                                if(strcmp(tmpSeq_ss->aname[i-1], "ss:")==0) { // also modify the secondary structure
                                        tmpSeq_ss->aseq[i-1][j] = '-';
                                }
                        }
                }
        }
        print_time_diff("adding_ss");

        // add similaraln
        if(!exclude_similar) {
                print_section_info("Below add similaraln");
                // old version is below, add one by one, very slow
                /*
                for(i=0;i<similaraln.size();i++) {
                        tmpSeq1 = merge_align_by_one_sequence_insert(tmpSeq, similaraln[i], repnames[i]);
                        delete tmpSeq;
                        tmpSeq = tmpSeq1;
                }
                for(i=0;i<similaraln.size();i++) {
                        tmpSeq1 = merge_align_by_one_sequence_insert(tmpSeq_ss, similaraln[i], repnames[i]);
                        delete tmpSeq_ss;
                        tmpSeq_ss = tmpSeq1;
                }
                */

                // new version, add them at one time
             if(similaraln.size()) {
                tmpSeq1 = merge_master_and_slaves(tmpSeq, similaraln, repnames);
                delete tmpSeq;
                tmpSeq = tmpSeq1;
	        //tmpSeq->printali(outfilename, blocksize);
                tmpSeq1 = merge_master_and_slaves(tmpSeq_ss, similaraln, repnames);
                delete tmpSeq_ss;
                tmpSeq_ss = tmpSeq1;
                print_time_diff("adding_similaraln");
             }
        }
        cout << "This is the place"<< endl;

        // when adding similaraln, some gap characters "-" could be introduced into the consensus_ss line, 
        // replace them with "."
        for(i=0;i<tmpSeq_ss->nal;i++) {
                if(strcmp(tmpSeq_ss->aname[i], "Consensus_ss:")==0) {
                        for(j=0;j<tmpSeq_ss->alilen;j++) {
                                if (tmpSeq_ss->aseq[i][j] == '-') tmpSeq_ss->aseq[i][j] = '.';
                        }
                }
        }
	tmpSeq->printali(outfilename, blocksize);
        print_section_info("Below write alignments in log file and alignment file");
        cout << "Final alignment with secondary structure info is below"<< endl;
        cout << endl << "  output file Name: " << outfilename << endl << endl;
	tmpSeq_ss->printali(blocksize);
        cout << "  program finished"  << endl << endl;
        print_time_diff("print_alignments");
        /*
        for(i=0;i<tmpSeq_ss->nal;i++) {
                cout << "i: " << i << " " << tmpSeq_ss->aseq[i] << endl;
                cout << "strlen: " << strlen(tmpSeq_ss->aseq[i]) << endl;
                if(tmpSeq_ss->aseq[i][tmpSeq_ss->alilen]=='\0') {
                        cout << "good ending" << endl;
                }
                else {
                        cout << "bad ending" << endl;
                }
        }
        */
	//delete tmpSeq;
}

template <typename TNODE> 
void btree<TNODE>::refine_align_new(subalign *x) {

	int i, j, k;

	alignrefine *y = new alignrefine(x, 100);

	double threshold_here = 0.2;
        int core_num = 3;

	y->set_gappy_threshold(threshold_here);
        y->set_CORE_POS_NUMBER(core_num);

	y->assign_gi();
	y->get_sw();
	y->calculate_weights();

        int debug_this=1;

	y->deal_with_gappy(1);
        if(debug_this>1){ cout << "after 1:" << endl; x->printali(80);}
	y->deal_with_gappy(2);
        if(debug_this>1){ cout << "after 2:" << endl; x->printali(80);}
	y->deal_with_gappy(3);
        if(debug_this>1){cout << "after 3:" << endl; x->printali(80);}
	y->deal_with_gappy(4);
        if(debug_this>1){cout << "after 4:" << endl; x->printali(80);}
	y->deal_with_gappy(1000);
        if(debug_this>1){cout << "after 1000:" << endl; x->printali(80);}
	y->deal_with_N();
        if(debug_this>1){cout << "after N" << endl; x->printali(80);}
	y->deal_with_C();
        if(debug_this>1){cout << "after C" << endl; x->printali(80);}
	y->delete_NC_terminal_gaps();
        if(debug_this>1){cout << "after NC" << endl; x->printali(80);}

	delete y;
	
	alignrefine *y1 = new alignrefine(x, 0);
	y1->treat_single_and_doublet();

	delete y1;

	/*alignrefine *y2 = new alignrefine(x, 0);
	y2->treat_single_and_doublet();
	*/

}

template <typename TNODE> 
void btree<TNODE>::printAlignmentFromAbs(TNODE *n) {

	int i, j,k,l;
	int seqNUM=0;
	subalign *tmpSeq = new subalign;
	char *seqStr;
	tmpSeq->nal = 0;
	tmpSeq->alilen = n->absSeqlength;
	tmpSeq->mnamelen = 0;
	subalign *tmpSeq1;
	cout << "printSeq here" << endl;

	if(n->abstractSeq.empty() ) {
		cout << "No abstract sequences, print the subalign" << endl;
		//refinegap(n->aln, 0.5, 1, 1, 1);
		//delete_complete_gap_positions(n->aln);
		
		n->aln->printali(blocksize);
		//n->aln->printali(blocksize);
		cout << "------------------" << endl;
	 	return;
	}

	assert(n->absSeqnum >0);

	for(i=1;i<=n->absSeqnum;i++) { seqNUM += preAligned[n->abstractSeq[i][0]]->aln->nal; }
	//cout << "seqNUM: " << seqNUM << endl;
	tmpSeq->aseq = new char *[seqNUM+1];//cmatrix(0, seqNUM-1, 0, n->absSeqlength);
	for(i=0;i<=seqNUM;i++) tmpSeq->aseq[i] = new char [n->absSeqlength+1];
	tmpSeq->aname = new char *[seqNUM+1];
	for(i=0;i<=seqNUM;i++) tmpSeq->aname[i] = new char [100];

	if(debug>1) 
	for(i=1;i<=n->absSeqnum;i++) {
	    cout << "i: " << i << endl;
	    for(j=1;j<=n->absSeqlength;j++) {
	   	cout << n->abstractSeq[i][j] << " ";
	    }
	    cout << endl;
	}
	int S=0; for(i=1;i<=n->absSeqnum;i++) { S+=preAligned[n->abstractSeq[i][0]]->aln->nal; } tmpSeq->nal = S;
	if(debug>1) cout << "Nal: " << tmpSeq->nal << endl;

	S=0;
	for(i=1;i<=n->absSeqnum;i++) {
	   if(tmpSeq->mnamelen<preAligned[n->abstractSeq[i][0]]->aln->mnamelen) tmpSeq->mnamelen = preAligned[n->abstractSeq[i][0]]->aln->mnamelen;
	   for(j=1;j<=preAligned[n->abstractSeq[i][0]]->aln->nal;j++) {
		seqStr = tmpSeq->aseq[S]; //new char [n->absSeqlength+1];
		//seqStr = new char [n->absSeqlength+1];
		l = 0;
		for(k=1;k<=n->absSeqlength;k++) {
			if(n->abstractSeq[i][k]==0) {
				seqStr[k-1] = '-';
			}
			else {
				seqStr[k-1] = preAligned[n->abstractSeq[i][0]]->aln->aseq[j-1][n->abstractSeq[i][k]-1];
				l++;
			}
		}
		seqStr[k-1]='\0';
		assert(k-1==n->absSeqlength);
		strcpy(tmpSeq->aname[S], preAligned[n->abstractSeq[i][0]]->aln->aname[j-1]);
		//cout << left << setw(20) << tmpSeq->aname[S] << tmpSeq->aseq[S] << endl;
		
		//cout << setw(20) << seqStr << endl;
		//(tmpSeq->nal)++;
		S++;
	   }
	}
	if(debug>1) cout << "before adding similar sequences:" << endl;
	if(debug>-1) tmpSeq->printali(blocksize+20);
        return;

	// add similar sequences
	int found = 0;
	//cout << found << endl;
	//for(i=0;i<preAligned.size();i++) {
	for(l=1;l<=n->absSeqnum;l++) {
		i = n->abstractSeq[l][0];
		//cout << "i: " << i << endl;
		if(preAligned[i]->similarSet) {
			//cout << "similarSet " << i << endl;
			found = 0;
			for(j=0;j<preAligned[i]->similarSet->nal;j++) {
				for(k=0;k<tmpSeq->nal;k++) {
					if(strcmp(preAligned[i]->similarSet->aname[j], tmpSeq->aname[k])==0) {
						found = 1;
						break;
					}
				}
				if(found==1) break;
			}
			//cout << tmpSeq->aname[k] << endl;
			if(found==0) {
				cout << "same name does not found" << endl;
				exit(0);
			}
			tmpSeq1 = merge_align_by_one_sequence(tmpSeq, preAligned[i]->similarSet, tmpSeq->aname[k]);	
			//delete tmpSeq;
			tmpSeq = tmpSeq1;
		}
	}

	tmpSeq->alignment = imatrix(tmpSeq->nal, tmpSeq->alilen);
	for(i=1;i<=tmpSeq->nal;i++) {
	    for(j=1;j<=tmpSeq->alilen;j++) {
		tmpSeq->alignment[i][j] = am2num(tmpSeq->aseq[i-1][j-1]);
	    }
	}

	
	//refinegap(tmpSeq, 0.5, 1, 1, 1);
	delete_complete_gap_positions(tmpSeq);
	//treat_single_residues(tmpSeq, 0.5);

	//gap_refine *gr = new gap_refine();
	//tmpSeq->aseq = gr->batch_refine(tmpSeq, 2);
	//tmpSeq->alilen = strlen(tmpSeq->aseq[0]);
	//delete gr;

	if(debug>1) cout << "after adding similar sequences:" << endl;
	if(debug>1) tmpSeq->printali(blocksize);

	//delete tmpSeq;
}

template <typename TNODE> 
void btree<TNODE>::get_descendants(TNODE *r) {

	int i, j, k;

	if(r->childL) {
		get_descendants(r->childL);
		get_descendants(r->childR);
		r->descendants = r->childL->descendants + r->childR->descendants;
	}
	else { 
		r->descendants = 1;
	}
}

template <typename TNODE> 
void btree<TNODE>::assign_weights(TNODE *r) {

	TNODE *tmp;
	if(r->childL) {
		r->seq_weight = 0;
		assign_weight(r->childL);
		assign_weight(r->childR);
	}

	else {
		r->seq_weight = 0;
		tmp = r;
		while(tmp->rootFlag) {
			r->seq_weight += tmp->branchlen()/tmp->descendants;
		}
	}
}

#endif
