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

static int Debug = 1;

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
	if(Debug>1) n->aln->printali(60);
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

	if( (!n->childL->aligned) || (!n->childR->aligned) ) { return; }

	// determine the distance from the TNODE to the leafs
	float dist=0;
	TNODE *k = n;
	while(k->childL!=0) {
		dist += k->childL->branchlen;
		k = k->childL;
	}
	//cout << "dist: " << dist << endl;
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
	n->aln = newalign;
	n->aligned = 1;
	if(Debug>1) n->aln->printali(60);
	if(Debug>1) cout << "------------" << endl;
}
			
template <typename TNODE>
void  btree<TNODE>::obtainPreAligned(TNODE *n) {

	int i;
	if(n==0) return;
	// if root is already aligned
	if( (n->aligned) && (n->rootFlag) ) {return;}
	if(n->aligned && (!n->parent->aligned) ) {
		TNODE *a = n;
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
		        if(Debug>1) newalign->printali(60);
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

template <typename TNODE> 
void btree<TNODE>::relaxConsistMatrix(double **dist_matrix, double max_dist_cutoff) {

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
int **btree<TNODE>::computeConsistMatrix(TNODE *a, TNODE *b) {

	int i,j,k,m;
	int lenx, leny;
	float **tmpMat;
	int **scoreMat;
	int index1, index2;
	
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
		scoreMat[i][j] = (int) (1000 * tmpMat[i][j]);
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
	galign.set_g(0); // gap open penalty
	galign.set_h(0); // gap extension penalty


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

	return path;
}

template <typename TNODE> 
void btree<TNODE>::printAlignmentFromAbs(TNODE *n, char *outfilename) {

	int i, j,k,l;
	int seqNUM=0;
	subalign *tmpSeq = new subalign;
	char *seqStr;
	tmpSeq->nal = 0;
	tmpSeq->alilen = n->absSeqlength;
	tmpSeq->mnamelen = 0;
	subalign *tmpSeq1;
	//cout << "printSeq here" << endl;

	if(n->abstractSeq.empty() ) {
		//cout << "No abstract sequences, print the subalign" << endl;
		//refinegap(n->aln, 0.5, 1, 1, 1);
		//delete_complete_gap_positions(n->aln);
		
		//n->aln->printali(60);
		n->aln->printali(outfilename, 60);
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

	// add similar sequences
	int found = 0;
	//cout << found << endl;
	for(i=0;i<preAligned.size();i++) {
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

	
	refinegap(tmpSeq, 0.5, 1, 1, 1);
	delete_complete_gap_positions(tmpSeq);
	tmpSeq->printali(outfilename, 60);

	delete tmpSeq;
	
	//cout <<"========"<<endl;
}

#endif
			
