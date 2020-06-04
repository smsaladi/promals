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
	if(Debug>1) newalign->printali(70);
	if(Debug>1) cout << "------------" << endl;
	gap_refine *gr = new gap_refine();
	newalign->aseq = gr->batch_refine(newalign, 2);
	delete_complete_gap_positions(newalign);
	newalign->alilen = strlen(newalign->aseq[0]);
	newalign->convertAseq2Alignment();
	n->aln = newalign;
	n->aligned = 1;
	delete gr;
	if(Debug>1) n->aln->printali(70);
	if(Debug>1) cout << "------------" << endl;
}
			
template <typename TNODE>
void  btree<TNODE>::obtainPreAligned(TNODE *n) {

	int i;
	if(n==0) return;
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
	for(i=0;i<relax_number;i++) {
	    //cout << " Before relax" << endl;
	    relaxConsistMatrix(dist_matrix, max_dist_cutoff);
	    fprintf(logfp, "        finished relaxing consistency matrix - round %d of %d\n", i+1, relax_number);
	    fflush(logfp);
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
			hmmProfPair->set_parameters(params, ss_dir_name, use_ss);
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
				hmmProfPair1->set_parameters(params, ss_dir_name, use_ss);
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

	// print the intermediate results
	printAlignmentFromAbs(a);
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

	cout << "before refine: "<< endl;
	tmpSeq->printali(80);

	
	refine_align_new(tmpSeq);
	cout << "after refine1: " << endl;
	tmpSeq->printali(80);
	refinegap(tmpSeq, 0.8, 1, 1, 1);
	delete_complete_gap_positions(tmpSeq);
	cout << "after refine2: " << endl;
	tmpSeq->printali(80);
	alignrefine *y = new alignrefine(tmpSeq, 0);
	y->treat_single_and_doublet();
	cout << "after refine3: " << endl;
	tmpSeq->printali(80);
	//treat_single_residues(tmpSeq, 0.5);

	//gap_refine *gr = new gap_refine();
	//tmpSeq->aseq = gr->batch_refine(tmpSeq, 2);
	//tmpSeq->alilen = strlen(tmpSeq->aseq[0]);
	//delete gr;

	tmpSeq->printali(outfilename, blocksize);

	// now re-order sequences and add predicted secondary structure information
	subalign *tmpSeq_ss = new subalign;
	tmpSeq_ss->alilen = tmpSeq->alilen;
	tmpSeq_ss->nal = tmpSeq->nal + preAligned.size() + 1;
	tmpSeq_ss->aseq = cmatrix(tmpSeq_ss->nal, tmpSeq_ss->alilen);
	tmpSeq_ss->aname = cmatrix(tmpSeq_ss->nal, tmpSeq->mnamelen+2);
	int tmpSeqnum = 0;
	int tmppos = 0;
	tmpSeq_ss->mnamelen = tmpSeq->mnamelen;
	char *consensus_ss = new char [tmpSeq_ss->alilen+1];
	int **consensus_counts = imatrix(tmpSeq_ss->alilen+1, 3);
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
		// get the representative 
		strcpy(tmpSeq_ss->aname[tmpSeqnum], tmpSeq->aname[j]);
		strcpy(tmpSeq_ss->aseq[tmpSeqnum], tmpSeq->aseq[j]);
		tmpSeqnum++;

		if(!preAligned[i]->similarSet) continue;
		// get sequences from similar set
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
	strcpy(tmpSeq_ss->aname[tmpSeqnum], "Consensus_ss:");
	strcpy(tmpSeq_ss->aseq[tmpSeqnum], consensus_ss);

	if(tmpSeq_ss->mnamelen<12) tmpSeq_ss->mnamelen=12;

	tmpSeq_ss->printali(blocksize);
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

	y->deal_with_gappy(1);
	y->deal_with_gappy(2);
	y->deal_with_gappy(3);
	y->deal_with_gappy(4);
	y->deal_with_gappy(1000);
	y->deal_with_N();
	y->deal_with_C();
	y->delete_NC_terminal_gaps();
	
	alignrefine *y1 = new alignrefine(x, 0);
	y1->treat_single_and_doublet();
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
	//cout << "printSeq here" << endl;

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
	if(debug>1) tmpSeq->printali(blocksize);

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
