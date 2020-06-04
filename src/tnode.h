#ifndef _tnode__
#define _tnode__

#include "util.h"
#include "amino.h"
#include "subalign.h"
#include "sequences.h"

// This tree node type only allows three children
class tnode {

    public:
	//char name[50];
	string name;
	tnode *childL;
	tnode *childR;
	tnode *childT; // third child
	tnode *parent;    
	bool rootFlag;
	double branchlen;
	int n;  // position in the array of pointers to tnodes
	int p;  // position in the vector of pointers to the pre-aligned groups
	int p_seq; // position in the vector of allseqs
	
	double dist;
	int seqIndex;
	int len;
	char *aseq;
	int *seq;

	subalign *aln;
	int aligned; // to check if there is an alignment

 	tnode();
	tnode(const tnode &a);
	~tnode();

	bool isInternal();
	bool isRoot();

	static int nodecount;

	void getSubalign();

	// abstract sequences
	vector<int *> abstractSeq;
	int absSeqnum;
	int absSeqlength;

	subalign *similarSet;
	subalign *consiston;

	subalign *aux_align;

	// secondary structure profile
	ss_prof *ssp;
	//void get_ssp() { ssp = aln->get_ssp(); }

	// weight for the branch above it
	double branch_weight;
	// weight for the sequence in a leaf node
	double seq_weight;
	// number of descendants 
	int descendants;
	
};

#endif
