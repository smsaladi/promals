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
	
};

#endif
