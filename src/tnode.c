#include "header_cpp.h"
#include "util.h"
#include "amino.h"
#include "tnode.h"


// functions of tnode

int tnode::nodecount = 0;

tnode::tnode() {
	//name[0]= '\0';
   	childL = 0;
	childR = 0;
	childT = 0;
	parent = 0;
	branchlen = 0;
	n = -1;
	p = -1;

	aseq = 0;
	seq = 0;
	len = 0;

	aln = 0;
	aligned = 0;

	rootFlag = 0;
	nodecount++;
}

tnode::tnode(const tnode &a) {
	int i;
	//strcpy(name, a.name);
	name = a.name;
	childL = a.childL;
	childR = a.childR;
	childT = a.childT;
	parent = a.parent;
	rootFlag = a.rootFlag;
	branchlen = a.branchlen;
	
	dist = a.dist;
	seqIndex = a.seqIndex;
	n = a.n;
	p = a.p;
     if(a.len>0) {
	len = a.len;
	aseq = cvector(len+1);
	seq = ivector(len+1);
	strcpy(aseq, a.aseq);
	for(i=1;i<=len;i++) seq[i] = a.seq[i];
     }
     else {
	len = 0;
	aseq = 0; 
	seq = 0;
     }
	// unfinished here

	aln = a.aln;
	aligned = a.aligned;

	nodecount++;
}

tnode::~tnode() {
	if(aseq) delete [] aseq;
	if(seq) delete [] seq;
        childL = 0;
        childR = 0;
        childT = 0;
        parent = 0;
        aseq = 0;
        seq = 0;

	if(aln!=0) { delete aln; aln = 0; }

	nodecount--;
}

bool tnode::isInternal() {
	if(childL || childR) { return true;}
	else return false;
}

bool tnode::isRoot() {
	if(rootFlag) return true;
	else return false;
}
	
// obtain the subalign from the one sequence and the one sequence name in tnode
void tnode::getSubalign() {

	char str[1000];
	
	if(!aseq) {
		cout << "sequence does not exist in getSubalign" << endl;
		exit(0);
	}
	if(name.empty() ) {
		cout << "sequence name does not exist in getSubalign" << endl;
		exit(0);
	}

	strcpy(str,name.c_str ()); 
	aln = oneSeq2subalign(aseq, str );

	aligned = 1;


}
