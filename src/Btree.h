#ifndef _tree__
#define _tree__

#include "util.h"
#include "amino.h"

// This tree node type only allows three children
class tnode {

    public:
	char name[50];
	tnode *childL;
	tnode *childR;
	tnode *childT; // third child
	tnode *parent;    
	bool rootFlag;
	double branchlen;
	int n;  // position in the array of pointers to tnodes
	
	double dist;
	int seqIndex;
	int len;
	char *aseq;
	int *seq;

 	tnode();
	tnode(const tnode &a);
	~tnode();

	bool isInternal();
	bool isRoot();

	static int nodecount;
	
};

class btree {

      public:
	tnode *root;

	// constructors
	btree();
	btree(const btree &a);
	~btree();

	int size;
	bool unrooted;
	char *treename;

	// array of pointers to the tnodes 
	tnode **v;

	// read tree from a file
	void readTree(char *treefile);
	void writeTree(char *outfile);
	void writeTopology(char *outfile);

	// reroot the tree
	void reroot(tnode *);
	// root the tree between r and its parent
	void rootTree(tnode *, double ); 
	void unrootTree();
        void leastSquareRoot(tnode **rootNode, double *dist, double *var, int &Num);
	
	static int btreecount;

     private: 
	void readTreefile(ifstream &tf, tnode *rt);
	void writeTreefile(ofstream &otf, tnode *rt);
	void writeTopologyfile(ofstream &otf, tnode *rt);
	void map_v();
	void map_v_recurse(tnode *r, int &index);
	
	void checkRoot();
	
};

#endif

