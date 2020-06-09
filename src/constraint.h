#ifndef __constraint__
#define __constraint__

#include "sequences.h"
#include "tnode.h"

// class constraint is inherited from class sequences
class constraint : public sequences {
 public:
  constraint() { nseqs = 0; }
  ~constraint();
  vector<tnode *> *prealigned;  // pointer to the vector of pre-aligned groups
  sequences *oseq;              // pointer to original sequences
  int *rep_index;               // prealigned group number for each sequence
  int *inside_index;            // the index within the prealigned group
  int *startp;  // starting position of the sequence in the original sequence

  // void read_seqs(char *file_list_name);
  void assign_prealigned(vector<tnode *> *prealn);
  void assign_originalseq(sequences *os);
  void allocate_index();
  void checkNamesSequences();
  string removegap(string a1, int touppercase);
  void printIndexSeq();
  float **get_constraint_matrix(int ci, int cj);
};

extern int *get_correspond(char *tmpseq1, char *tmpseq2, int &len);
extern vector<constraint *> read_multiple_constraint(char *filename);

#endif
