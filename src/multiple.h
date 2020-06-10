#ifndef __multiple__
#define __multiple__

#include "blastpgp.h"
#include "btree_template.h"
#include "constraint.h"
#include "header_cpp.h"
#include "seq_str_aln.h"
#include "sequences.h"
#include "subalign.h"

extern FILE *logfp;

class replacedletters {
 public:
  char *name;
  int *pos;
  short len;
};

class multiple {
 public:
  multiple(char *filename);

  char inputfileName[200];

  double distance_cutoff_similar;
  double distance_cutoff_div;

  sequences allseqs;

  btree<tnode> alltree;

  void build_tree();
  void build_cluster();
  void alignSimilar();


  void store_similar(tnode *r);
  void store_similar_henikoff(tnode *r);

  void map_allseqs_pos_to_tnode();

  int N_small;
  double **dist_matrix_preAligned;
  void get_distance_matrix_for_preAligned(int N_smallest);

  sequences aligned_seqs;

  void set_distance_cutoff_similar(double dist_cutoff);
  void set_distance_cutoff_similar(double dist_cutoff, int Ngroup);
  void distance2leaf(tnode *r, double *array);
  void output_alignment();

  void addSimilar();

  void alignDivergent_psipred(int x);

  vector<seq_str_aln *> ssaln;
  vector<seq_str_aln *> *get_seq_str_alns();
  vector<seq_str_aln *> *get_seq_str_alns1();
  void combine_structure_alignment();

  vector<subalign *> similaraln;
  vector<char *> repnames;

  void get_kcenter_clusters(int dim, int nc, int maxelem, int maxiter,
                            int maxinit);
  int max_cluster_elem;
  int max_iteration;
  int max_initiation;
  int *clusterindex;
  vector<subalign *> prealn;
  vector<char *> prealn_repnames;
  sequences repseqs;  // sequence object of representative sequences
  void cluster2tree_and_similarSet();

  vector<replacedletters> rl;
};

int **get_salnpos(char *name1, char *name2, int &numalignedp);
float **map_structure(seq_str_aln *a1, seq_str_aln *a2);
int **get_salnpos_fast(char *name1, char *name2, int &numalignedp);
float **map_structure_fast(seq_str_aln *a1, seq_str_aln *a2);
int **get_salnpos_tmalign(char *name1, char *name2, int &numalignedp);
float **map_structure_tmalign(seq_str_aln *a1, seq_str_aln *a2);
void combine_structmat(float **mat1, float **mat2, int d1, int d2,
                       float weight);
void combine_structure_alignments(btree<tnode> &alltree,
                                  vector<seq_str_aln *> &ssaln);
void combine_structure_alignments3(btree<tnode> &alltree,
                                   vector<seq_str_aln *> &ssaln);
void combine_constraint(btree<tnode> &alltree, constraint &cons,
                        float constraint_w);

float **map_structure_nopromals(seq_str_aln *a1, seq_str_aln *a2,
                                const char *prog_name);
float **map_structure_promals(subalign *auxa, subalign *auxb, seq_str_aln *a1,
                              seq_str_aln *a2, const char *prog_name);

subalign *merge_align_by_one_sequence(subalign *a, subalign *b, char *seq_name);
subalign *merge_align_by_one_sequence_insert(subalign *a, subalign *b,
                                             char *seq_name);
subalign *merge_master_and_slaves(subalign *a, vector<subalign *> &b,
                                  vector<char *> &seq_name);

extern char ssint2ss(int i);

// for file filter_similar
void do_filter_similar(sequences &allseqs, char *filename, float thr,
                       vector<subalign *> &similaraln);
subalign *get_similar_aln(vector<char *> seqnames, sequences &allseqs,
                          char *filefilter);
void get_rep_names(sequences allseqs, vector<subalign *> &similaraln,
                   vector<char *> &repnames);

int *kcenter(int **simmat, int dim, int nc, int maxelem, int maxiter,
             int maxinit);

void combine_structure_alignments3_updated(btree<tnode> &alltree,
                                           vector<seq_str_aln *> &ssaln);
int **get_salnpos_updated(seq_str_aln *name1, seq_str_aln *name2, int ind1,
                          int ind2, int &numalignedp, char *prog_name);
float **map_structure_nopromals_updated(seq_str_aln *a1, seq_str_aln *a2,
                                        const char *prog_name);
float **map_structure_promals_updated(subalign *auxa, subalign *auxb,
                                      seq_str_aln *a1, seq_str_aln *a2,
                                      char *prog_name);
float **map_structure_promals_updated(subalign *auxa, subalign *auxb,
                                      seq_str_aln *a1, seq_str_aln *a2,
                                      int dali, int fast, int tmalign);
#endif
