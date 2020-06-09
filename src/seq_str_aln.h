#ifndef __seq_str_aln__
#define __seq_str_aln__

#include "header_cpp.h"
#include "subalign.h"

// sequence alignments between a query sequence and hits with known structures
// from SCOP40

class seq_str_aln {
 public:
  seq_str_aln(char *q);
  seq_str_aln();
  ~seq_str_aln();

  vector<int *> aln;             // from 1..nhits
  vector<char *> id;             // from 1..nhits
  vector<int> start;             // from 1..nhits
  vector<int> end;               // from 1..nhits
  vector<int> str_start;         // from 1..nhits
  vector<int> str_end;           // from 1..nhits
  vector<subalign *> prof;       // from 1..nhits
  vector<subalign *> oneseqaln;  // from 1..nhits
  vector<int> slen;              // from 1..nhits, subject length

  int nhits;
  int subject_N_C_extension;
  void set_subject_N_C_extension(int a) { subject_N_C_extension = a; }

  char *query;
  char query_name[50];
  int len;
  char blastoutput[200];
  int id_cutoff;
  int below_id_cutoff;

  void set_id_cutoff(int a);
  void set_below_id_cutoff(int a);
  void run_blast(char *blastexe, char *fasta_file, char *database, char *suffix,
                 char *chkfile, char *options);
  void read_blast_results();
  int checkboundary(int qstart, int qend);
  void get_prof();
  void get_prof_update();
  void print_result();
};

// void add_structural_info(btree *tree);

#endif
