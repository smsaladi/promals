#ifndef __param_
#define __param_

#include "header_cpp.h"

extern double ave_grp_thr;
extern double minProb;
extern string outFile;
extern int probconsBLOSUM;

extern void getParameter(int argc, char **argv);
extern void printParameters();
extern void printHelp();
extern int solv;
extern int ss;
extern int unaligned;
extern int relax_number;
extern int reverse_align_order;
extern double id_thr;

extern char program_dir[500];

// for psipred
extern int psipred_env_number;
extern char psipred_dir[500];
extern char psipred_parameter_file[200];
extern char runpsipred_command[500];

// for database homolog
extern char blast_dir[500];
extern char blastpgp_cmd[500];
extern char uniref90_file[500];
extern char blastpgp_command[500];

// for using ss-dependent amino acid frequencies
extern int use_ss_freq;  // 0, 1, or 2

extern float ss_w;
extern float score_w;
extern float score_shift;
extern int adjust_weight;

// blastpgp options
extern int iteration_number;
extern double evalue;
extern char filter_query[3];

// purge blast output option
extern int max_num_sequences;
extern double low_id_thr;

// output alignment
extern int blocksize;

// clean blastpgp and psipred intermediate results
extern int clean_blast_after;
extern int clean_blast_before;

extern int relax_count;  // correponding to N_small in multiple.[ch]

// maximum number of pre-aligned groups
extern int max_group_number;

//
extern int struct_id_cutoff;
extern int below_id_cutoff;

//
extern int before_relax_combine;

//
extern float struct_weight;
extern float sequence_weight;
extern float user_constraint_weight;
extern float pdb_weight;
extern double minDaliZ;

// use different structural alignments
extern int use_dali;
extern int use_fast;
extern int use_tmalign;

// realign psiblast alignment
extern int realign_psiblast;

// constraint file
extern char constraint_file[500];
extern char user_constraint[500];

// weight end penalty
extern double weight_end_penalty;

//
extern char mafft[500];

// filter similar sequences before alignment
extern int filter_similar;
extern char cdhit[500];
extern float cdhit_c_option;
extern int exclude_similar;

extern int use_updated_database;

extern int max_struct_number;
#endif
