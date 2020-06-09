#ifndef __param_
#define __param_


#include "all.h"

extern double ave_grp_thr;
extern double minProb;
extern string outFile;
extern int probconsBLOSUM;

extern void getParameter(int argc, char **argv, int prog);
extern void printParameters();
extern void printParameters1();
extern void printHelp(int prog);
extern int useLocal;
extern float weightG;
extern int solv;
extern int ss;
extern int unaligned;
extern char parameter_file[];
extern char parameter_file1[];
extern char parameter_file2[];
extern int relax_number;
extern int reverse_align_order;
extern double id_thr;

extern char program_dir[500];

// for psipred
extern int psipred_env_number;
extern char psipred_dir[500];
extern char psipred_parameter_file[200];
extern char runpsipred_command[500];
extern char runpsipred1_command[500];

// for database homolog
extern char blast_dir[500];
extern char blastpgp_cmd[500];
extern char uniref90_file[500];
extern char blastpgp_command[500];
extern char blosum62_file[500];

// for using ss-dependent amino acid frequencies
extern int use_ss_freq; // 0, 1, or 2

extern int use_single_sequence;
extern float ss_w;
extern float score_w;
extern float score_shift;
extern int adjust_weight;

// blastpgp options
extern int iteration_number;
extern double evalue;

// purge blast output option
extern int max_num_sequences;
extern double low_id_thr;

// output alignment
extern int blocksize;

//clean blastpgp and psipred intermediate results
extern int clean_blast_after;
extern int clean_blast_before;

extern int relax_count; // correponding to N_small in multiple.[ch]

// maximum number of pre-aligned groups
extern int max_group_number;

#endif
