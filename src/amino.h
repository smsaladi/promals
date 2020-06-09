#ifndef __AMINO_
#define __AMINO_

#include <cctype>
#include <string>
#include "header_cpp.h"
int a3let2num(char *let);
int am2num_c(int c);
int am2num(int c);
int am2num_dayhoff(int c);
int am2numBZX(int c);
int am2nm(int c);
int a3let2num(char *let);
void read_aa_matrix(char *filename, int **arr);

extern double dayhoff_freq[];
extern double robinson_freq[];
extern float log_robinson_freq[];
extern void get_log_robinson_freq();
extern int dayhoff_mutab[];
extern double hydrophobicity[];
extern double amino_charge[];
extern const char *am;
extern const char *am_dayhoff;
extern const char *am3[];
extern double q_blosum62[21][21];
extern float log_q_blosum62[21][21];
void get_log_q_blosum62();
extern float log_q_blosum62_ratio[21][21];
void get_log_q_blosum62_ratio();

// probcons parameters
void useProbconsFrequencies();
extern float emitPairsDefault[20][20];
extern string alphabetDefault;
extern float emitSingleDefault[20];
	    
#endif
