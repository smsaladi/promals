#ifndef __regularizer__
#define __regularizer__

#include "header_cpp.h"


extern void read_regularizer(char *filename);
extern void print_regularizer();

// transition probabilities from or to begin and end states
extern float rsn;
extern float rsb;
extern float rnb;
extern float rnn;
extern float rec;
extern float ret;
extern float rcc;
extern float rct;

// transition probabilities from or to match states
extern float rbm[12];
extern float rme[12];

// transition probabilities among match, insert, delete states
extern float **rtrans0;
extern float ***rtrans;
//extern float rtrans[10][4][4];

// background residue pair emission probabilites
extern float **rqmatrix0;
extern float ***rqmatrix;
//extern float rqmatrix[10][21][21];

// background residue emission probabilites
extern float *rbfreq0;
extern float **rbfreq;
extern float *log_rbfreq0, **log_rbfreq;
//extern float rbfreq[10][21];

// predicted secondary structure background probabilities
extern float rssfreq[4], log_rssfreq[4];

extern int done_read_regularizer;

#endif
