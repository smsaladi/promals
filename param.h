#ifndef __param_
#define __param_


#include "all.h"

extern double ave_grp_thr;
extern double minProb;
extern string outFile;
extern int probconsBLOSUM;

extern void getParameter(int argc, char **argv, int prog);
extern void printParameters();
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

#endif
