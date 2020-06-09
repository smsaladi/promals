#ifndef __blastpgp_
#define __blastpgp_

#include "header_cpp.h"
#include "subalign.h"

extern subalign *get_blastpgp_alignment(char *basename, char *query_name,
                                        char *query_seq);

extern char uniref90_file[500];

extern void set_uniref90(char *file1);

extern int get_blastpgp_cmd();

extern char blastpgp_options[2000];

extern void set_blastpgp_options(char *options);

extern subalign *run_blastpgp(char *query_name, char *query_seq);

extern subalign *read_blastpgp_result(const char *blastpgp_result,
                                      char *query_name, char *query_seq);

extern void purge_alignment_by_first_sequence(subalign *x);

extern subalign *get_blastpgp_result(const char *blastpgp_result,
                                     char *query_name, char *query_seq);

extern subalign *read_blastpgp_alignment(const char *blastpgp_alignment,
                                         char *query_name, char *query_seq);

extern void clean_blast_psipred(char *seqname);

#endif
