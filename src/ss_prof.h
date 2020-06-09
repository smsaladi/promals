#ifndef __ss_profile__
#define __ss_profile__

#include "header_cpp.h"
#include "util.h"

// predicted secondary structure of PSIPRED
class ss_prof {

    public:
	int sslen;
	int *sstype;
	float **ssfreq;
	int *ssrel;
	int *alphabet1;
	char *seq;

	int valid_file;

	ss_prof(char *base_name);
	~ss_prof();

	void read_ss_files(char *base_name);

	void get_alphabet1();

	void print_ss_info();

	int done_prof;

	int check_aa_seq(char *target_seq);
};
	
#endif
