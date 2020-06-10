#ifndef __subalign_
#define __subalign_

#include <unistd.h>

#include "amino.h"
#include "regularizer.h"
#include "ss_prof.h"
#include "util.h"
#include "param.h"

class subalign {
 public:
  subalign();
  subalign(char *filename);
  subalign(string filename);
  // copy constructor
  subalign(const subalign &);
  subalign(char *filename, const char *format_name, int max_name_len);
  ~subalign();

  // assign
  const subalign &operator=(const subalign &);

  void readali(char *filename);
  void printali(int blocksize);
  void printali(int blocksize, int header);
  void printali(const char *filename, int blocksize);
  void printProfile();

  double *getEff_num_seq();

  void convertAseq2Alignment();

  subalign *sub2align(int *mark);
  subalign *sub2align(int *mark, int *mark1);

 public:
  char **aname, **aseq; /* aseq: [0,nal-1][0,alilen-1] */
  int **alignment;      /* alignment: [1,nal][1,alilen] */
  int nal, alilen, *astart, *alen;
  int mnamelen;

  // double *eff_num_seq;
  // int *ngap;

 private:
  void *mymalloc(int size);
  char *strsave(char *str);
  char *strnsave(char *str, int l);
  char **incbuf(int n, char **was);
  int *incibuf(int n, int *was);
  // convert the allocations by malloc to allocations by new
  void reassign();

  static int subaligncount;

 public:
  int *apos_filtr;
  double **n_effAa;
  double **pseudoCnt;
  double *sum_eff_let;
  int *maskgapRegion;
  // double qmatrix[21][21];
  double average_sum_eff_let;

  double *gap_content;

  int alilen_mat;
  double n_eff;
  double gap_threshold;
  double gapt_threshold;
  double beta;

  void profile();
  int done_profile;
  void neffsForEachCol_maskGapReg(int **ali, int n, int len, double effgapmax,
                                  double effgapRegionMin, double **n_effAa,
                                  double *sum_eff_let, int *maskgapRegion,
                                  int *apos_filtr, int *len_lowgaps,
                                  double *nef);
  double effective_number_nogaps(int **ali, int *marks, int n, int start,
                                 int end);
  void pseudoCounts(double **matrix, double n_eff, int len, double **pseudoCnt);

  void h_weight_all();
  double *hwt_all;

 public:
  double **distMat;

 public:
  // purge an alignment by
  // 1. removing highly similar sequences (id > high_id_thr)
  // 2. removing highly divergent sequences (id < low_id_thr)
  // 3. removing sequence fragment (gap_fraction > gap_fraction_thr)
  // 4. selecting only a subset if remaining num > max_num_kept
  subalign *purge_align(double low_id_thr, double high_id_thr, int max_num_kept,
                        double gap_fraction_thr);

  // my new version of profile - effective position determined by Henikoff gap
  // weight
 public:
  int *prof_pos;
  double **prof_effn;
  double **prof_freq;
  double *prof_sum_eff;
  double *prof_hwt_all;
  double *prof_gap_content;
  double prof_nef;
  double prof_average_sum_eff_let;
  int prof_len;
  double prof_raw_gap_threshold;
  double prof_gap_threshold;
  void prof_h_weight_all(double raw_gap_fraction_cutoff);
  void prof_positions(double prof_gap_threshold_here);
  void set_prof_raw_gap_threshold(double gapthr);
  void set_prof_gap_threshold(double gapthr);
  // only getting the effective counts portion
  void prof();
  void prof(double tmp_gap_threshold);
  void prof_get_effn(int **ali, int n, int len, int prof_len, double **n_effAa,
                     double *sum_eff_let, int *prof_pos, double *nef);
  int done_prof;
  int done_prof_freq;
  void pseudoCounts(double **matrix, double n_eff, int len, double **pseudoCnt,
                    float ***input_matrices, float **input_freqs,
                    int *mat_selections);
  void pseudoCounts(double **matrix, double n_eff, int len, double **pseudoCnt,
                    float **input_matrix, float *input_bfreq);
  void log_pseudoCounts();
  // get_freq, 0: blosum62; 1: general SCOP; 2: alphabet1-dependent SCOP
  void get_prof_freq(int get_freq, int take_log);

  // secondary structure profile
 public:
  ss_prof *ss;  // it is currently depending on only one representative sequence
  // 1. select representative
  char *repres_name;
  subalign *oneSeqAln;
  void select_representative();
  int select_representative_henikoff(float core_position_nongap_freq_cutoff);
  void get_oneSeqAln(int tmp_index);

  // 2. get secondary structure profile for the representative
  void get_ss_prof(char *dir_name, char *runpsipred_command);
  void get_ss_prof1(char *dir_name, char *query_name,
                    char *runpsipred1_command);

  // 3. map between profile positions and secondary structure profile positions
  int *prof_map_ss;  // mapping the profile positions to ss profile positions
  void get_prof_map_ss(char *rep_name);  // representative name

  // 4. map alphabet1 to profile positions
  int *prof_alphabet1;
  void get_prof_alphabet1();

  // 5. get alphabet1 dependent frequencies if needed
  // get_freq, 0: blosum62; 1: general SCOP; 2: alphabet1-dependent SCOP
  // void get_prof_freq(int get_freq);

 public:
  float *score_bg_aa;
  float *score_bg_ss;
  void get_score_bg(int bg_type);
  void get_score_bg_mine(float **aa_loop, int bg_type);
  void get_score_bg_mine(float *aa_loop0, float *ss_loop0, int bg_type);
  void get_score_bg_mine2(float *aa_loop0, float *ss_loop0, int bg_type);
  int done_get_score_bg;
  int bg_type_here;
  subalign *purge_align_one_seq_name(char *seqname);
  int done_score_bg;
  int done_score_bg2;
};

subalign *oneSeq2subalign(char *seq, char *name);

#endif
