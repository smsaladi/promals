#include "time.h"
#include "multiple.h"

int main(int argc, char **argv) {
  time(&timestart);
  timecount = 0;

  getParameter(argc, argv, 1);
  printParameters();
  printParameters1();

  get_log_robinson_freq();
  get_log_q_blosum62_ratio();
  get_log_q_blosum62();

  if (probconsBLOSUM) useProbconsFrequencies();

  multiple *my_align = new multiple(argv[1]);

  if (my_align->allseqs.nseqs <= my_align->max_cluster_elem) {
    my_align->build_tree();
    my_align->set_distance_cutoff_similar(1 - id_thr, max_group_number);
    my_align->alignSimilar();
  } else {
    my_align->build_cluster();
  }

  // remove filter files
  char command[500];
  sprintf(command, "rm -f %s_fiLtEr*", argv[1]);
  system(command);

  my_align->alignDivergent_psipred(2);
}

