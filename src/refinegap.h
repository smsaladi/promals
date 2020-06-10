#include "consv1.h"
#include "mathfunc.h"

void refinegap(subalign *aln, double gap_thr_lowest, int use_hwt,
               int remove_frag, int extend_gap_region);

void adjust_gap(subalign *aln, int start, int end);

void delete_complete_gap_positions(subalign *aln);

