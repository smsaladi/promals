#include "header_cpp.h"
#include "amino.h"
#include "consv1.h"
#include "hmm_profpair.h"
#include "hmm_profpair1.h"
#include "hmm_local.h"
#include "kmer_dist.h"
#include "util.h"
#include "mathfunc.h"
#include "smallnumber.h"
#include "subalign.h"
#include "sequences.h"
#include "tnode.h"
//#include "btree_template.h"
//#include "progressiveAlignHMM.h"
#include "sparsematrix.h"
#include "time.h"
#include "mm.h"
#include "ScoreType.h"
#include "param.h"
#include "hmm_multim.h"
#include <vector>
#include <map>
#include "refinegap.h"
#include "ss_prof.h"
//#include "multiple.h"
#include "profilehmm.h"
#include "regularizer.h"
#include "blastpgp.h"
#include "hmm_psipred.h"
#include <algorithm>
#include "refine.h"
#include "seq_str_aln.h"
//#include "constraint.h"

using namespace std;

extern struct tms      tmsstart, tmsend;
