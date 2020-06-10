#include <map>
#include <vector>
#include <algorithm>

#include "amino.h"
#include "consv1.h"
#include "header_cpp.h"
#include "hmm_local.h"
#include "hmm_profpair.h"
#include "kmer_dist.h"
#include "mathfunc.h"
#include "sequences.h"
#include "smallnumber.h"
#include "subalign.h"
#include "tnode.h"
#include "util.h"
#include "hmm_multim.h"
#include "mm.h"
#include "param.h"
#include "refinegap.h"
#include "sparsematrix.h"
#include "ss_prof.h"
#include "time.h"
#include "blastpgp.h"
#include "hmm_psipred.h"
#include "profilehmm.h"
#include "refine.h"
#include "regularizer.h"
#include "seq_str_aln.h"

#include "ScoreType.h"

using namespace std;

extern struct tms tmsstart, tmsend;

