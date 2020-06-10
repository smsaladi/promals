#include <algorithm>
#include <queue>
#include <stack>

#include "gap_refining.h"

static int debug_here = 1;

gap_refine::gap_refine() {
  num_end_letters = 5;
  min_cb_size = 4;
  problem_gr_size = 20;

  seq = 0;
  // cb = 0 ;
  // gr = 0;

  gr_begin = NULL;
  gr_end = NULL;
  done_get_pb62 = 0;
}

gap_refine::~gap_refine() {
  int i;

  if (seq) {
    delete[] seq;
  }
  if (gr_begin) delete[] gr_begin;
  if (gr_end) delete[] gr_end;
  gr_group.clear();
}

void gap_refine::setup_seq(char **aseq, int n, int len) {
  int i;
  nseqs = n;
  alilen = len;

  seq = new string[n];
  for (i = 0; i < n; i++) seq[i] = aseq[i];
}

void gap_refine::setup_seq(subalign *x) {
  setup_seq(x->aseq, x->nal, x->alilen);
}

void gap_refine::treat_end_letters() {
  int i, j, k, m, n;
  string gap1, gap2;
  string tmpstr1, tmpstr2;

  // cout << alilen << "\t" << seq[1].length() << endl;

  // 1. work with N-terminal discrete letters
  // string:      AAA---A---AAA
  // position:    0123456789012
  // variable:       j  km  n
  for (i = 0; i < nseqs; i++) {
    for (j = 0; j < alilen; j++) {
      if (seq[i][j] == '-') break;
    }
    if (j == alilen) continue;
    if (j <= num_end_letters) {
      tmpstr1 = string(seq[i], 0, j);
    } else
      continue;
    for (k = j; k < alilen; k++) {
      if (seq[i][k] != '-') break;
    }
    gap1 = seq[i].substr(j, k - j);
    for (m = k; m < alilen; m++) {
      if (seq[i][m] == '-') break;
    }
    tmpstr2 = seq[i].substr(k, m - k);
    for (n = m; n < alilen; n++) {
      if (seq[i][n] != '-') break;
    }
    gap2 = seq[i].substr(m, n - m);

    if (tmpstr1.length() + tmpstr2.length() <= num_end_letters) {
      seq[i] = gap1 + gap2 + tmpstr1 + tmpstr2 + seq[i].substr(n);
    } else {
      seq[i] = gap1 + tmpstr1 + seq[i].substr(k);
    }

    // cout << i << "  " << tmpstr1 << " " << gap1 << " " << tmpstr2 << " " <<
    // gap2 << endl;
  }

  // work with C-terminal discrete letters
  for (i = 0; i < nseqs; i++) {
    for (j = alilen - 1; j >= 0; j--) {
      if (seq[i][j] == '-') break;
    }
    if (j == 0) continue;
    if (alilen - 1 - j <= num_end_letters) {
      tmpstr1 = seq[i].substr(j + 1);
    } else
      continue;
    for (k = j; k >= 0; k--) {
      if (seq[i][k] != '-') break;
    }
    gap1 = seq[i].substr(k + 1, j - k);
    for (m = k; m >= 0; m--) {
      if (seq[i][m] == '-') break;
    }
    tmpstr2 = seq[i].substr(m + 1, k - m);
    for (n = m; n >= 0; n--) {
      if (seq[i][n] != '-') break;
    }
    gap2 = seq[i].substr(n + 1, m - n);

    if (tmpstr1.length() + tmpstr2.length() <= num_end_letters) {
      seq[i] = seq[i].substr(0, n + 1) + tmpstr2 + tmpstr1 + gap1 + gap2;
    } else {
      seq[i] = seq[i].substr(0, k + 1) + tmpstr1 + gap1;
    }
    // cout << i << "  " << tmpstr1 << " " << gap1 << " " << tmpstr2 << " " <<
    // gap2 << endl;
  }

  // 3. now change N, C-terminal gaps to small letters 'a'
  for (i = 0; i < nseqs; i++) {
    for (j = 0; j < alilen; j++) {
      if (seq[i][j] != '-') break;
      seq[i][j] = 'a';
    }

    for (j = alilen - 1; j >= 0; j--) {
      if (seq[i][j] != '-') break;
      seq[i][j] = 'a';
    }
  }
}

void gap_refine::define_cb_gr() {
  int i, j, k;

  int *cb_mark = new int[alilen];

  // 1. mark the positions
  stack<int> pstack;
  for (i = 0; i < alilen; i++) {
    for (j = 0; j < nseqs; j++) {
      if (seq[j][i] == '-') {
        cb_mark[i] = 0;
        // cout << i << " " << cb_mark[i] << endl;
        if ((i > 0) && (cb_mark[i - 1])) {
          // cout << "pstack: " << pstack.size() << " " << min_cb_size << endl;
          if (pstack.size() < min_cb_size) {
            while (!pstack.empty()) {
              k = pstack.top();
              cb_mark[k] = 0;
              pstack.pop();
            }
          }
        }
        break;
      }
    }
    if (j != nseqs) continue;
    if (cb_mark[i - 1] == 0) {
      while (!pstack.empty()) pstack.pop();
    }
    cb_mark[i] = 1;
    // cout << i << " " << cb_mark[i] << endl;
    pstack.push(i);
  }
  /*
  cout << endl;
  //for(i=0;i<alilen;i++) { cout << cb_mark[i]; }
  cout << endl;
  printseq(80, cb_mark);
  */
  // 2. find the boundaries
  num_gr = 0;
  for (i = 1; i < alilen; i++) {
    if (cb_mark[i] && (!cb_mark[i - 1])) {
      num_gr++;
    }
  }
  if (gr_begin) delete[] gr_begin;
  if (gr_end) delete[] gr_end;
  gr_begin = new int[num_gr];
  gr_end = new int[num_gr];

  int gr_count = 0;
  for (i = 1; i < alilen; i++) {
    if (cb_mark[i] && (!cb_mark[i - 1])) {
      gr_end[gr_count] = i;
      gr_count++;
    }
    if (cb_mark[i - 1] && (!cb_mark[i])) {
      gr_begin[gr_count] = i;
    }
  }
  for (i = 0; i < num_gr; i++) {
    if (debug_here > 1)
      fprintf(stdout, "%d\t%d\t%d\n", i, gr_begin[i], gr_end[i]);
  }
  delete[] cb_mark;
}

// treat a continuous gap in a gappy region
void gap_refine::treat_continuous_gaps() {
  int i, j, k;
  int cgap = 0;

  for (i = 0; i < num_gr; i++) {
    for (j = 0; j < nseqs; j++) {
      cgap = 0;
      if (debug_here > 1)
        cout << seq[j].substr(gr_begin[i], gr_end[i] - gr_begin[i]) << endl;
      for (k = gr_begin[i]; k < gr_end[i]; k++) {
        if (seq[j][k] != '-') {
          cgap = 1;
          break;
        }
      }
      if (!cgap) {
        // cout << "==========" << endl;
        for (k = gr_begin[i]; k < gr_end[i]; k++) {
          seq[j][k] = 'a';
        }
      }
    }
  }
  get_log_q_blosum62_ratio();
}

void gap_refine::do_refine(int mode) {
  int i, j, k, m, n;

  string *new_seq = new string[nseqs];
  int core_begin = 0, core_end;

  int *cb_mark;
  cb_mark = new int[seq[0].length() + 1];
  int tmpp = 0;

  for (i = 0; i < nseqs; i++) new_seq[i] = "";

  if (debug_here > 1) cout << "num_gr: " << num_gr << endl;
  if (num_gr == 0) {
    for (j = 0; j < nseqs; j++) {
      for (i = 0; i < seq[j].length(); i++) {
        if (seq[j][i] == 'a') seq[j][i] = '-';
      }
    }
    return;
  }
  for (i = 0; i < num_gr; i++) {
    // 1. for the core block
    if (i == 0)
      core_begin = 0;
    else
      core_begin = gr_end[i - 1];
    core_end = gr_begin[i];
    if (debug_here > 1)
      cout << "cores: " << core_begin << "  " << core_end << endl;
    for (j = 0; j < nseqs; j++) {
      new_seq[j] += seq[j].substr(core_begin, core_end - core_begin);
      if (debug_here > 1) cout << new_seq[j] << endl;
    }
    for (j = 0; j < core_end - core_begin; j++) {
      cb_mark[tmpp] = 1;
      tmpp++;
    }

    // 2. for the gappy region
    if (mode == 0) refine_group_push_aside(gr_begin[i], gr_end[i]);
    if (mode == 1) refine_group_by_seqnum(gr_begin[i], gr_end[i]);
    if (mode == 2) refine_group(gr_begin[i], gr_end[i]);
    for (k = 0; k < gr_group[0].gseq.size(); k++) {
      new_seq[gr_group[0].seqnum[k]] += gr_group[0].gseq[k];
    }
    for (j = 0; j < nseqs; j++) {
      if (debug_here > 1) cout << new_seq[j] << endl;
    }
    for (j = 0; j < gr_group[0].length; j++) {
      cb_mark[tmpp] = 0;
      tmpp++;
    }
  }
  core_begin = gr_end[i - 1];
  core_end = seq[0].length();
  if (debug_here > 1) cout << "nseqs: " << nseqs << endl;
  for (j = 0; j < nseqs; j++) {
    new_seq[j] += seq[j].substr(core_begin, core_end - core_begin);
    if (debug_here > 1) cout << "j: " << j << " " << new_seq[j] << endl;
  }
  for (j = 0; j < core_end - core_begin; j++) {
    cb_mark[tmpp] = 1;
    tmpp++;
  }

  for (j = 0; j < nseqs; j++) {
    for (i = 0; i < new_seq[j].length(); i++) {
      if (new_seq[j][i] == 'a') new_seq[j][i] = '-';
    }
  }
  if (debug_here > 1) cout << endl;
  for (j = 0; j < nseqs; j++) {
    seq[j] = new_seq[j];
    if (debug_here > 1) cout << new_seq[j] << endl;
  }
  if (debug_here > 1) cout << endl << endl;
  alilen = seq[0].length();

  if (debug_here > 1) printseq(60, cb_mark);

  delete[] new_seq;
  delete[] cb_mark;
}

void gap_refine::get_pb62() {
  int i, j;
  double minvalue = 100;

  if (done_get_pb62) return;

  double bg_freq_gap[21];
  for (i = 1; i <= 20; i++) {
    bg_freq_gap[i] = 0;
    for (j = 1; j <= 20; j++) {
      bg_freq_gap[i] += gap_matrix[i][j];
    }
    if (debug_here > 1) cout << "i: " << i << " " << bg_freq_gap[i] << endl;
  }
  for (i = 1; i <= 20; i++) {
    for (j = 1; j <= 20; j++) {
      pb62[i][j] = log(gap_matrix[i][j] / bg_freq_gap[i] / bg_freq_gap[j]);
    }
  }

  for (i = 1; i <= 20; i++) {
    for (j = 1; j <= 20; j++) {
      if (pb62[i][j] < minvalue) {
        minvalue = pb62[i][j];
      }
    }
  }
  for (i = 1; i <= 20; i++) {
    for (j = 1; j <= 20; j++) {
      pb62[i][j] = pb62[i][j] - minvalue;
      pb62[i][j] += 0.001;
    }
  }

  for (i = 1; i <= 20; i++) {
    for (j = 1; j <= 20; j++) {
      if (debug_here > 1) cout << pb62[i][j] << " ";
    }
    if (debug_here > 1) cout << endl;
  }

  done_get_pb62 = 1;
}

void gap_refine::printseq(int blocksize) {
  int i, j, k;

  int blocknum = alilen / blocksize;

  for (i = 0; i < blocknum; i++) {
    for (j = 0; j < nseqs; j++) {
      cout << seq[j].substr(i * blocksize, blocksize);
      cout << endl;
    }
    cout << endl;
  }

  if (alilen == blocknum * blocksize) return;

  for (j = 0; j < nseqs; j++) {
    cout << seq[j].substr(blocknum * blocksize, alilen);
    cout << endl;
  }
  cout << endl;
}

void gap_refine::printseq(int blocksize, int *mark) {
  int i, j, k;

  int blocknum = alilen / blocksize;

  for (i = 0; i < blocknum; i++) {
    for (j = i * blocksize; j < i * blocksize + blocksize; j++) {
      cout << mark[j];
    }
    cout << endl;
    for (j = 0; j < nseqs; j++) {
      cout << seq[j].substr(i * blocksize, blocksize);
      cout << endl;
    }
    cout << endl;
  }

  if (alilen == blocknum * blocksize) return;

  for (i = blocknum * blocksize; i < alilen; i++) {
    cout << mark[i];
  }
  cout << endl;
  for (j = 0; j < nseqs; j++) {
    cout << seq[j].substr(blocknum * blocksize, alilen);
    cout << endl;
  }
  cout << endl;
}

// alignment order of the groups by the sequence length
// align by substitution matrix
void gap_refine::refine_group(int begin, int end) {
  int i, j, k;

  gr_group.clear();

  string tmpstr;

  int new_group = 0;
  if (debug_here > 1) cout << begin << "\t" << end << endl;

  // 1. define groups
  for (i = 0; i < nseqs; i++) {
    tmpstr = "";
    new_group = 1;
    for (j = begin; j < end; j++) {
      if ((seq[i][j] != '-') and (seq[i][j] != 'a')) {
        tmpstr += seq[i][j];
      }
    }
    if (debug_here > 1) cout << i << ": " << tmpstr << endl;
    for (j = 0; j < gr_group.size(); j++) {
      if (tmpstr.size() == gr_group[j].length) {
        gr_group[j].gseq.push_back(tmpstr);
        gr_group[j].seqnum.push_back(i);
        new_group = 0;
        break;
      }
    }
    if (!new_group) continue;
    group ngr;
    ngr.length = tmpstr.size();
    ngr.seqnum.push_back(i);
    ngr.gseq.push_back(tmpstr);
    if (gr_group.size() == 0) {
      gr_group.push_back(ngr);
      continue;
    }
    vector<group>::iterator tmpiter;
    for (tmpiter = gr_group.begin(); tmpiter != gr_group.end(); tmpiter++) {
      if (ngr.length > tmpiter->length) {
        break;
      }
    }
    gr_group.insert(tmpiter, ngr);
  }

  // print the groups
  for (i = 0; i < gr_group.size(); i++) {
    if (debug_here > 1) cout << "group: " << i << endl;
    for (j = 0; j < gr_group[i].gseq.size(); j++) {
      if (debug_here > 1)
        cout << gr_group[i].seqnum[j] << " " << gr_group[i].gseq[j] << endl;
    }
    if (debug_here > 1) cout << endl;
  }

  // 2. align every group to the first reference group with the longest length
  for (i = 1; i < gr_group.size(); i++) {
    gap_align_two_groups(gr_group[0], gr_group[i]);
  }
}

// just split the residues and push them aside
void gap_refine::refine_group_push_aside(int begin, int end) {
  int i, j, k;

  gr_group.clear();

  // first find the length of longest sequence
  group ngr;
  ngr.length = 0;
  string tmpstr;
  for (i = 0; i < nseqs; i++) {
    tmpstr = "";
    for (j = begin; j < end; j++) {
      if ((seq[i][j] != '-') and (seq[i][j] != 'a')) {
        tmpstr += seq[i][j];
      }
    }
    ngr.gseq.push_back(tmpstr);
    ngr.seqnum.push_back(i);
    if (ngr.length < ngr.gseq[i].length()) {
      ngr.length = tmpstr.length();
    }
    tmpstr.clear();
  }
  for (i = 0; i < nseqs; i++) {
    tmpstr = "";
    tmpstr.insert(0, ngr.length - ngr.gseq[i].length(), '-');
    ngr.gseq[i].insert(ngr.gseq[i].length() / 2, tmpstr);
  }
  gr_group.push_back(ngr);

  for (i = 0; i < gr_group[0].gseq.size(); i++) {
    if (debug_here > 1) cout << "gseq: " << gr_group[0].gseq[i] << endl;
  }
}

struct comparison_group {
  bool operator()(const group &a, const group &b) {
    if (a.gseq.size() > b.gseq.size()) {
      return true;
    } else if (a.gseq.size() < b.gseq.size()) {
      return false;
    } else {
      if (a.length > b.length)
        return false;
      else
        return true;
    }
  }
};

// alignment order of the groups by the number of sequences
void gap_refine::refine_group_by_seqnum(int begin, int end) {
  int i, j, k;

  gr_group.clear();

  string tmpstr;

  int new_group = 0;
  if (debug_here > 1) cout << begin << "\t" << end << endl;

  // 1. define groups
  for (i = 0; i < nseqs; i++) {
    tmpstr = "";
    new_group = 1;
    for (j = begin; j < end; j++) {
      if ((seq[i][j] != '-') and (seq[i][j] != 'a')) {
        tmpstr += seq[i][j];
      }
    }
    if (debug_here > 1) cout << i << ": " << tmpstr << endl;
    for (j = 0; j < gr_group.size(); j++) {
      if (tmpstr.size() == gr_group[j].length) {
        gr_group[j].gseq.push_back(tmpstr);
        gr_group[j].seqnum.push_back(i);
        new_group = 0;
        break;
      }
    }
    if (!new_group) continue;
    group ngr;
    ngr.length = tmpstr.size();
    ngr.seqnum.push_back(i);
    ngr.gseq.push_back(tmpstr);
    if (gr_group.size() == 0) {
      gr_group.push_back(ngr);
      continue;
    }
    /*
    vector<group>::iterator tmpiter;
    for(tmpiter=gr_group.begin();tmpiter!=gr_group.end();tmpiter++) {
            // Here is the difference from refing_group
            if(ngr.gseq.size()>tmpiter->gseq.size()) {
                    break;
            }
            else if(ngr.gseq.size()==tmpiter->gseq.size()) {
                    if(ngr.length>tmpiter->length) break;
            }
    }
    */
    gr_group.insert(gr_group.begin(), ngr);
  }

  // sort the groups by 1: number of sequences 2: sequence length
  sort(gr_group.begin(), gr_group.end(), comparison_group());

  // print the groups
  for (i = 0; i < gr_group.size(); i++) {
    if (debug_here > 1) cout << "group by seqnum: " << i << endl;
    for (j = 0; j < gr_group[i].gseq.size(); j++) {
      if (debug_here > 1)
        cout << gr_group[i].seqnum[j] << " " << gr_group[i].gseq[j] << endl;
    }
    if (debug_here > 1) cout << endl;
  }

  // 2. align every group to the first reference group with the longest length
  for (i = 1; i < gr_group.size(); i++) {
    if (gr_group[0].length > gr_group[i].length)
      gap_align_two_groups(gr_group[0], gr_group[i]);
    else
      gap_align_two_groups(gr_group[i], gr_group[0]);
  }
}

void gap_refine::gap_align_two_groups(group &a, group &b) {
  int i, j, k, m, n;

  assert(a.length > b.length);

  // cout << a.length << " " << b.length << endl;

  double total_score, max_score;
  int max_pos;

  // 1. calculate initial score
  //  group a:    AAAAAAAAAA
  //  group b:    ------AAAA
  total_score = 0;
  queue<double> scores;
  double pscore;
  for (i = a.length - b.length; i < a.length; i++) {
    m = i - a.length + b.length;
    pscore = 0;
    for (j = 0; j < a.gseq.size(); j++) {
      if ((a.gseq[j][i] == '-') || (a.gseq[j][i] == 'a')) continue;
      for (k = 0; k < b.gseq.size(); k++) {
        if ((b.gseq[k][m] == '-') || (b.gseq[k][m] == 'a')) continue;
        pscore += pb62[am2num(a.gseq[j][i])][am2num(b.gseq[k][m])];
      }
    }
    if (debug_here > 1) cout << i << "\t" << pscore << endl;
    total_score += pscore;
    scores.push(pscore);
  }

  if (debug_here > 1) cout << "i: 0 " << total_score << endl;

  // 2. find the max score and position
  // group a: AAAAAAAAAA   group a: AAAAAAAAAA    group a: AAAAAAAAAA  group a:
  // AAAAAAAAAA group b: A------AAA   group b: AA------AA    group b: AAA------A
  // group b: AAAA------
  max_score = total_score;
  max_pos = 0;
  for (i = 0; i < b.length; i++) {
    total_score -= scores.front();
    if (debug_here > 1) cout << "i: " << i + 1 << " " << scores.front() << endl;
    scores.pop();
    pscore = 0;
    for (j = 0; j < a.gseq.size(); j++) {
      if ((a.gseq[j][i] == '-') || (a.gseq[j][i] == 'a')) continue;
      for (k = 0; k < b.gseq.size(); k++) {
        if ((b.gseq[k][i] == '-') || (b.gseq[k][i] == 'a')) continue;
        pscore += pb62[am2num(a.gseq[j][i])][am2num(b.gseq[k][i])];
      }
    }
    total_score += pscore;
    if (debug_here > 1)
      cout << "i: " << i + 1 << " " << total_score << "\t" << pscore << endl;
    if (max_score < total_score) {
      max_score = total_score;
      max_pos = i + 1;
    }
  }

  // 3. add aligned group b to group a
  // also add group a to group b
  string tmpstring = "";
  string gapstring = "";
  for (j = 0; j < a.length - b.length; j++) gapstring += '-';
  if (debug_here > 1) cout << gapstring << endl;
  if (debug_here > 1)
    cout << "a, b sizes before: " << a.gseq.size() << " " << b.gseq.size()
         << endl;
  for (i = 0; i < b.gseq.size(); i++) {
    b.gseq[i].insert(max_pos, gapstring);
    tmpstring = b.gseq[i];
    if (debug_here > 1) cout << "tmpstring b: " << tmpstring << endl;
    a.gseq.push_back(tmpstring);
    a.seqnum.push_back(b.seqnum[i]);
  }
  b.length = a.length;
  if (debug_here > 1)
    cout << "a, b sizes middle: " << a.gseq.size() << " " << b.gseq.size()
         << endl;
  int size_diff = a.gseq.size() - b.gseq.size();
  for (i = 0; i < size_diff; i++) {
    b.gseq.push_back(a.gseq[i]);
    if (debug_here > 1) cout << "tmpstring a: " << a.gseq[i] << endl;
    b.seqnum.push_back(a.seqnum[i]);
  }
  if (debug_here > 1)
    cout << "a, b sizes after: " << a.gseq.size() << " " << b.gseq.size()
         << endl;
}

char **gap_refine::outalign() {
  int i, j;
  char **tmpseqs = cmatrix(nseqs, seq[0].length() + 1);

  for (i = 0; i < nseqs; i++) {
    strcpy(tmpseqs[i], seq[i].c_str());
  }

  return tmpseqs;
}

char **gap_refine::batch_refine(subalign *x, int mode) {
  int i, j;

  // set up the sequences
  setup_seq(x);

  // 1. treat the N, C-terminal gaps
  treat_end_letters();

  // 2. define cb and gr
  define_cb_gr();

  // 3. treat long continues gaps
  treat_continuous_gaps();

  // 4. re-define cb and gr
  define_cb_gr();

  // 5. get_substitution matrix
  get_pb62();

  // 6. make refinement
  do_refine(mode);

  // 7. return the new alignment sequences
  char **tmpseq = outalign();
  return tmpseq;
}

double gap_matrix[21][21] = {
    {0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
     0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
     0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
     0.0000000, 0.0000000, 0.0000000},
    {0.0000000, 0.0039102, 0.0009895, 0.0009106, 0.0002025, 0.0007992,
     0.0003218, 0.0003036, 0.0004195, 0.0000562, 0.0003790, 0.0003136,
     0.0002894, 0.0003746, 0.0002176, 0.0002015, 0.0002529, 0.0003380,
     0.0002490, 0.0002850, 0.0002441},
    {0.0000000, 0.0009806, 0.0093634, 0.0030962, 0.0008752, 0.0038506,
     0.0017161, 0.0017384, 0.0013605, 0.0002646, 0.0011010, 0.0009473,
     0.0009255, 0.0009333, 0.0007588, 0.0005395, 0.0006596, 0.0007912,
     0.0006950, 0.0006817, 0.0007330},
    {0.0000000, 0.0009153, 0.0030620, 0.0083944, 0.0005381, 0.0019937,
     0.0010081, 0.0011893, 0.0012907, 0.0002615, 0.0011547, 0.0008818,
     0.0008850, 0.0011598, 0.0010253, 0.0006292, 0.0010629, 0.0008921,
     0.0010811, 0.0009122, 0.0009383},
    {0.0000000, 0.0002155, 0.0008923, 0.0005241, 0.0027075, 0.0026945,
     0.0013123, 0.0011412, 0.0010719, 0.0001984, 0.0007109, 0.0005955,
     0.0007601, 0.0007725, 0.0005126, 0.0005459, 0.0005859, 0.0005765,
     0.0002948, 0.0006433, 0.0007173},
    {0.0000000, 0.0007611, 0.0037210, 0.0020004, 0.0025958, 0.0215940,
     0.0058807, 0.0053059, 0.0034445, 0.0006311, 0.0020432, 0.0024349,
     0.0026034, 0.0020798, 0.0014238, 0.0014475, 0.0014411, 0.0017492,
     0.0010517, 0.0018499, 0.0020733},
    {0.0000000, 0.0003171, 0.0016034, 0.0009634, 0.0011657, 0.0055896,
     0.0074624, 0.0053806, 0.0018831, 0.0003600, 0.0008793, 0.0013952,
     0.0015445, 0.0010606, 0.0007074, 0.0006656, 0.0007195, 0.0009276,
     0.0005066, 0.0008821, 0.0009174},
    {0.0000000, 0.0002996, 0.0015607, 0.0011659, 0.0010785, 0.0050183,
     0.0054465, 0.0102301, 0.0036039, 0.0005863, 0.0014553, 0.0021405,
     0.0026700, 0.0016971, 0.0010124, 0.0009070, 0.0011967, 0.0014877,
     0.0007006, 0.0012821, 0.0015561},
    {0.0000000, 0.0003785, 0.0013177, 0.0013388, 0.0011030, 0.0033671,
     0.0019645, 0.0038093, 0.0154740, 0.0011544, 0.0057505, 0.0040874,
     0.0035110, 0.0054485, 0.0024492, 0.0018763, 0.0031292, 0.0034138,
     0.0010006, 0.0023390, 0.0031318},
    {0.0000000, 0.0000751, 0.0002631, 0.0002744, 0.0002053, 0.0006754,
     0.0003892, 0.0006560, 0.0011539, 0.0056248, 0.0006675, 0.0004046,
     0.0006189, 0.0006993, 0.0003631, 0.0001968, 0.0004223, 0.0002427,
     0.0002024, 0.0002540, 0.0003221},
    {0.0000000, 0.0004538, 0.0012553, 0.0013713, 0.0007802, 0.0024373,
     0.0011417, 0.0018999, 0.0062077, 0.0006938, 0.0643787, 0.0036281,
     0.0034441, 0.0059827, 0.0058053, 0.0024807, 0.0059180, 0.0043755,
     0.0016898, 0.0030195, 0.0043169},
    {0.0000000, 0.0004031, 0.0011191, 0.0010898, 0.0006935, 0.0028569,
     0.0017141, 0.0026543, 0.0045373, 0.0004610, 0.0036157, 0.0318834,
     0.0030144, 0.0038882, 0.0021975, 0.0018043, 0.0035972, 0.0036702,
     0.0010193, 0.0023853, 0.0035844},
    {0.0000000, 0.0003469, 0.0009863, 0.0009368, 0.0008234, 0.0027450,
     0.0017033, 0.0029373, 0.0036905, 0.0006125, 0.0032197, 0.0027034,
     0.0140884, 0.0067152, 0.0031120, 0.0017477, 0.0034164, 0.0028614,
     0.0009885, 0.0021706, 0.0030688},
    {0.0000000, 0.0004095, 0.0010500, 0.0012813, 0.0008041, 0.0023311,
     0.0012853, 0.0019296, 0.0058520, 0.0007498, 0.0057888, 0.0036767,
     0.0068276, 0.0162028, 0.0043574, 0.0022513, 0.0052003, 0.0039520,
     0.0014514, 0.0028768, 0.0035961},
    {0.0000000, 0.0002603, 0.0008230, 0.0011902, 0.0005569, 0.0016867,
     0.0008442, 0.0013436, 0.0027315, 0.0004054, 0.0058527, 0.0021284,
     0.0032685, 0.0043414, 0.0142724, 0.0021682, 0.0063034, 0.0032330,
     0.0016839, 0.0023615, 0.0035282},
    {0.0000000, 0.0002150, 0.0005555, 0.0006852, 0.0005245, 0.0014762,
     0.0007176, 0.0010218, 0.0019436, 0.0001926, 0.0022167, 0.0015456,
     0.0016943, 0.0020878, 0.0019345, 0.0044682, 0.0022127, 0.0029511,
     0.0010001, 0.0019546, 0.0027177},
    {0.0000000, 0.0002792, 0.0007523, 0.0011546, 0.0006195, 0.0017459,
     0.0008638, 0.0014479, 0.0033505, 0.0004177, 0.0058217, 0.0034438,
     0.0035901, 0.0052527, 0.0062217, 0.0023651, 0.0275551, 0.0068429,
     0.0014860, 0.0024433, 0.0037737},
    {0.0000000, 0.0003375, 0.0008045, 0.0009397, 0.0005413, 0.0017689,
     0.0009625, 0.0016448, 0.0033382, 0.0002328, 0.0038283, 0.0031835,
     0.0027152, 0.0036591, 0.0029089, 0.0029062, 0.0062593, 0.0129229,
     0.0012254, 0.0025687, 0.0042063},
    {0.0000000, 0.0002651, 0.0007972, 0.0010679, 0.0003276, 0.0011707,
     0.0005536, 0.0007539, 0.0010863, 0.0001992, 0.0015908, 0.0008733,
     0.0009925, 0.0014060, 0.0016196, 0.0010111, 0.0014525, 0.0013132,
     0.0071331, 0.0011943, 0.0012950},
    {0.0000000, 0.0002985, 0.0006650, 0.0009733, 0.0006493, 0.0018883,
     0.0009427, 0.0014059, 0.0023518, 0.0002720, 0.0027567, 0.0021077,
     0.0020494, 0.0026649, 0.0021965, 0.0019347, 0.0022274, 0.0025055,
     0.0011810, 0.0122971, 0.0048944},
    {0.0000000, 0.0002371, 0.0007228, 0.0009528, 0.0007126, 0.0021754,
     0.0010560, 0.0016327, 0.0031605, 0.0002890, 0.0038439, 0.0031339,
     0.0029675, 0.0033850, 0.0032651, 0.0027362, 0.0035221, 0.0042092,
     0.0011946, 0.0049149, 0.0120942}};
