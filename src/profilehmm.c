#include "all.h"
#include "profilehmm.h"

static int debug = 0;
profilehmm::profilehmm(subalign *x1) {
  x = x1;
  // lenx = x->alilen;

  nal = x->nal;
  alilen = x->alilen;

  use_position_specific_regularizer = 0;
  reg_weight = 1;

  // gapopen_end = 0;

  y = 0;

  ss_weight = 0.5;
  use_picasso = 0;
}

/*
// set up HMMER regularizer
void profilehmm::set_regularizer() {

        int i;

        reg_m[1] = 0.7939;
        reg_m[2] = 0.0278;
        reg_m[3] = 0.0135;
        reg_i[1] = 0.1551;
        reg_i[2] = 0.1331;
        reg_i[3] = 0.001;  // added
        reg_d[1] = 0.9002;
        reg_d[3] = 0.563;
        reg_d[2] = 0.001;   // added

        sreg_m = NULL;
        sreg_i = NULL;
        sreg_d = NULL;
        sreg_weight = 1;
        use_position_specific_regularizer = 0;
}
*/

void profilehmm::set_align(subalign *y1) {
  y = y1;

  if (y->done_prof)
    leny = y->prof_len;
  else
    leny = 0;

  // Vm=Vx=Vy=0; Tm=Tx=Ty=0;
  // path=path1=path2=0;
  // Fm=Fx=Fy=Bm=Bx=By=0;
  // probMat = 0;
}

profilehmm::~profilehmm() {
  free_gmatrix<float>(Fm, lenx, leny);
  free_gmatrix<float>(Fi, lenx, leny);
  free_gmatrix<float>(Fd, lenx, leny);
  free_gvector<float>(Fn);
  free_gvector<float>(Fb);
  free_gvector<float>(Fe);
  free_gvector<float>(Fc);
  free_gmatrix<float>(Bm, lenx, leny);
  free_gmatrix<float>(Bi, lenx, leny);
  free_gmatrix<float>(Bd, lenx, leny);
  free_gvector<float>(Bn);
  free_gvector<float>(Bb);
  free_gvector<float>(Be);
  free_gvector<float>(Bc);

  free_gmatrix<float>(score_matrix, lenx, leny);
  free_gvector<float>(score_bg);

  free_gmatrix<float>(probMat, lenx, leny);

  free_gmatrix<float>(t_m, lenx, 4);
  free_gmatrix<float>(t_i, lenx, 3);
  free_gmatrix<float>(t_i, lenx, 3);
  free_gvector<float>(t_bm);
}

void profilehmm::set_parameters(char *file_name, char *ss_dir_name,
                                int use_ss) {
  int i, j, k;

  // read the paramter file
  if (!done_read_regularizer) read_regularizer(file_name);

  // calculate the profile of x
  // including information about frequencies and counts, independent gap
  // contents default effective gap threshold is 0.5
  if (!x->done_prof) x->prof();
  if ((use_ss == 1) || (use_ss == 0)) {
    if (!x->done_prof) x->get_prof_freq(use_ss, 0);
  } else if (use_ss == 2) {
    if (!x->ss) {
      x->select_representative();
      x->get_ss_prof(ss_dir_name, runpsipred_command);
      // cout << x->repres_name << endl;
      x->get_prof_map_ss(x->repres_name);
      x->get_prof_alphabet1();
    }
    if (!x->done_prof) x->get_prof_freq(use_ss, 0);
    use_position_specific_regularizer = 1;
  } else {
    cout << "Error: use_ss must be 0, 1, or 2" << endl;
    exit(0);
  }

  nal = x->nal;
  alilen = x->alilen;
  alignment = x->alignment;
  // apos_filtr = x->apos_filtr;
  prof_len = x->prof_len;
  prof_pos = x->prof_pos;
  // ssx = x->ss;
  hwt_all = x->prof_hwt_all;
  prof_sum_eff = x->prof_sum_eff;

  lenx = x->prof_len;

  // transform sequences into sequential numbers
  // sequence example: -AAC---DEF--G-
  // number sign     : -+++---+++--+-
  // seq number      : 01233334566677
  int **seq_num = imatrix(nal, alilen);
  for (i = 1; i <= x->nal; i++) {
    int tmp_aa_count = 0;
    for (j = 1; j <= alilen; j++) {
      if (alignment[i][j] != 0) {
        tmp_aa_count += 1;
        seq_num[i][j] = tmp_aa_count;
      } else {
        seq_num[i][j] = 0 - tmp_aa_count;
      }
    }
  }

  // estimate transition probabilities
  // float t_sn, t_sb, t_nn, t_nb;
  // float *t_bm, t_bd1;
  // float **t_m, **t_i; **t_d;
  // float t_ec, t_et, t_ct, t_cc;
  t_bm = gvector<float>(lenx);
  t_m = gmatrix<float>(lenx, 4);  // 1: m;  2: i; 3: d; 4: e
  t_i = gmatrix<float>(lenx, 3);  // 1: m;  2: i; 3: d;
  t_d = gmatrix<float>(lenx, 3);  // 1: m;  2: i; 3: d;

  for (i = 1; i <= lenx; i++) {
    t_bm[i] = 0;
    for (j = 1; j <= 3; j++) {
      t_m[i][j] = 0;
      t_i[i][j] = 0;
      t_d[i][j] = 0;
    }
    t_m[i][4] = 0;
  }

  // if(debug>-1) for(i=1;i<=x->nal;i++) { cout << i << " " << hwt_all[i] <<
  // endl; }

  int pos_diff;

  cout << "t_i: " << t_i[6][1] << " " << t_i[6][2] << " " << t_i[6][3] << endl;
  for (i = 1; i <= x->nal; i++) {
    // from 1 to prof_len-1
    for (j = 1; j < prof_len; j++) {
      // both are letters
      if ((seq_num[i][prof_pos[j]] > 0) && (seq_num[i][prof_pos[j + 1]] > 0)) {
        pos_diff = seq_num[i][prof_pos[j + 1]] - seq_num[i][prof_pos[j]];
        if (pos_diff == 1) {  // no insertion(s) in between
          t_m[j][1] += hwt_all[i];
        } else {  // insertions in between
          t_m[j][2] += hwt_all[i];
          t_i[j][1] += hwt_all[i];
          t_i[j][2] += hwt_all[i] * (pos_diff - 2);
        }
      }
      // first is letter, second is gap (D)
      if ((seq_num[i][prof_pos[j]] > 0) && (seq_num[i][prof_pos[j + 1]] <= 0)) {
        pos_diff = 0 - seq_num[i][prof_pos[j + 1]] - seq_num[i][prof_pos[j]];
        if (pos_diff == 0) {
          t_m[j][3] += hwt_all[i];
        } else {
          t_m[j][2] += hwt_all[i];
          t_i[j][3] += hwt_all[i];
          t_i[j][2] += hwt_all[i] * (pos_diff - 1);
        }
      }
      // first is gap, second is letter
      if ((seq_num[i][prof_pos[j]] <= 0) && (seq_num[i][prof_pos[j + 1]] > 0)) {
        pos_diff = seq_num[i][prof_pos[j + 1]] + seq_num[i][prof_pos[j]];
        if (pos_diff == 1) {
          t_d[j][1] += hwt_all[i];
        } else {
          t_d[j][2] += hwt_all[i];
          t_i[j][1] += hwt_all[i];
          t_i[j][2] += hwt_all[i] * (pos_diff - 2);
        }
      }
      // first is gap, second is gap
      if ((seq_num[i][prof_pos[j]] <= 0) &&
          (seq_num[i][prof_pos[j + 1]] <= 0)) {
        pos_diff = seq_num[i][prof_pos[j]] - seq_num[i][prof_pos[j + 1]];
        if (pos_diff == 0) {
          t_d[j][3] += hwt_all[i];
        } else {
          t_d[j][2] += hwt_all[i];
          t_i[j][3] += hwt_all[i];
          t_i[j][2] += hwt_all[i] * (pos_diff - 1);
        }
      }
    }

    /* wing retraction disables these transition estimations
    // for the transition to the first position
    if(seq_num[i][prof_pos[1]]==0) t_bd1 += hwt_all[i];
    else t_bm[1] += hwt_all[i];

    // for the transition from the last position to th ending position
    if(seq_num[i][prof_pos[prof_len]]==0) t_me += hwt_all[i];
    else t_de += hwt_all[i];
    */
  }

  free_imatrix(seq_num, nal, alilen);

  // transitions from the start and the N state
  t_sb = rsb;
  t_sn = rsn;
  t_nn = rnn;
  t_nb = rnb;
  if (debug > 1)
    cout << t_sb << " " << t_sn << " " << t_nn << " " << t_nb << endl;

  // transitions from the B state to the m states
  // here adopt the wing retraction
  t_bd1 = 0;
  if (prof_len <= 20) {
    t_bm[1] = 0.5;
    for (i = 2; i <= prof_len; i++) {
      t_bm[i] = 0.5 / (prof_len - 1);
    }
  } else {
    for (i = 1; i <= 10; i++) {
      t_bm[i] = rbm[i];
      if (debug > 1) cout << "t_bm: " << i << " " << t_bm[i] << endl;
    }
    for (i = 11; i <= prof_len; i++) {
      t_bm[i] = rbm[11] / (prof_len - 10);
      if (debug > 1) cout << "t_bm: " << i << " " << t_bm[i] << endl;
    }
  }

  // transitions from the m states to the E state
  t_m[lenx][4] = 1;
  for (i = 1; i < prof_len; i++) {
    t_m[i][4] = 0.5 / (prof_len - 1);
    if (debug > 1) cout << "t_m4: " << i << " " << t_m[i][4] << endl;
  }

  // transitions from the E state
  // transitions from the C state
  t_ec = rec;
  t_et = ret;
  t_cc = rcc;
  t_ct = rct;

  // transitions among the M, I and D states
  double count_mm, count_mi, count_md, count_im, count_ii, count_id, count_dm,
      count_di, count_dd;
  double total_count;

  if (!use_position_specific_regularizer) {
    for (i = 1; i < prof_len; i++) {
      if (debug > -1) cout << i << " prof_sum_eff " << prof_sum_eff[i] << endl;
      if (debug > -1)
        cout << i << " m " << t_m[i][1] << " " << t_m[i][2] << " " << t_m[i][3]
             << " " << t_m[i][4] << endl;
      if (debug > -1)
        cout << prof_sum_eff[i] * t_m[i][1] << " " << rtrans0[1][1] * reg_weight
             << endl;
      if (debug > -1)
        cout << prof_sum_eff[i] * t_m[i][2] << " " << rtrans0[1][2] * reg_weight
             << endl;
      if (debug > -1)
        cout << prof_sum_eff[i] * t_m[i][3] << " " << rtrans0[1][3] * reg_weight
             << endl;
      count_mm = prof_sum_eff[i] * t_m[i][1] + rtrans0[1][1] * reg_weight;
      count_mi = prof_sum_eff[i] * t_m[i][2] + rtrans0[1][2] * reg_weight;
      count_md = prof_sum_eff[i] * t_m[i][3] + rtrans0[1][3] * reg_weight;
      total_count = count_mm + count_mi + count_md;
      t_m[i][1] = (1 - t_m[i][4]) * count_mm / total_count;
      t_m[i][2] = (1 - t_m[i][4]) * count_mi / total_count;
      t_m[i][3] = (1 - t_m[i][4]) * count_md / total_count;
      if (debug > -1)
        cout << i << " m " << t_m[i][1] << " " << t_m[i][2] << " " << t_m[i][3]
             << " " << t_m[i][4] << endl;

      if (i > 1) {
        if (debug > -1)
          cout << i << " d " << t_d[i][1] << " " << t_d[i][2] << " "
               << t_d[i][3] << endl;
        count_dm = prof_sum_eff[i] * t_d[i][1] + rtrans0[3][1] * reg_weight;
        count_di = prof_sum_eff[i] * t_d[i][2] + rtrans0[3][2] * reg_weight;
        count_dd = prof_sum_eff[i] * t_d[i][3] + rtrans0[3][3] * reg_weight;
        total_count = count_dm + count_di + count_dd;
        t_d[i][1] = count_dm / total_count;
        t_d[i][2] = count_di / total_count;
        t_d[i][3] = count_dd / total_count;
        if (debug > -1)
          cout << i << " d " << t_d[i][1] << " " << t_d[i][2] << " "
               << t_d[i][3] << endl;
      }

      if (debug > -1)
        cout << i << " i " << t_i[i][1] << " " << t_i[i][2] << " " << t_i[i][3]
             << endl;
      count_im = prof_sum_eff[i] * t_i[i][1] + rtrans0[2][1] * reg_weight;
      count_ii = prof_sum_eff[i] * t_i[i][2] + rtrans0[2][2] * reg_weight;
      count_id = prof_sum_eff[i] * t_i[i][3] + rtrans0[2][3] * reg_weight;
      total_count = count_im + count_ii + count_id;
      t_i[i][1] = count_im / total_count;
      t_i[i][2] = count_ii / total_count;
      t_i[i][3] = count_id / total_count;
      if (debug > -1)
        cout << i << " i " << t_i[i][1] << " " << t_i[i][2] << " " << t_i[i][3]
             << endl;
      if (debug > -1) cout << endl;
    }
  }

  // predicted secondary structure dependent prior for transition state
  // probabilities
  else {
    for (i = 1; i < prof_len; i++) {
      count_mm = prof_sum_eff[i] * t_m[i][1] +
                 rtrans[x->prof_alphabet1[i]][1][1] * reg_weight;
      count_mi = prof_sum_eff[i] * t_m[i][2] +
                 rtrans[x->prof_alphabet1[i]][1][2] * reg_weight;
      count_md = prof_sum_eff[i] * t_m[i][3] +
                 rtrans[x->prof_alphabet1[i]][1][3] * reg_weight;
      total_count = count_mm + count_mi + count_md;
      t_m[i][1] = (1 - t_m[i][4]) * count_mm / total_count;
      t_m[i][2] = (1 - t_m[i][4]) * count_mi / total_count;
      t_m[i][3] = (1 - t_m[i][4]) * count_md / total_count;
      if (debug > -1)
        cout << i << " m " << t_m[i][1] << " " << t_m[i][2] << " " << t_m[i][3]
             << " " << t_m[i][4] << endl;

      if (i > 1) {
        if (debug > -1)
          cout << i << " d " << t_d[i][1] << " " << t_d[i][2] << " "
               << t_d[i][3] << endl;
        count_dm = prof_sum_eff[i] * t_d[i][1] +
                   rtrans[x->prof_alphabet1[i]][3][1] * reg_weight;
        count_di = prof_sum_eff[i] * t_d[i][2] +
                   rtrans[x->prof_alphabet1[i]][3][2] * reg_weight;
        count_dd = prof_sum_eff[i] * t_d[i][3] +
                   rtrans[x->prof_alphabet1[i]][3][3] * reg_weight;
        total_count = count_dm + count_di + count_dd;
        t_d[i][1] = count_dm / total_count;
        t_d[i][2] = count_di / total_count;
        t_d[i][3] = count_dd / total_count;
        if (debug > -1)
          cout << i << " d " << t_d[i][1] << " " << t_d[i][2] << " "
               << t_d[i][3] << endl;
      }

      if (debug > -1)
        cout << i << " i " << t_i[i][1] << " " << t_i[i][2] << " " << t_i[i][3]
             << endl;
      count_im = prof_sum_eff[i] * t_i[i][1] +
                 rtrans[x->prof_alphabet1[i]][2][1] * reg_weight;
      count_ii = prof_sum_eff[i] * t_i[i][2] +
                 rtrans[x->prof_alphabet1[i]][2][2] * reg_weight;
      count_id = prof_sum_eff[i] * t_i[i][3] +
                 rtrans[x->prof_alphabet1[i]][2][3] * reg_weight;
      total_count = count_im + count_ii + count_id;
      t_i[i][1] = count_im / total_count;
      t_i[i][2] = count_ii / total_count;
      t_i[i][3] = count_id / total_count;
      if (debug > -1)
        cout << i << " i " << t_i[i][1] << " " << t_i[i][2] << " " << t_i[i][3]
             << endl;
      if (debug > -1) cout << endl;
    }
  }

  // without the wing
  t_m[lenx - 1][1] += t_m[lenx - 1][3];
  t_m[lenx - 1][3] = 0;

  t_i[lenx - 1][1] += t_i[lenx - 1][3];
  t_i[lenx - 1][3] = 0;

  t_d[lenx - 1][1] += t_d[lenx - 1][3];
  t_d[lenx - 1][3] = 0;
}
void profilehmm::print_transitions() {
  int i, j, k;

  cout << "Transition probabilities";

  cout << "t_sb: " << t_sb << endl;
  cout << "t_sn: " << t_sn << endl;
  cout << "t_nb: " << t_nb << endl;
  cout << "t_nn: " << t_nn << endl;

  cout << "t_bm: " << endl;
  for (i = 1; i <= lenx; i++) {
    cout << t_bm[i] << endl;
  }
  cout << "t_m, t_i, t_d: " << endl;
  for (i = 1; i <= lenx; i++) {
    cout << "# " << i << endl;
    for (j = 1; j <= 4; j++) {
      cout << t_m[i][j] << " ";
    }
    cout << endl;
    for (j = 1; j <= 3; j++) {
      cout << t_i[i][j] << " ";
    }
    cout << endl;
    for (j = 1; j <= 3; j++) {
      cout << t_d[i][j] << " ";
    }
    cout << endl;
  }

  cout << "t_ec: " << t_ec << endl;
  cout << "t_et: " << t_et << endl;
  cout << "t_cc: " << t_cc << endl;
  cout << "t_ct: " << t_ct << endl;
}

void profilehmm::log_convert() {
  int i, j, k;

  assert(t_sn > 0);
  assert(t_sb > 0);
  assert(t_nn > 0);
  assert(t_nb > 0);
  t_sn = log(t_sn);
  t_sb = log(t_sb);
  t_nn = log(t_nn);
  t_nb = log(t_nb);

  for (i = 1; i <= lenx; i++) {
    assert(t_bm[i] > 0);
    t_bm[i] = log(t_bm[i]);
  }

  for (i = 1; i <= lenx; i++) {
    for (j = 1; j <= 3; j++) {
      if (t_m[i][j] < 0) {
        cout << "t_m is less than 0: " << i << " " << j << " " << t_m[i][j]
             << endl;
        exit(0);
      }
      if (t_d[i][j] < 0) {
        cout << "t_d is less than 0: " << i << " " << j << " " << t_d[i][j]
             << endl;
        exit(0);
      }
      if (t_i[i][j] < 0) {
        cout << "t_i is less than 0: " << i << " " << j << " " << t_i[i][j]
             << endl;
        exit(0);
      }
      if (t_m[i][j] == 0)
        t_m[i][j] = LOG_ZERO;
      else
        t_m[i][j] = log(t_m[i][j]);
      if (t_d[i][j] == 0)
        t_d[i][j] = LOG_ZERO;
      else
        t_d[i][j] = log(t_d[i][j]);
      if (t_i[i][j] == 0)
        t_i[i][j] = LOG_ZERO;
      else
        t_i[i][j] = log(t_i[i][j]);
    }
    assert(t_m[i][4] > 0);
    t_m[i][4] = log(t_m[i][4]);
  }

  assert(t_ec > 0);
  assert(t_et > 0);
  assert(t_ct > 0);
  assert(t_cc > 0);
  t_ec = log(t_ec);
  t_et = log(t_et);
  t_ct = log(t_ct);
  t_cc = log(t_cc);
}

void profilehmm::get_scores(int bg_type) {
  int i, j, k;

  assert(bg_type >= 0);
  assert(bg_type <= 2);
  assert(y != 0);
  if (y->bg_type_here != bg_type) {
    y->get_score_bg(bg_type);
  }
  score_bg = gvector<float>(y->prof_len);
  if (debug > -1) cout << "score_bgs" << endl;
  for (i = 1; i <= y->prof_len; i++) {
    score_bg[i] = y->score_bg_aa[i] + y->score_bg_ss[i] * ss_weight;
    if (debug > -1)
      cout << i << " " << y->score_bg_aa[i] << " " << y->score_bg_ss[i] << " "
           << score_bg[i] << endl;
  }

  score_matrix = gmatrix<float>(lenx, leny);
  if (debug > -1) cout << "score matrix" << endl;
  for (i = 1; i <= lenx; i++) {
    for (j = 1; j <= leny; j++) {
      score_matrix[i][j] = 0;
      for (k = 1; k <= 20; k++) {
        score_matrix[i][j] += y->prof_effn[j][k] * log(x->prof_freq[i][k]);
        if (debug > -1)
          cout << i << " " << j << " " << k << " " << y->prof_effn[j][k] << " "
               << log(x->prof_freq[i][k]) << endl;
      }
      for (k = 1; k <= 3; k++) {
        score_matrix[i][j] += ss_weight * y->prof_sum_eff[j] *
                              y->ss->ssfreq[j][k] * log(x->ss->ssfreq[i][k]);
        if (debug > -1)
          cout << i << " " << j << " " << k << " " << y->prof_sum_eff[j] << " "
               << y->ss->ssfreq[j][k] << " " << log(x->ss->ssfreq[i][k])
               << endl;
      }
      score_matrix[i][j] /= y->prof_sum_eff[j];
      if (debug > 1)
        cout << i << " " << j << " " << score_matrix[i][j] / y->prof_sum_eff[j]
             << endl;
      if (debug > -1)
        cout << i << " " << j << " " << score_matrix[i][j] << endl;
    }
  }
}

void profilehmm::forward() {
  int i, j, k;

  Fn = gvector<float>(leny);
  Fb = gvector<float>(leny);
  Fe = gvector<float>(leny);
  Fc = gvector<float>(leny);

  Fm = gmatrix<float>(lenx, leny);
  Fi = gmatrix<float>(lenx, leny);
  Fd = gmatrix<float>(lenx, leny);

  // Fn
  Fn[1] = t_sn + score_bg[1];
  if (debug > -1) cout << "Fn" << endl;
  for (i = 2; i <= leny; i++) {
    Fn[i] = t_nn + Fn[i - 1] + score_bg[i];
    if (debug > -1) cout << i << " " << Fn[i] << endl;
  }

  // Fb
  Fb[0] = t_sb;
  if (debug > -1) cout << "Fb" << endl;
  for (i = 1; i <= leny; i++) {
    Fb[i] = t_nb + Fn[i];
    if (debug > -1) cout << i << " " << Fb[i] << endl;
  }

  // Fm, Fi, Fd at position 1
  for (i = 1; i <= lenx; i++) {
    Fm[i][0] = LOG_ZERO;
    Fi[i][0] = LOG_ZERO;
    Fd[i][0] = LOG_ZERO;
  }
  for (j = 1; j <= leny; j++) {
    Fm[1][j] = Fb[j - 1] + t_bm[1] + score_matrix[1][j];
    if (debug > -1)
      cout << "Fm "
           << "1 " << j << " " << Fm[1][j] << endl;
    Fi[1][j] = Fm[1][j - 1] + t_m[1][2] + score_bg[j];
    Fd[1][j] = LOG_ZERO;
  }
  // Fm
  for (i = 2; i <= lenx; i++) {
    for (j = 1; j <= leny; j++) {
      Fm[i][j] =
          LOG_ADD(Fm[i - 1][j - 1] + t_m[i - 1][1],
                  Fi[i - 1][j - 1] + t_i[i - 1][1],
                  Fd[i - 1][j - 1] + t_d[i - 1][1], Fb[j - 1] + t_bm[i]) +
          score_matrix[i][j];
      if (debug > -1) cout << "Fm " << i << " " << j << " " << Fm[i][j] << endl;
      Fi[i][j] = LOG_ADD(Fm[i][j - 1] + t_m[i][2], Fi[i][j - 1] + t_i[i][2],
                         Fd[i][j - 1] + t_d[i][2]) +
                 score_bg[j];
      if (debug > -1) cout << "Fi " << i << " " << j << " " << Fi[i][j] << endl;
      Fd[i][j] =
          LOG_ADD(Fm[i - 1][j] + t_m[i - 1][3], Fi[i - 1][j] + t_i[i - 1][3],
                  Fd[i - 1][j] + t_d[i - 1][3]);
      if (debug > -1) cout << "Fd " << i << " " << j << " " << Fd[i][j] << endl;
    }
  }

  // Fe
  if (debug > -1) cout << "Fe" << endl;
  for (j = 1; j <= leny; j++) {
    Fe[j] = LOG_ZERO;
    for (i = 1; i <= lenx; i++) {
      Fe[j] = LOG_ADD(Fe[j], Fm[i][j] + t_m[i][4]);
    }
    if (debug > -1) cout << j << " " << Fe[j] << endl;
  }

  // Fc
  if (debug > -1) cout << "Fc" << endl;
  for (j = 1; j <= leny; j++) {
    if (j == 1) {
      Fc[j] = LOG_ZERO;
      continue;
    }
    Fc[j] = LOG_ADD(Fc[j - 1] + t_cc, Fe[j - 1] + t_ec) + score_bg[j];
    if (debug > -1) cout << j << " " << Fc[j] << endl;
  }

  Ft = LOG_ADD(Fc[leny] + t_ct, Fe[leny] + t_et);

  cout << "Ft: " << Ft << endl;
}

void profilehmm::backward() {
  int i, j, k;

  Bn = gvector<float>(leny);
  Bb = gvector<float>(leny);
  Be = gvector<float>(leny);
  Bc = gvector<float>(leny);

  Bm = gmatrix<float>(lenx, leny + 1);
  Bi = gmatrix<float>(lenx, leny + 1);
  Bd = gmatrix<float>(lenx, leny + 1);

  // Bc and Be
  Bc[leny] = t_ct;
  Be[leny] = t_et;
  if (debug > -1) cout << "Bc" << endl;
  for (j = leny - 1; j >= 1; j--) {
    Bc[j] = Bc[j + 1] + t_cc + score_bg[j + 1];
    if (debug > -1) cout << j << " " << Bc[j] << endl;
  }
  if (debug > -1) cout << "Be" << endl;
  for (j = leny - 1; j >= 1; j--) {
    Be[j] = Bc[j + 1] + t_ec + score_bg[j + 1];
    if (debug > -1) cout << j << " " << Be[j] << endl;
  }

  // Bm, Bi, Bd
  for (j = leny; j >= 1; j--) {
    Bm[lenx][j] = Be[j] + t_m[lenx][4];
    Bi[lenx][j] = LOG_ZERO;
    Bd[lenx][j] = LOG_ZERO;
    cout << "Bm: " << lenx << " " << j << " " << Bm[lenx][j] << endl;
    cout << "Bi: " << lenx << " " << j << " " << Bi[lenx][j] << endl;
    cout << "Bd: " << lenx << " " << j << " " << Bd[lenx][j] << endl;
  }
  for (i = lenx; i >= 1; i--) {
    Bm[i][leny] = Be[leny] + t_m[i][4];
    Bi[i][leny] = LOG_ZERO;
    Bd[i][leny] = LOG_ZERO;
    cout << "Bm: " << i << " " << leny << " " << Bm[i][leny] << endl;
    cout << "Bi: " << i << " " << leny << " " << Bi[i][leny] << endl;
    cout << "Bd: " << i << " " << leny << " " << Bd[i][leny] << endl;
  }

  for (i = lenx - 1; i >= 1; i--) {
    for (j = leny - 1; j >= 1; j--) {
      Bm[i][j] =
          LOG_ADD(Bm[i + 1][j + 1] + t_m[i][1] + score_matrix[i + 1][j + 1],
                  Bi[i][j + 1] + t_m[i][2] + score_bg[j + 1],
                  Bd[i + 1][j] + t_m[i][3], Be[j] + t_m[i][4]);
      Bi[i][j] = LOG_ADD(
          Bm[i + 1][j + 1] + t_i[i][1] + score_matrix[i + 1][j + 1],
          Bi[i][j + 1] + t_i[i][2] + score_bg[j + 1], Bd[i + 1][j] + t_i[i][3]);
      Bd[i][j] = LOG_ADD(
          Bm[i + 1][j + 1] + t_d[i][1] + score_matrix[i + 1][j + 1],
          Bi[i][j + 1] + t_d[i][2] + score_bg[j + 1], Bd[i + 1][j] + t_d[i][3]);
      cout << "Bm: " << i << " " << j << " " << Bm[i][j] << endl;
      cout << "Bi: " << i << " " << j << " " << Bi[i][j] << endl;
      cout << "Bd: " << i << " " << j << " " << Bd[i][j] << endl;
    }
  }

  // Bb, Bn
  cout << "Bb" << endl;
  for (j = 0; j <= leny; j++) {
    Bb[j] = LOG_ZERO;
    if (j == leny) break;
    for (i = 1; i <= lenx; i++) {
      Bb[j] = LOG_ADD(Bm[i][j + 1] + t_bm[i] + score_matrix[i][j + 1], Bb[j]);
    }
    cout << j << " " << Bb[j] << endl;
  }

  Bn[leny] = Bb[leny] + t_nb;
  cout << "Bn" << endl;
  for (j = leny - 1; j >= 1; j--) {
    Bn[j] = LOG_ADD(Bb[j] + t_nb, Bn[j + 1] + t_nn + score_bg[j + 1]);
    cout << j << " " << Bn[j] << endl;
  }

  Bs = LOG_ADD(Bb[0] + t_sb, Bn[1] + t_sn + score_bg[1]);

  cout << "Bs: " << Bs << endl;
}

void profilehmm::get_match_prob() {
  int i, j, k;

  probMat = gmatrix<float>(lenx, leny);
  cout << "probMat:" << endl;
  for (i = 1; i <= lenx; i++) {
    for (j = 1; j <= leny; j++) {
      probMat[i][j] = exp(Fm[i][j] + Bm[i][j] - Ft);
      cout << i << " " << j << " " << probMat[i][j] << endl;
    }
  }
}

/*
void profilehmm::viterbi() {

        int i,j,k;
        int xi, yj, a;

        float mm, mxy, xym, xy, E;
        int maxalnlength = lenx + leny +1;
        int *reverse_path, *reverse_path1, *reverse_path2;

        mm = log(1-2*delta-tau);
        mxy = log(1-epsilon-tau);
        xym = log(delta);
        xy = log(epsilon);
        //E = log(tau);
        E = 0;

        gapopen_end = 1;

        // initialize the matrices
        Vm = gmatrix<float>(lenx, leny);
        Vx = gmatrix<float>(lenx, leny);
        Vy = gmatrix<float>(lenx, leny);
        Tm = imatrix(lenx, leny);
        Tx = imatrix(lenx, leny);
        Ty = imatrix(lenx, leny);
        path = ivector(maxalnlength);
        path1 = ivector(maxalnlength);
        path2 = ivector(maxalnlength);
        reverse_path = ivector(maxalnlength);
        reverse_path1 = ivector(maxalnlength);
        reverse_path2 = ivector(maxalnlength);

        for(i=0;i<=lenx;i++) Vm[i][0] = Vx[i][0] = Vy[i][0] = -1000000;
        for(i=0;i<=leny;i++) Vm[0][i] = Vx[0][i] = Vy[0][i] = -1000000;
        //Vm[0][0] = 0; //Vx[0][0] = Vy[0][0] = 0;
        Vm[0][0] = log(theta) - xy; //Vx[0][0] = Vy[0][0] = 0;
        for(i=1;i<=lenx;i++) {
                if(!gapopen_end) Vx[i][0] = max2(xy+Vm[i-1][0], xy+Vx[i-1][0],
Tx[i][0]);
                // for end gap, waive the gap opening penalty: change
                else Vx[i][0] = max2(mxy+Vm[i-1][0], xy+Vx[i-1][0], Tx[i][0]);
                if(debug>1) cout << "Vx i 0: " << i << " " << Vx[i][0] << endl;
        }
        for(j=1;j<=leny;j++) {
                if(!gapopen_end) Vy[0][j] = max2(xy+Vm[0][j-1], xy+Vy[0][j-1],
Ty[0][j]); else Vy[0][j] = max2(mxy+Vm[0][j-1], xy+Vy[0][j-1], Ty[0][j]);
                if(debug>1) cout << "Vy 0 j: " << j << " " << Vy[0][j] << endl;
        }

        Vm[0][0] = log(1-2*theta);
        for(i=1;i<=lenx;i++) {
           for(j=1;j<=leny;j++) {
                // testing
                // Vm[i][j] = log_odds_score(i,j) + max3(mm+Vm[i-1][j-1],
mxy+Vx[i-1][j-1], mxy+Vy[i-1][j-1], Tm[i][j]);
                // Vx[i][j] = log_odds_x_d(i) + max2(mxy+Vm[i-1][j],
xy+Vx[i-1][j], Tx[i][j]);
                // Vy[i][j] = log_odds_y_d(j) + max2(mxy+Vm[i][j-1],
xy+Vy[i][j-1], Ty[i][j]);

                Vm[i][j] = log_odds_score(i,j) + max3(mm+Vm[i-1][j-1],
mxy+Vx[i-1][j-1], mxy+Vy[i-1][j-1], Tm[i][j]);
                // for end gap, waive the gap opening penalty: change
                if(j==leny) {
                     if(!gapopen_end) Vx[i][j] = max2(xy+Vm[i-1][j],
xy+Vx[i-1][j], Tx[i][j]); else Vx[i][j] = max2(mxy+Vm[i-1][j], xy+Vx[i-1][j],
Tx[i][j]);
                }
                else Vx[i][j] = max2(mxy+Vm[i-1][j], xy+Vx[i-1][j], Tx[i][j]);
                // for end gap, waive the gap opening penalty: change
                if(i==lenx) {
                     if(!gapopen_end) Vy[i][j] = max2(xy+Vm[i][j-1],
xy+Vy[i][j-1], Ty[i][j]); else Vy[i][j] = max2(mxy+Vm[i][j-1], xy+Vy[i][j-1],
Ty[i][j]);
                }
                else Vy[i][j] = max2(mxy+Vm[i][j-1], xy+Vy[i][j-1], Ty[i][j]);

                if(debug>1) fprintf(stdout, "%d %d %f %f %f %d %d %d %f\n", i,
j, Vm[i][j], Vx[i][j], Vy[i][j], Tm[i][j], Tx[i][j], Ty[i][j],
log_odds_score(i,j));
           }
        }

        VE = E + max3(Vm[lenx][leny]+LOG(1-2*theta), Vx[lenx][leny]+LOG(theta),
Vy[lenx][leny]+LOG(theta), TE);

        // trace back
        xi = lenx;
        yj = leny;
        int mark_mxy;
        if(TE==1) {
                reverse_path[0] = 0;
                reverse_path1[0] = 0;
                reverse_path2[0] = 0;
                mark_mxy = 1;
                //xi--; yj--;
        }
        else if(TE==2) {
                reverse_path[0] = 1;
                reverse_path1[0] = 0;
                reverse_path2[0] = 1;
                mark_mxy = 2;
                //xi--;
        }
        else {
                reverse_path[0] = -1;
                reverse_path1[0] = 1;
                reverse_path2[0] = 0;
                mark_mxy = 3;
                //yj--;
        }
        i = 1;
        //fprintf(stdout, "=================\n");
        while( (xi!=0) || (yj!=0) ) {
                if(debug>1) fprintf(stdout, "%d %d    xi: %d yj: %d Tm %d Tx %d
Ty %d\n", i, mark_mxy, xi, yj, Tm[xi][yj], Tx[xi][yj], Ty[xi][yj]);
                if(mark_mxy==1) {
                        a = Tm[xi][yj];
                        if(a==1) {
                                //reverse_path[i+1] = 0;
                                //xi--; yj--;
                                mark_mxy = 1;
                        }
                        else if(a==2) {
                                //reverse_path[i+1] = 1;
                                //xi--;
                                mark_mxy = 2;
                        }
                        else if(a==3) {
                                //reverse_path[i+1] = -1;
                                //yj--;
                                mark_mxy = 3;
                        }
                        reverse_path[i] = 0;
                        reverse_path1[i] = 0; // 0 means amino acid
                        reverse_path2[i] = 0;
                        xi--; yj--;
                        i++;
                        continue;
                }
                if(mark_mxy==2) {
                        a = Tx[xi][yj];
                        if(a==1) {
                                //reverse_path[i+1] = 0;
                                // xi--; yj--;
                                mark_mxy = 1;
                        }
                        else if(a==2) {
                                //reverse_path[i+1] = 1;
                                // xi--;
                                mark_mxy = 2;
                        }
                        reverse_path[i] = 1;
                        reverse_path1[i] = 0;
                        reverse_path2[i] = 1; // 1 means gap
                        xi--;
                        i++;
                        continue;
                }
                if(mark_mxy==3) {
                        a = Ty[xi][yj];
                        if(a==1) {
                                //reverse_path[i+1] = 0;
                                // xi--; yj--;
                                mark_mxy = 1;
                        }
                        else if(a==2) {
                                //reverse_path[i+1] = -1;
                                // yj--;
                                mark_mxy = 3;
                        }
                        reverse_path[i] = -1;
                        reverse_path1[i] = 1; // 1 means gap
                        reverse_path2[i] = 0;
                        yj--;
                        i++;
                        continue;
                }
        }
        //fprintf(stdout, "=================\n");

        v_alilen = i-1;
        if(debug>1) fprintf(stdout, "v_alilen: %d\n", v_alilen);
        for(i=1;i<=v_alilen;i++) {
                path[i] = reverse_path[v_alilen-i+1];
                path1[i] = reverse_path1[v_alilen-i+1];
                path2[i] = reverse_path2[v_alilen-i+1];
                if(debug>1) fprintf(stdout, "%d %d %d\n", i, reverse_path1[i],
reverse_path2[i]);
        }
        if(debug>1) fprintf(stdout, "\n");
        delete [] reverse_path;
        delete [] reverse_path1;
        delete [] reverse_path2;

}

// this is correponds to sum-of-pairs of the log-odds RATIO scores of weighted
amino acid pairs float profilehmm::log_odds_score(int xi, int yj) {

        int i,j;

        float score = 0;

        for(i=1;i<=20;i++) {
           if(!x->pseudoCnt[xi][i]) continue;
           for(j=1;j<=20;j++) {
                ////score += (x->pseudoCnt[xi][i] * y->pseudoCnt[yj][j] *
log(q_blosum62[i][j]/robinson_freq[i]/robinson_freq[j]) );
                if(!y->pseudoCnt[yj][j]) continue;
                score += (x->pseudoCnt[xi][i] * y->pseudoCnt[yj][j] *
log_q_blosum62_ratio[i][j]);
                // testing
                // score += (x->pseudoCnt[xi][i] * y->pseudoCnt[yj][j] *
log(q_blosum62[i][j]) );
           }
        }

        // if(xi==lenx) if(yj==leny) { for(i=1;i<=20;i++) {
        //		fprintf(stdout, "*  %d %f  %f\n", i,
x->pseudoCnt[xi][i], y->pseudoCnt[yj][i]);
        //	}
        //}

        return score;
}

float profilehmm::max3(float a,float b, float c, int &t) {

        float max;

        if(a>=b) {
                if(a>=c) {
                        max = a;
                        t = 1;
                }
                else {
                        max = c;
                        t = 3;
                }
        }
        else {
                if(b>c) {
                        max = b;
                        t = 2;
                }
                else {
                        max = c;
                        t = 3;
                }
        }
        return max;
}


float profilehmm::max2(float a, float b, int &t) {

        if(a>=b) {
                t= 1;
                return a;
        }
        else {
                t = 2;
                return b;
        }
}

void profilehmm::printHmmAlign(int blocksize) {

        int i,j,k;

        int nblocks;

        int max_name_len;

        if(x->mnamelen > y->mnamelen) max_name_len = x->mnamelen;
        else max_name_len =  y->mnamelen;

        fprintf(stdout, "%d %d\n", x->mnamelen, y->mnamelen);
        fprintf(stdout, "\n");
        fprintf(stdout, "The final alignment: \n\n");

        // number of blocks
        nblocks = (v_alilen-1)/blocksize + 1;

        int a = -1;
        int b = -1;
        int a_prev, b_prev, a_after, b_after;
        for(i=1; i<=nblocks; i++) {
                //fprintf(stdout, "=================\n\n");
                a_prev = a;
                for(j=1; j<=x->nal; j++) {
                        cout << setiosflags(ios::left) << setw(max_name_len+3)
<< x->aname[j-1]; for(k=(i-1)*blocksize+1; k<=i*blocksize; k++) { if(k>v_alilen)
break; if(path[k]!=-1) { a++; cout << x->aseq[j-1][a];
                                }
                                else cout << "*";
                        }
                        fprintf(stdout, "\n");
                        a_after = a;
                        a = a_prev;
                }
                a = a_after;
                fprintf(stdout, "\n");

                b_prev = b;
                for(j=1; j<=y->nal; j++) {
                        cout << setiosflags(ios::left) << setw(max_name_len+3)
<< y->aname[j-1]; for(k=(i-1)*blocksize+1; k<=i*blocksize; k++) { if(k>v_alilen)
break; if(path[k]!=1) { b++; cout << y->aseq[j-1][b];
                                }
                                else cout << "*";
                        }
                        fprintf(stdout, "\n");
                        b_after = b;
                        b = b_prev;
                }
                b = b_after;
                fprintf(stdout, "\n\n");
        }
}

// generate an alignment of the two alignments
subalign *profilehmm::productViterbiAlign() {

        int i,j,k;

        //subalign *newAlign = new subalign();
        subalign *a = new subalign(); // newAlign;

        a->nal = x->nal + y->nal;
        a->alilen = v_alilen;
        if(x->mnamelen > y->mnamelen) a->mnamelen = x->mnamelen;
        else a->mnamelen = y->mnamelen;
        if(debug>1) cout << "a->mnamelen: " << a->mnamelen <<endl;
        if(debug>1) cout << "a->alilen: " << a->alilen <<endl;
        if(debug>1) cout << "a->nal: " << a->nal <<endl;
        a->aname = cmatrix(a->nal, a->mnamelen+1);
        a->aseq = cmatrix(a->nal, a->alilen+1);
        a->alignment = imatrix(a->nal, a->alilen);

        // copy the names
        for(i=0;i<x->nal;i++) {
                strcpy(a->aname[i], x->aname[i]);
        }
        for(i=x->nal;i<a->nal;i++) {
                strcpy(a->aname[i], y->aname[i - x->nal]);
        }

        if(debug>1) for(i=1;i<=v_alilen;i++) { cout << i << "\t" << path[i] <<
endl; }

        // derive the sequences from alignment path
        int xi=0, yi=0;
        for(i=0;i<a->alilen;i++) {
                for(j=0;j<x->nal;j++) {
                        if(path[i+1]>=0) {
                                a->aseq[j][i] = x->aseq[j][xi];
                        }
                        else {
                                a->aseq[j][i] = '-';
                        }
                }
                if(path[i+1]>=0) xi++;
                for(j=x->nal;j<a->nal;j++) {
                        if(path[i+1]<=0) {
                                a->aseq[j][i] = y->aseq[j - x->nal][yi];
                        }
                        else {
                                a->aseq[j][i] = '-';
                        }
                }
                if(path[i+1]<=0) yi++;
                if(debug>1) cout << "i: " << i<< endl;
        }
        for(i=0;i<a->nal;i++) {
                a->aseq[i][a->alilen]='\0';
        }

        // convert the letters to numbers
        a->convertAseq2Alignment();

        return a;
}

void profilehmm::forward() {

        int i,j,k;

        ScoreType mm, mxy, xym, xy, E;

        mm = LOG( 1-2*delta-tau );
        mxy = LOG( 1-epsilon-tau );
        xym = LOG( delta );
        xy = LOG( epsilon );
        //E = LOG( tau );
        E = LOG_ONE;

        //cout<<"Forward algorithm:"<<endl;

        Fm = new ScoreType * [lenx+1];
        Fx = new ScoreType * [lenx+1];
        Fy = new ScoreType * [lenx+1];

        for(i=0;i<=lenx;i++) {
                Fm[i] = new ScoreType [leny+1];
                Fx[i] = new ScoreType [leny+1];
                Fy[i] = new ScoreType [leny+1];
        }

        if(debug>1) fprintf(stdout, "===================\n");
        Fm[0][0] = LOG_ONE; Fx[0][0] = LOG_ZERO; Fy[0][0] = LOG_ZERO;
        if(debug>1) fprintf(stdout, "===================\n");
        for(i=0;i<=lenx;i++) {
                for(j=0;j<=leny;j++) {
                        if( (i==0)&&(j==0) ) {; }
                        //else if( (i==0)&&(j==1) ) {
                        else if( (i==0)&&(j!=0) ) {
                                Fm[i][j] = LOG_ZERO;
                                Fx[i][j] = LOG_ZERO;
                                //Fy[i][j] = xym * Fm[i][j-1] + xy * Fy[i][j-1];
                                // waive gap open penalty for ends
                                ////Fy[i][j] = xy * Fm[i][j-1] + xy *
Fy[i][j-1];
                                ////Fy[i][j] = Fy[i][j] * log_odds_y(j);
                                if(j==1) {
                                        Fy[i][j] = LOG(theta) + log_odds_y(j);
                                }
                                else Fy[i][j] = LOG_ADD(xy+Fm[i][j-1],
xy+Fy[i][j-1])+log_odds_y(j);
                        }
                        //else if( (i==1)&&(j==0) ) {
                        else if( (i!=0)&&(j==0) ) {
                                Fm[i][j] = LOG_ZERO;
                                Fy[i][j] = LOG_ZERO;
                                //Fx[i][j] = xym * Fm[i-1][j] + xy * Fx[i-1][j];
                                // waive gap open penalty for ends
                                //Fx[i][j] = xy * Fm[i-1][j] + xy * Fx[i-1][j];
                                //Fx[i][j] = Fx[i][j] * log_odds_x(i);
                                if(i==1) {
                                        Fx[i][j] = LOG(theta) + log_odds_x(i);
                                }
                                else Fx[i][j] = LOG_ADD(xy+Fm[i-1][j],
xy+Fx[i-1][j])+log_odds_x(i);
                        }
                        else {
                                ////Fm[i][j] = mm * Fm[i-1][j-1] + mxy *
Fx[i-1][j-1] + mxy * Fy[i-1][j-1];
                                ////Fm[i][j] = Fm[i][j] * log_odds(i,j);
                                if( (i==1) && (j==1) ) {
                                        Fm[i][j] = LOG(1-2*theta) +
log_odds(i,j); } else Fm[i][j] = LOG_ADD(mm+Fm[i-1][j-1],mxy+Fx[i-1][j-1],
mxy+Fy[i-1][j-1])+log_odds(i,j);
                                // waive gap open penalty for ends
                                ////if(j==leny)Fx[i][j] = xy * Fm[i-1][j] + xy *
Fx[i-1][j];
                                ////else Fx[i][j] = xym * Fm[i-1][j] + xy *
Fx[i-1][j];
                                //Fx[i][j] = xym * Fm[i-1][j] + xy * Fx[i-1][j];
                                ////Fx[i][j] = Fx[i][j] * log_odds_x(i);

                                //** waive gap open penalty for ends
                                //**if(j==leny)Fx[i][j] = LOG_ADD(xy+Fm[i-1][j],
xy+Fx[i-1][j]);
                                //**else Fx[i][j] = LOG_ADD(xym+Fm[i-1][j],
xy+Fx[i-1][j]);
                                // not waiving gap open penalty for ends
                                Fx[i][j] = LOG_ADD(xym+Fm[i-1][j],
xy+Fx[i-1][j]); Fx[i][j] = Fx[i][j] + log_odds_x(i);

                                // waive gap open penalty for ends
                                ////if(i==lenx)Fy[i][j] = xy * Fm[i][j-1] + xy *
Fy[i][j-1];
                                ////else Fy[i][j] = xym * Fm[i][j-1] + xy *
Fy[i][j-1];
                                //Fy[i][j] = xym * Fm[i][j-1] + xy * Fy[i][j-1];
                                ////Fy[i][j] = Fy[i][j] * log_odds_y(j);

                                //** waive gap open penalty for ends
                                //** if(i==lenx)Fy[i][j] =
LOG_ADD(xy+Fm[i][j-1],xy+Fy[i][j-1]);
                                //** else Fy[i][j] =
LOG_ADD(xym+Fm[i][j-1],xy+Fy[i][j-1]);
                                // not waiving gap open penalty for ends
                                Fy[i][j] =
LOG_ADD(xym+Fm[i][j-1],xy+Fy[i][j-1]); Fy[i][j] = Fy[i][j] + log_odds_y(j);
                        }
                        ////if(debug>1) fprintf(stdout, "%d %d: %d %f   %d %f %d
%f\n",
i,j,Fm[i][j].INT,Fm[i][j].DEC,Fx[i][j].INT,Fx[i][j].DEC,Fy[i][j].INT,Fy[i][j].DEC);
                        if(debug>1) fprintf(stdout, "%d %d: %f   %f   %f\n",
i,j,Fm[i][j],Fx[i][j],Fy[i][j]); fflush(stdout);
                }
        }

        ////FE = E * (Fm[lenx][leny] + Fx[lenx][leny] + Fy[lenx][leny]);
        FE = E +
LOG_ADD(Fm[lenx][leny]+LOG(1-2*theta),Fx[lenx][leny]+LOG(theta),Fy[lenx][leny]+LOG(theta));
        if(debug>1) fprintf(stdout, "FE: %f\n", FE);

}

// this is correponds to sum-of-pairs of the log-odds scores of weighted amino
acid pairs ScoreType profilehmm::log_odds(int xi, int yj) {

        int i,j;

        ////double score = 0;
        float score = 0;

        ////if(xi>lenx) return ScoreType(0);
        ////if(yj>leny) return ScoreType(0);
        if(xi>lenx) return LOG_ZERO;
        if(yj>leny) return LOG_ZERO;

        for(i=1;i<=20;i++) {
           if(!x->pseudoCnt[xi][i]) continue;
           for(j=1;j<=20;j++) {
                ////score += (x->pseudoCnt[xi][i] * y->pseudoCnt[yj][j] *
log(q_blosum62[i][j]) ); if(!y->pseudoCnt[yj][j]) continue; score +=
(x->pseudoCnt[xi][i] * y->pseudoCnt[yj][j] * log_q_blosum62[i][j] );
           }
        }

        // if(xi==lenx) if(yj==leny) { for(i=1;i<=20;i++) {
        //		fprintf(stdout, "*  %d %f  %f\n", i,
x->pseudoCnt[xi][i], y->pseudoCnt[yj][i]);
        //	}
        //}

        ////return ScoreType(score,1);
        return score;
}

// this is correponds to log-odds scores of weighted amino acid frequencies
ScoreType profilehmm::log_odds_x(int xi) {

        int i,j;

        ////double score = 0;
        float score = 0;

        ////if(xi>lenx) return ScoreType(0);
        if(xi>lenx) return LOG_ZERO;

        for(i=1;i<=20;i++) {
                if(!x->pseudoCnt[xi][i]) continue;
                ////score += (x->pseudoCnt[xi][i] * log(robinson_freq[i]) );
                score += (x->pseudoCnt[xi][i] * log_robinson_freq[i] );
        }

        // if(xi==lenx) if(yj==leny) { for(i=1;i<=20;i++) {
        //		fprintf(stdout, "*  %d %f  %f\n", i,
x->pseudoCnt[xi][i], y->pseudoCnt[yj][i]);
        //	}
        //}

        ////return ScoreType(score,1);
        return score;
}

// this is correponds to log-odds scores of weighted amino acid frequencies
ScoreType profilehmm::log_odds_y(int yj) {

        int i,j;

        ////double score = 0;
        float score = 0;

        ////if(yj>leny) return ScoreType(0);
        if(yj>leny) return LOG_ZERO;

        for(i=1;i<=20;i++) {
                if(!y->pseudoCnt[yj][i]) continue;
                ////score += (y->pseudoCnt[yj][i] * log(robinson_freq[i]) );
                score += (y->pseudoCnt[yj][i] * log_robinson_freq[i] );
        }

        // if(xi==lenx) if(yj==leny) { for(i=1;i<=20;i++) {
        //		fprintf(stdout, "*  %d %f  %f\n", i,
x->pseudoCnt[xi][i], y->pseudoCnt[yj][i]);
        //	}
        //}

        ////return ScoreType(score,1);
        return score;
}

// this is correponds to log-odds scores of weighted amino acid frequencies
float profilehmm::log_odds_x_d(int xi) {

        int i,j;

        ////double score = 0;
        float score = 0;

        //if(xi>lenx) return ScoreType(0);

        for(i=1;i<=20;i++) {
                ////score += (x->pseudoCnt[xi][i] * log(robinson_freq[i]) );
                score += (x->pseudoCnt[xi][i] * log_robinson_freq[i]);
        }

        // if(xi==lenx) if(yj==leny) { for(i=1;i<=20;i++) {
        //		fprintf(stdout, "*  %d %f  %f\n", i,
x->pseudoCnt[xi][i], y->pseudoCnt[yj][i]);
        //	}
        //}

        //return ScoreType(score,1);
        return score;
}

// this is correponds to log-odds scores of weighted amino acid frequencies
float profilehmm::log_odds_y_d(int yj) {

        int i,j;

        float score = 0;

        //if(yj>leny) return ScoreType(0);

        for(i=1;i<=20;i++) {
                ////score += (y->pseudoCnt[yj][i] * log(robinson_freq[i]) );
                score += (y->pseudoCnt[yj][i] * log_robinson_freq[i] );
        }

        // if(xi==lenx) if(yj==leny) { for(i=1;i<=20;i++) {
        //		fprintf(stdout, "*  %d %f  %f\n", i,
x->pseudoCnt[xi][i], y->pseudoCnt[yj][i]);
        //	}
        //}

        //return ScoreType(score,1);
        return score;
}

void profilehmm::backward() {

        int i,j,k;

        ScoreType mm, mxy, xym, xy, E;

        mm = LOG( 1-2*delta-tau );
        mxy = LOG( 1-epsilon-tau );
        xym = LOG( delta );
        xy = LOG( epsilon );
        //E = LOG( tau );
        E = LOG_ONE;

        Bm = new ScoreType * [lenx+2];
        Bx = new ScoreType * [lenx+2];
        By = new ScoreType * [lenx+2];

        for(i=0;i<=lenx+1;i++) {
                Bm[i] = new ScoreType [leny+2];
                Bx[i] = new ScoreType [leny+2];
                By[i] = new ScoreType [leny+2];
        }

        if(debug>1) fprintf(stdout, "===================\n");
        //Bm[lenx][leny] = Bx[lenx][leny] = By[lenx][leny] = E;
        Bm[lenx][leny] = LOG(1-2*theta);
        Bx[lenx][leny] = By[lenx][leny] = LOG(theta);
        for(i=0;i<=lenx+1;i++) Bm[i][leny+1] = Bx[i][leny+1] = By[i][leny+1] =
LOG_ZERO; for(i=0;i<=leny+1;i++) Bm[lenx+1][i] = Bx[lenx+1][i] = By[lenx+1][i] =
LOG_ZERO; if(debug>1) fprintf(stdout, "===================\n");
        for(i=lenx;i>=0;i--) {
                for(j=leny;j>=0;j--) {
                        if( (i==lenx) && (j==leny) ) continue;
                        //** waive gap open penalty at ends
                        ////Bm[i][j] = mm * log_odds(i+1,j+1) * Bm[i+1][j+1] +
xym * (log_odds_x(i+1) * Bx[i+1][j] + log_odds_y(j+1) * By[i][j+1]); Bm[i][j] =
LOG_ADD(mm + log_odds(i+1,j+1) + Bm[i+1][j+1], xym + LOG_ADD(log_odds_x(i+1) +
Bx[i+1][j] , log_odds_y(j+1) + By[i][j+1]));

                        if( (j==0)&&(i==0) ){
                                ////Bx[i][j] = xy * log_odds(i+1,j+1) *
Bm[i+1][j+1] + xy * log_odds_x(i+1) * Bx[i+1][j]; Bx[i][j] = LOG_ADD( xy +
log_odds(i+1,j+1) + Bm[i+1][j+1] , xy + log_odds_x(i+1) + Bx[i+1][j]);
                        }
                        else
                        ////Bx[i][j] = mxy * log_odds(i+1,j+1) * Bm[i+1][j+1] +
xy * log_odds_x(i+1) * Bx[i+1][j]; Bx[i][j] = LOG_ADD( mxy + log_odds(i+1,j+1) +
Bm[i+1][j+1] , xy + log_odds_x(i+1) + Bx[i+1][j]); if( (j==0) &&(i==0) ){
                                ////By[i][j] = xy * log_odds(i+1,j+1) *
Bm[i+1][j+1] + xy * log_odds_y(j+1) * By[i][j+1]; By[i][j] = LOG_ADD( xy +
log_odds(i+1,j+1) + Bm[i+1][j+1] , xy + log_odds_y(j+1) + By[i][j+1]);
                        }
                        else
                        ////By[i][j] = mxy * log_odds(i+1,j+1) * Bm[i+1][j+1] +
xy * log_odds_y(j+1) * By[i][j+1]; By[i][j] = LOG_ADD(mxy + log_odds(i+1,j+1) +
Bm[i+1][j+1] , xy + log_odds_y(j+1) + By[i][j+1]);
                        ////if(debug>1) fprintf(stdout, "%d %d: %d %f   %d %f %d
%f\n",
i,j,Bm[i][j].INT,Bm[i][j].DEC,Bx[i][j].INT,Bx[i][j].DEC,By[i][j].INT,By[i][j].DEC);
                        if(debug>1) fprintf(stdout, "%d %d: %f   %f   %f\n",
i,j,Bm[i][j],Bx[i][j],By[i][j]);
                }
        }

        ////BE = mm * Bm[1][1] * log_odds(1,1) + xy * Bx[1][0] * log_odds_x(1) +
xy * By[0][1] * log_odds_y(1);
        ////fprintf(stdout, "BE: %d %f Bm 0 0: %d %f\n", BE.INT, BE.DEC,
Bm[0][0].INT, Bm[0][0].DEC);
        //BE = LOG_ADD(mm + Bm[1][1] + log_odds(1,1) , xy + Bx[1][0] +
log_odds_x(1) , xy + By[0][1] + log_odds_y(1) ); if(debug>1) BE =
LOG_ADD(LOG(1-2*theta)+Bm[1][1] + log_odds(1,1) , LOG(theta) + Bx[1][0] +
log_odds_x(1) , LOG(theta) + By[0][1] + log_odds_y(1) );
        //fprintf(stdout, "BE: %f Bm 0 0: %f\n", BE, Bm[0][0]);

        // adjust BE:
        BE = (FE + BE)/2;

        //fprintf(stdout, "================\n");

        // not the full probabities
        ////ScoreType full_prob;
        ////for(i=1;i<=lenx;i++) { for(j=1;j<=leny;j++) { full_prob = Fm[i][j] *
Bm[i][j] + Fx[i][j] * Bx[i][j] + Fy[i][j] * By[i][j]; //fprintf(stdout, "%d %d:
%d %f\n", i, j, full_prob.INT, full_prob.DEC); } }

        //fprintf(stdout, "================\n");

        //ScoreType full_prob;
        //cout << "Obtain posterior probabilities " << endl;
        //ScoreType **aligned_pair_prob = gmatrix<ScoreType>(lenx, leny);
        probMat = gmatrix<float>(lenx, leny);
        for(i=1;i<=lenx;i++) {
                for(j=1;j<=leny;j++) {
                        ////aligned_pair_prob[i][j] = Fm[i][j] * Bm[i][j] / BE;
                        //aligned_pair_prob[i][j] = Fm[i][j] + Bm[i][j] - BE;
                        //cout << aligned_pair_prob[i][j] << endl;
                        ////probMat[i][j] = aligned_pair_prob[i][j].real();
                        ////probMat[i][j] = EXP(aligned_pair_prob[i][j]);
                        ////why EXP gives negative values?
                        probMat[i][j] = exp(Fm[i][j] + Bm[i][j] - BE);
                        ////if(debug>2) fprintf(stdout, "%d %d: %d %f %f\n", i,
j, aligned_pair_prob[i][j].INT, aligned_pair_prob[i][j].DEC,
aligned_pair_prob[i][j].real()); if(debug>2) fprintf(stdout, "%d %d: %f\n", i,
j, probMat[i][j]);
                }
        }

        //free_gmatrix<ScoreType>(aligned_pair_prob, lenx, leny);

        //fprintf(stdout, "================\n");

}



void profilehmm::getPosterior() {

        ; // nothing here now

}

void profilehmm::getTransitions() {

        ep = log(epsilon);
        ep1 = log(1-epsilon);
        d = log(delta);
        d2 = log(1-2*delta);
        th = log(theta);
        th2 = log(1-2*theta);

}

*/
