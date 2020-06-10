#include "hmm_psipred.h"
#include "mm.h"
#include "param.h"

static int debug = 1;

void hmm_psipred_parameters::read_parameters(char *filename, int num,
                                             int log_mark_emission) {
  int i, j, k;
  char line[500];
  char tmpstr[500];

  env_number = num;

  // allocations
  aa_pair = new float **[env_number + 1];
  for (i = 0; i <= env_number; i++) {
    aa_pair[i] = gmatrix<float>(20, 20);
  }
  aa_bg1 = gmatrix<float>(env_number, 20);
  aa_loop = gmatrix<float>(env_number, 20);

  ss_pair = gmatrix<float>(env_number, env_number);
  aa_pair1 = new float **[env_number + 1];
  for (i = 0; i <= env_number; i++) {
    aa_pair1[i] = gmatrix<float>(20, 20);
  }
  aa_loop1 = gmatrix<float>(env_number, 20);
  ss_pair1 = gmatrix<float>(env_number, env_number);
  ss_loop = gmatrix<float>(env_number, env_number);
  ss_loop1 = gmatrix<float>(env_number, env_number);

  tran_begin = gvector<float>(3);
  tran_match = gmatrix<float>(env_number, 4);
  tran_x = gmatrix<float>(env_number, 3);
  tran_y = gmatrix<float>(env_number, 3);

  // read the file
  ifstream fp(filename, ios::in);
  if (!fp) {
    cout << "Error: profile hmm (with psipred predictions) paramter file "
         << filename << " not readable" << endl;
    exit(0);
  }

  int begin_count = 0;
  while (fp.good()) {
    fp.getline(line, 500);
    if (strncmp(line, "tran_begin", 9) == 0) {
      begin_count++;
      if (begin_count == 2) break;
    }
  }
  for (i = 1; i <= 3; i++) {
    fp >> tran_begin[i];
    tran_begin[i] = mylog(tran_begin[i]);
  }
  fp >> tmpstr;
  for (i = 0; i <= env_number; i++) {
    for (j = 1; j <= 4; j++) {
      fp >> tran_match[i][j];
      tran_match[i][j] = mylog(tran_match[i][j]);
    }
  }
  fp >> tmpstr;
  for (i = 0; i <= env_number; i++) {
    for (j = 1; j <= 3; j++) {
      fp >> tran_x[i][j];
      tran_x[i][j] = mylog(tran_x[i][j]);
    }
  }
  fp >> tmpstr;
  for (i = 0; i <= env_number; i++) {
    for (j = 1; j <= 3; j++) {
      fp >> tran_y[i][j];
      tran_y[i][j] = mylog(tran_y[i][j]);
    }
  }

  for (i = 0; i <= env_number; i++) {
    fp >> tmpstr;
    fp >> tmpstr;
    for (j = 1; j <= 20; j++) {
      for (k = 1; k <= 20; k++) {
        fp >> aa_pair1[i][j][k];
        // if(log_mark_emission)
        aa_pair[i][j][k] = mylog(aa_pair1[i][j][k]);
      }
    }
    fp >> tmpstr;
    fp >> tmpstr;
    for (j = 1; j <= 20; j++) {
      fp >> aa_loop1[i][j];
      // if(log_mark_emission)
      aa_loop[i][j] = mylog(aa_loop1[i][j]);
    }
    fp >> tmpstr;
    fp >> tmpstr;
    for (j = 1; j <= env_number; j++) {
      fp >> ss_pair1[i][j];
      // if(log_mark_emission)
      ss_pair[i][j] = mylog(ss_pair1[i][j]);
    }
    fp >> tmpstr;
    fp >> tmpstr;
    for (j = 1; j <= env_number; j++) {
      fp >> ss_loop1[i][j];
      ss_loop[i][j] = mylog(ss_loop1[i][j]);
    }
  }

  for (i = 0; i <= env_number; i++) {
    for (j = 1; j <= 20; j++) {
      aa_bg1[i][j] = 0;
      for (k = 1; k <= 20; k++) {
        aa_bg1[i][j] += aa_pair1[i][j][k];
      }
    }
  }

  fp.close();
}

inline float hmm_psipred_parameters::mylog(float a) {
  if (a == 0) return log(0.0000001);
  if (a < 0) {
    cout << "cannot take logarithm of a negative" << endl;
    exit(0);
  }
  return log(a);
}

void hmm_psipred_parameters::print_parameters() {
  int i, j, k;
  cout << fixed;
  cout << setprecision(3);
  cout << setw(7);
  cout << left;

  cout << endl;

  cout << "tran_begin:" << endl;
  for (i = 1; i <= 3; i++) {
    cout << tran_begin[i] << endl;
  }
  cout << "tran_match:" << endl;
  for (i = 0; i <= env_number; i++) {
    for (j = 1; j <= 4; j++) {
      cout << setw(9) << tran_match[i][j];
    }
    cout << endl;
  }
  cout << "tran_x:" << endl;
  for (i = 0; i <= env_number; i++) {
    for (j = 1; j <= 3; j++) {
      cout << setw(9) << tran_x[i][j];
    }
    cout << endl;
  }
  cout << "tran_y:" << endl;
  for (i = 0; i <= env_number; i++) {
    for (j = 1; j <= 3; j++) {
      cout << setw(9) << tran_y[i][j];
    }
    cout << endl;
  }

  for (i = 0; i <= env_number; i++) {
    cout << "aa_pair " << i << endl;
    for (j = 1; j <= 20; j++) {
      for (k = 1; k <= 20; k++) {
        cout << setw(7) << aa_pair[i][j][k];
      }
      cout << endl;
    }
    cout << "aa_loop " << i << endl;
    for (j = 1; j <= 20; j++) {
      cout << aa_loop[i][j] << endl;
    }
    cout << "aa_bg1 " << i << endl;
    for (j = 1; j <= 20; j++) {
      cout << aa_bg1[i][j] << endl;
    }
    cout << "ss_pair" << i << endl;
    for (j = 1; j <= env_number; j++) {
      cout << ss_pair[i][j] << endl;
    }
    cout << "ss_loop" << i << endl;
    for (j = 1; j <= env_number; j++) {
      cout << ss_loop[i][j] << endl;
    }
  }
}

void hmm_psipred::set_parameters(hmm_psipred_parameters *pointer1) {
  p = pointer1;
  aa_pair = p->aa_pair;
  aa_loop = p->aa_loop;
  ss_pair = p->ss_pair;
  aa_pair1 = p->aa_pair1;
  aa_loop1 = p->aa_loop1;
  ss_pair1 = p->ss_pair1;
  ss_loop = p->ss_loop;
  ss_loop1 = p->ss_loop1;
  tran_begin = p->tran_begin;
  tran_match = p->tran_match;
  tran_x = p->tran_x;
  tran_y = p->tran_y;
  env_number = p->env_number;
  // cout << "env_number: " <<  env_number << endl;

  if (env_number == 9) {
    x_alphabet = x->ss->alphabet1;
    y_alphabet = y->ss->alphabet1;
  } else if (env_number == 3) {
    x_alphabet = x->ss->sstype;
    y_alphabet = y->ss->sstype;
  }
}

hmm_psipred::~hmm_psipred() {
  int i;

  if (Fm) free_gmatrix<float>(Fm, lenx, leny);
  if (Fx) free_gmatrix<float>(Fx, lenx, leny);
  if (Fy) free_gmatrix<float>(Fy, lenx, leny);

  if (Bm) free_gmatrix<float>(Bm, lenx, leny);
  if (Bx) free_gmatrix<float>(Bx, lenx, leny);
  if (By) free_gmatrix<float>(By, lenx, leny);
  if (probMat) {
    free_gmatrix<float>(probMat, lenx, leny);
  }
  if (score_matrix) free_gmatrix<float>(score_matrix, lenx, leny);
  if (score_bg_x) delete[] score_bg_x;
  if (score_bg_y) delete[] score_bg_y;
}

/*
void hmm_psipred::viterbi(int lor) {

        int i,j,k,m;
        int *ax = x->alignment[1];
        int *ay = y->alignment[1];
        vector<float> vv;
        vector<float>::iterator it;

        // allocation
        V = new float **[p->num_match_states +3];
        for(i=1;i<=p->num_match_states +2;i++) {
                V[i] = gmatrix<float>(lenx, leny);
        }
        T = new int **[p->num_match_states +3];
        for(i=1;i<=p->num_match_states +2;i++) {
                T[i] = gmatrix<int>(lenx, leny);
        }

        //cout<< "IN viterbi here =============" << endl;
        cout << "Viterbi algorithm:" << endl;

        for(i=0;i<=lenx;i++) {
           for(j=0;j<=leny;j++) {
                // match states
                //cout << "======== " << i << "  " << j << endl;
                for(k=1;k<=num_match_states;k++) {
                    if( (i==0) || (j==0) ) {
                        V[k][i][j] = LOG_ZERO;
                        T[k][i][j] = -1000;
                        if(debug>1) fprintf(stdout, "%d %d %d:  %7.3f    %d", k,
i, j, V[k][i][j], T[k][i][j]); if(debug>1) cout << endl; continue;
                    }
                    if( (i==1)&&(j==1) ) {
                        V[k][i][j] = p->t_log[0][k];
                        T[k][i][j] = 0;
                        if(lor==1) V[k][i][j] += p->m_lor[k][ax[i]][ay[j]];
                        else V[k][i][j] += p->m_log[k][ax[i]][ay[j]];
                        if(debug>1) fprintf(stdout, "%d %d %d:  %7.3f    %d", k,
i, j, V[k][i][j], T[k][i][j]); if(debug>1) cout << endl; continue;
                    }
                    vv.clear();
                    for(m=1;m<=num_match_states+2;m++) {
                        vv.push_back(p->t_log[m][k] + V[m][i-1][j-1]);
                    }
                    it = max_element(vv.begin(), vv.end());
                    if(lor==1) V[k][i][j]=p->m_lor[k][ax[i]][ay[j]]+(*it);
                    else V[k][i][j]=p->m_log[k][ax[i]][ay[j]]+(*it);
                    vv.erase(++it, vv.end());
                    T[k][i][j]=vv.size();
                    if(debug>1) fprintf(stdout, "%d %d %d:  %7.3f    %d", k, i,
j, V[k][i][j], T[k][i][j]); if(debug>1) cout << endl;
              }
                // insertion state 1
             //cout << "+++++++++++" << endl;
             if(i==0) { V[p1][i][j] = LOG_ZERO; }
             else {
                if( (j==0)&&(i==1) ) {
                        V[p1][i][j] = p->t_log[0][p1];
                        T[p1][i][j] = 0;
                        if(lor==1) V[p1][i][j] += p->indel_lor1[ax[i]];
                        else V[p1][i][j] += p->indel_log1[ax[i]];
                        if(debug>1) fprintf(stdout, "%d %d %d:  %7.3f    %d",
p1, i, j, V[p1][i][j], T[p1][i][j]); if(debug>1) cout << endl;
                }
                else {
                        vv.clear();
                        for(m=1;m<=num_match_states+2;m++) {
                                vv.push_back(p->t_log[m][p1] + V[m][i-1][j]);
                                //cout << p1 << " " << m << " " << i-1 << " " <<
j  << " " <<  p->t_log[m][p1] << " " << V[m][i-1][j] << endl;
                        }
                        it = max_element(vv.begin(), vv.end());
                        if(lor==1) V[p1][i][j]=p->indel_lor1[ax[i]] + (*it);
                        else V[p1][i][j]=p->indel_log1[ax[i]] + (*it);
                        vv.erase(++it, vv.end());
                        T[p1][i][j]=vv.size();
                        //fprintf(stdout, "%d %d %d:  %7.3f    %d", p1, i, j,
V[p1][i][j], T[p1][i][j]);
                        //cout << endl;
                }
             }
             //cout << "+++++++++++" << endl;

                // insertion state 2
             if(j==0) { V[p2][i][j] = LOG_ZERO; }
             else {
                if( (i==0)&&(j==1) ) {
                        V[p2][i][j] = p->t_log[0][p2];
                        if(lor==1) V[p2][i][j] += p->indel_lor2[ay[j]];
                        else V[p2][i][j] += p->indel_log2[ay[j]];
                        T[p2][i][j] = -1000;
                        if(debug>1) fprintf(stdout, "%d %d %d:  %7.3f    %d",
p2, i, j, V[p2][i][j], T[p2][i][j]); if(debug>1) cout << endl;
                }
                else {
                        vv.clear();
                        for(m=1;m<=num_match_states+2;m++) {
                                vv.push_back(p->t_log[m][p2]+V[m][i][j-1]);
                                //cout << m << " " << p2 << " " << i << " " <<
j-1  << " " <<  p->t_log[m][p2] << " " << V[m][i][j-1] << endl;
                        }
                        //for(it=vv.begin();it<vv.end();it++) { cout << *it << "
"; }
                        //cout << endl;
                        it=max_element(vv.begin(), vv.end());
                        if(lor==1) V[p2][i][j]=p->indel_lor2[ay[j]] + (*it);
                        else V[p2][i][j]=p->indel_log2[ay[j]] + (*it);
                        vv.erase(++it, vv.end());
                        T[p2][i][j]=vv.size();
                        if(debug>1) fprintf(stdout, "%d %d %d:  %7.3f    %d",
p2, i, j, V[p2][i][j], T[p2][i][j]); if(debug>1) cout << endl;
                }
            }
          }
        }
        if(debug>1) cout<< "IN viterbi here =============" << endl;

        // the ending state
        vv.clear();
        for(m=1;m<=num_match_states+2;m++) {
                vv.push_back(p->t_log[m][p_end]+V[m][lenx][leny]);
                if(debug>1) cout << m << " " << V[m][lenx][leny] << " " <<
p->t_log[m][p_end] << " " << V[m][lenx][leny]+p->t_log[m][p_end] << endl;
        }
        it=max_element(vv.begin(), vv.end());
        VE = (*it);
        vv.erase(++it, vv.end());
        TE = vv.size();
        vv.clear();

        if(debug>1) cout << "VE and TE: " << VE << "  " << TE << endl;

        // trace back
        int *trace = new int[lenx+leny+2];
        trace[0] = TE;
        int px = lenx;
        int py = leny;
        //cout << "Tracing result: " << endl;
        for(i=1;i<=lenx+leny;i++) {
                trace[i] = T[trace[i-1]][px][py];
                if (trace[i]==-1000) { break;}
                if(trace[i-1]<=num_match_states) {
                        px--;
                        py--;
                }
                else if(trace[i-1]==num_match_states+1) {
                        px--;
                }
                else if(trace[i-1]==num_match_states+2) {
                        py--;
                }
                if(debug>1) cout << i << "  " << trace[i] << endl;
        }
        int trace_len = i;
        string sx, sy;
        sx = "";
        sy = "";
        k = 0; m = 0;
        for(i=0;i<trace_len;i++) {
                j = trace_len - i-1;
                if(trace[j]<=num_match_states) {
                        sx.append(1, x->aseq[0][k]);
                        sy.append(1, y->aseq[0][m]);
                        k++;
                        m++;
                }
                if(trace[j]==num_match_states+1) {
                        sx.append(1, x->aseq[0][k]);
                        sy.append(1, '-');
                        k++;
                }
                if(trace[j]==num_match_states+2) {
                        sx.append(1, '-');
                        sy.append(1, y->aseq[0][m]);
                        m++;
                }
        }
        cout << left << setw(20) << x->aname[0] << "  " << sx << endl;
        cout << left << setw(20) << y->aname[0] << "  " << sy << endl;
        cout << endl;

        // clear up
        for(i=1;i<=p->num_match_states +2;i++) {
                free_gmatrix<float>(V[i], lenx, leny);
                free_gmatrix<int>(T[i], lenx, leny);
        }
        delete [] V;
        delete [] T;
        delete [] trace;
}

*/

// added score_shift
void hmm_psipred::get_scores(double ss_weight) {
  int i, j, k, l;

  float tmp_aa_score;
  float tmp_ss_score;

  score_matrix = gmatrix<float>(lenx, leny);
  score_bg_x = gvector<float>(lenx);
  score_bg_y = gvector<float>(leny);

  if (!x->done_score_bg2) {
    x->get_score_bg_mine2(aa_loop[0], ss_loop1[0], env_number);
  }
  if (!y->done_score_bg2) {
    y->get_score_bg_mine2(aa_loop[0], ss_loop1[0], env_number);
  }
  if (debug > 1) cout << x->done_score_bg << " " << y->done_score_bg << endl;

  if (debug > 1) cout << "score_bg_x" << endl;
  for (i = 1; i <= lenx; i++) {
    score_bg_x[i] = x->score_bg_aa[i] + ss_weight * x->score_bg_ss[i];
    if (debug > 1)
      cout << i << " " << x->score_bg_aa[i] << " " << x->score_bg_ss[i] << " "
           << score_bg_x[i] << endl;
  }
  if (debug > 1) cout << "score_bg_y" << endl;
  for (i = 1; i <= leny; i++) {
    score_bg_y[i] = y->score_bg_aa[i] + ss_weight * y->score_bg_ss[i];
    if (debug > 1)
      cout << i << " " << y->score_bg_aa[i] << " " << y->score_bg_ss[i] << " "
           << score_bg_y[i] << endl;
  }

  if (debug > 1) cout << "score_matrix" << endl;
  for (i = 1; i <= lenx; i++) {
    for (j = 1; j <= leny; j++) {
      tmp_aa_score = 0;
      for (k = 1; k <= 20; k++) {
        tmp_aa_score += x->prof_effn[i][k] * y->prof_freq[j][k];
      }
      score_matrix[i][j] = tmp_aa_score / x->prof_sum_eff[i];
      tmp_aa_score = 0;
      for (k = 1; k <= 20; k++) {
        tmp_aa_score += y->prof_effn[j][k] * x->prof_freq[i][k];
      }
      score_matrix[i][j] =
          score_matrix[i][j] + tmp_aa_score / y->prof_sum_eff[j];
      /*
      tmp_aa_score = 0;
      for(k=1;k<=20;k++) {
              if(!x->prof_effn[i][k]) continue; // skipping residues with zero
      counts for(l=1;l<=20;l++) { if(!y->prof_effn[j][l]) continue; // skipping
      residues with zero counts tmp_aa_score += x->prof_effn[i][k] *
      y->prof_effn[j][l] * aa_pair1[x_alphabet[i]][k][l]; // sum-of-pairs
      measure
              }
      }
      score_matrix[i][j] =
      log(tmp_aa_score/x->prof_sum_eff[i]/y->prof_sum_eff[j]); // divided by
      sum_eff
      */
      tmp_ss_score = log(
          ss_pair1[x_alphabet[i]][y_alphabet[j]]);  // the probability of
                                                    // emitting a LE pair in the
      // environment of the position of the first sequence
      // LE: local environment; for three secondary structure types, a pair is
      // [HEC][HEC]

      if (debug > 1)
        cout << i << " " << j << " " << score_matrix[i][j] << " "
             << tmp_ss_score << " ";
      score_matrix[i][j] += ss_weight * tmp_ss_score;
      if (debug > 1) cout << score_matrix[i][j] << endl;
    }
  }
}

void hmm_psipred::get_scores_sum_of_pairs(double ss_weight) {
  int i, j, k, l;

  float tmp_aa_score;
  float tmp_ss_score;

  score_matrix = gmatrix<float>(lenx, leny);
  score_bg_x = gvector<float>(lenx);
  score_bg_y = gvector<float>(leny);

  if (!x->done_score_bg) {
    x->get_score_bg_mine(aa_loop1[0], ss_loop1[0], env_number);
  }
  if (!y->done_score_bg) {
    y->get_score_bg_mine(aa_loop1[0], ss_loop1[0], env_number);
  }
  if (debug > 1) cout << x->done_score_bg << " " << y->done_score_bg << endl;

  if (debug > 1) cout << "score_bg_x" << endl;
  for (i = 1; i <= lenx; i++) {
    score_bg_x[i] = x->score_bg_aa[i] + ss_weight * x->score_bg_ss[i];
    if (debug > 1)
      cout << i << " " << x->score_bg_aa[i] << " " << x->score_bg_ss[i] << " "
           << score_bg_x[i] << endl;
  }
  if (debug > 1) cout << "score_bg_y" << endl;
  for (i = 1; i <= leny; i++) {
    score_bg_y[i] = y->score_bg_aa[i] + ss_weight * y->score_bg_ss[i];
    if (debug > 1)
      cout << i << " " << y->score_bg_aa[i] << " " << y->score_bg_ss[i] << " "
           << score_bg_y[i] << endl;
  }

  if (debug > 1) cout << "score_matrix" << endl;
  for (i = 1; i <= lenx; i++) {
    for (j = 1; j <= leny; j++) {
      tmp_aa_score = 0;
      for (k = 1; k <= 20; k++) {
        if (!x->prof_effn[i][k])
          continue;  // skipping residues with zero counts
        for (l = 1; l <= 20; l++) {
          if (!y->prof_effn[j][l])
            continue;  // skipping residues with zero counts
          tmp_aa_score +=
              x->prof_effn[i][k] * y->prof_effn[j][l] *
              aa_pair1[x_alphabet[i]][k][l];  // sum-of-pairs measure
        }
      }
      score_matrix[i][j] = log(tmp_aa_score / x->prof_sum_eff[i] /
                               y->prof_sum_eff[j]);  // divided by sum_eff

      tmp_ss_score = log(
          ss_pair1[x_alphabet[i]][y_alphabet[j]]);  // the probability of
                                                    // emitting a LE pair in the
      // environment of the position of the first sequence
      // LE: local environment; for three secondary structure types, a pair is
      // [HEC][HEC]

      if (debug > 1)
        cout << i << " " << j << " " << score_matrix[i][j] << " "
             << tmp_ss_score << " ";
      score_matrix[i][j] += ss_weight * tmp_ss_score;
      if (debug > 1) cout << score_matrix[i][j] << endl;
    }
  }
}

void hmm_psipred::get_scores(double ss_weight, double score_weight) {
  int i, j, k, l;

  float tmp_aa_score;
  float tmp_ss_score;

  score_matrix = gmatrix<float>(lenx, leny);
  score_bg_x = gvector<float>(lenx);
  score_bg_y = gvector<float>(leny);

  // aa_loop: with logarithm taken; ss_loop1: raw frequencies
  if (!x->done_score_bg) {
    x->get_score_bg_mine2(aa_loop[0], ss_loop1[0], env_number);
  }
  if (!y->done_score_bg) {
    y->get_score_bg_mine2(aa_loop[0], ss_loop1[0], env_number);
  }
  if (debug > 1) cout << x->done_score_bg << " " << y->done_score_bg << endl;

  if (debug > 1) cout << "score_bg_x" << endl;
  for (i = 1; i <= lenx; i++) {
    score_bg_x[i] =
        score_weight * x->score_bg_aa[i] + ss_weight * x->score_bg_ss[i];
    if (debug > 1)
      cout << i << " " << x->score_bg_aa[i] << " " << x->score_bg_ss[i] << " "
           << score_bg_x[i] << endl;
  }
  if (debug > 1) cout << "score_bg_y" << endl;
  for (i = 1; i <= leny; i++) {
    score_bg_y[i] =
        score_weight * y->score_bg_aa[i] + ss_weight * y->score_bg_ss[i];
    if (debug > 1)
      cout << i << " " << y->score_bg_aa[i] << " " << y->score_bg_ss[i] << " "
           << score_bg_y[i] << endl;
  }

  if (debug > 1) cout << "score_matrix" << endl;
  for (i = 1; i <= lenx; i++) {
    for (j = 1; j <= leny; j++) {
      tmp_aa_score = 0;
      for (k = 1; k <= 20; k++) {
        // if((i==233)&&(j==202)) { cout << "k" << " " << k << " " <<
        // x->prof_effn[i][k] << " " << y->prof_freq[j][k] << endl; }
        tmp_aa_score += x->prof_effn[i][k] * y->prof_freq[j][k];
      }
      score_matrix[i][j] = tmp_aa_score / x->prof_sum_eff[i];
      tmp_aa_score = 0;
      for (k = 1; k <= 20; k++) {
        // if((i==233)&&(j==202)) { cout << k << " " << y->prof_effn[j][k] << "
        // " << x->prof_freq[i][k] << endl; }
        tmp_aa_score += y->prof_effn[j][k] * x->prof_freq[i][k];
      }
      score_matrix[i][j] =
          score_matrix[i][j] + tmp_aa_score / y->prof_sum_eff[j];
      /*
      for(k=1;k<=20;k++) {
              if(!x->prof_effn[i][k]) continue; // skipping residues with zero
      counts for(l=1;l<=20;l++) { if(!y->prof_effn[j][l]) continue; // skipping
      residues with zero counts tmp_aa_score += x->prof_effn[i][k] *
      y->prof_effn[j][l] * aa_pair1[x_alphabet[i]][k][l]; // sum-of-pairs
      measure
              }
      }
      score_matrix[i][j] =
      log(tmp_aa_score/x->prof_sum_eff[i]/y->prof_sum_eff[j]); // divided by
      sum_eff
      */

      tmp_ss_score = log(
          ss_pair1[x_alphabet[i]][y_alphabet[j]]);  // the probability of
                                                    // emitting a LE pair in the
      // environment of the position of the first sequence
      // LE: local environment; for three secondary structure types, a pair is
      // [HEC][HEC]

      if (debug > 1)
        cout << i << " " << j << " "
             << " " << x->prof_sum_eff[i] << " " << y->prof_sum_eff[j] << " "
             << score_matrix[i][j] << " " << tmp_ss_score << " ";
      score_matrix[i][j] =
          score_weight * score_matrix[i][j] + ss_weight * tmp_ss_score;
      // score_matrix[i][j] += ss_weight * tmp_ss_score;
      // score_matrix[i][j] *= score_weight;
      score_matrix[i][j] += score_shift;
      if (debug > 1) cout << score_matrix[i][j] << endl;
    }
  }
}

void hmm_psipred::get_scores_sum_of_pairs(double ss_weight,
                                          double score_weight) {
  int i, j, k, l;

  float tmp_aa_score;
  float tmp_ss_score;

  score_matrix = gmatrix<float>(lenx, leny);
  score_bg_x = gvector<float>(lenx);
  score_bg_y = gvector<float>(leny);

  if (!x->done_score_bg) {
    x->get_score_bg_mine(aa_loop1[0], ss_loop1[0], env_number);
  }
  if (!y->done_score_bg) {
    y->get_score_bg_mine(aa_loop1[0], ss_loop1[0], env_number);
  }
  if (debug > 1) cout << x->done_score_bg << " " << y->done_score_bg << endl;

  if (debug > 1) cout << "score_bg_x" << endl;
  for (i = 1; i <= lenx; i++) {
    score_bg_x[i] =
        (x->score_bg_aa[i] + ss_weight * x->score_bg_ss[i]) * score_weight;
    if (debug > 1)
      cout << i << " " << x->score_bg_aa[i] << " " << x->score_bg_ss[i] << " "
           << score_bg_x[i] << endl;
  }
  if (debug > 1) cout << "score_bg_y" << endl;
  for (i = 1; i <= leny; i++) {
    score_bg_y[i] =
        (y->score_bg_aa[i] + ss_weight * y->score_bg_ss[i]) * score_weight;
    if (debug > 1)
      cout << i << " " << y->score_bg_aa[i] << " " << y->score_bg_ss[i] << " "
           << score_bg_y[i] << endl;
  }

  if (debug > 1) cout << "score_matrix" << endl;
  for (i = 1; i <= lenx; i++) {
    for (j = 1; j <= leny; j++) {
      tmp_aa_score = 0;
      for (k = 1; k <= 20; k++) {
        if (!x->prof_effn[i][k])
          continue;  // skipping residues with zero counts
        for (l = 1; l <= 20; l++) {
          if (!y->prof_effn[j][l])
            continue;  // skipping residues with zero counts
          tmp_aa_score +=
              x->prof_effn[i][k] * y->prof_effn[j][l] *
              aa_pair1[x_alphabet[i]][k][l];  // sum-of-pairs measure
        }
      }
      score_matrix[i][j] = log(tmp_aa_score / x->prof_sum_eff[i] /
                               y->prof_sum_eff[j]);  // divided by sum_eff

      tmp_ss_score = log(
          ss_pair1[x_alphabet[i]][y_alphabet[j]]);  // the probability of
                                                    // emitting a LE pair in the
      // environment of the position of the first sequence
      // LE: local environment; for three secondary structure types, a pair is
      // [HEC][HEC]

      if (debug > 1)
        cout << i << " " << j << " " << score_matrix[i][j] << " "
             << tmp_ss_score << " ";
      score_matrix[i][j] += ss_weight * tmp_ss_score;
      score_matrix[i][j] *= score_weight;
      if (debug > 1) cout << score_matrix[i][j] << endl;
    }
  }
}

void hmm_psipred::forward() {
  int i, j, k, m;
  int *ax = x->alignment[1];
  int *ay = y->alignment[1];

  int alphax, alphay;
  int alphax_1, alphay_1;
  int axi, ayj;

  Fm = gmatrix<float>(lenx, leny);
  Fx = gmatrix<float>(lenx, leny);
  Fy = gmatrix<float>(lenx, leny);

  if (debug > 1) cout << "Forward algorithm:" << endl;

  for (i = 0; i <= lenx; i++) {
    for (j = 0; j <= leny; j++) {
      if ((i == 0) && (j == 0)) {
        Fm[i][j] = Fx[i][j] = Fy[i][j] = LOG_ZERO;
        continue;
      }

      if ((i == 0) && (j == 1)) {
        // cout << tran_begin[3]<<" "<<y_alphabet[j]<<" "<< ay[j]<<" " <<
        // aa_loop[y_alphabet[j]][ay[j]]<<endl; Fy[i][j] = tran_begin[3] +
        // aa_loop[y_alphabet[j]][ay[j]];
        Fy[i][j] = tran_begin[3] + aa_loop[y_alphabet[j]][ay[j]];
        // cout << i << " " << j << " " << y->aseq[0][j-1] << " " << Fy[i][j] <<
        // endl;
        Fx[i][j] = Fm[i][j] = LOG_ZERO;
        continue;
      }
      if ((i == 0) && (j > 1)) {
        // cout << "*****" << Fy[i][j-1] << " " << tran_y[y_alphabet[j-1]][2] <<
        // " " << aa_loop[y_alphabet[j]][ay[j]] << endl; Fy[i][j] = Fy[i][j-1] +
        // tran_y[y_alphabet[j-1]][2] + aa_loop[y_alphabet[j]][ay[j]]; cout <<
        // y_alphabet[j] << " " << ay[j]  << endl;
        Fy[i][j] = Fy[i][j - 1] + tran_y[y_alphabet[j - 1]][2] +
                   aa_loop[y_alphabet[j]][ay[j]];
        Fx[i][j] = Fm[i][j] = LOG_ZERO;
        // cout << i << " " << j << " " << y->aseq[0][j-1] << " " << Fy[i][j] <<
        // endl;
        continue;
      }
      if ((i == 1) && (j == 0)) {
        // cout << tran_begin[2] << " " << aa_loop[x_alphabet[i]][ax[i]] <<
        // endl; Fx[i][j] = tran_begin[2] + aa_loop[x_alphabet[i]][ax[i]];
        Fx[i][j] = tran_begin[2] + aa_loop[x_alphabet[i]][ax[i]];
        Fy[i][j] = Fm[i][j] = LOG_ZERO;
        // cout << i << " " << j << " " << Fx[i][j] << endl;
        continue;
      }
      if ((i > 1) && (j == 0)) {
        // Fx[i][j] = Fx[i-1][j] + tran_x[x_alphabet[i-1]][2] +
        // aa_loop[x_alphabet[i]][ax[i]];
        Fx[i][j] = Fx[i - 1][j] + tran_x[x_alphabet[i - 1]][2] +
                   aa_loop[x_alphabet[i]][ax[i]];
        // cout << i << " " << j << " " << Fx[i][j] << endl;
        Fy[i][j] = Fm[i][j] = LOG_ZERO;
        continue;
      }
      // alphax = x_alphabet[i];
      // alphay = y_alphabet[j];
      alphax = x_alphabet[i];
      alphay = y_alphabet[j];
      axi = ax[i];
      ayj = ay[j];
      if ((i == 1) && (j == 1)) {
        Fm[i][j] = tran_begin[1] + aa_pair[alphax][axi][ayj] +
                   ss_w * ss_pair[alphax][alphay];
        // cout << i << " " << j << " " << Fm[i][j] << endl;
        Fx[i][j] = Fy[i][j] = LOG_ZERO;
        continue;
      }

      // alphax_1 = x_alphabet[i-1];
      // alphay_1 = y_alphabet[j-1];
      alphax_1 = x_alphabet[i - 1];
      alphay_1 = y_alphabet[j - 1];
      if (i == 1) alphax_1 = 1;
      if (j == 1) alphay_1 = 3;
      Fm[i][j] = LOG_ADD(tran_match[alphax_1][1] + Fm[i - 1][j - 1],
                         tran_x[alphax_1][1] + Fx[i - 1][j - 1],
                         tran_y[alphay_1][1] + Fy[i - 1][j - 1]) +
                 p->aa_pair[alphax][axi][ayj] +
                 ss_w * p->ss_pair[alphax][alphay];
      Fx[i][j] = LOG_ADD(tran_match[alphax_1][2] + Fm[i - 1][j],
                         tran_x[alphax_1][2] + Fx[i - 1][j]) +
                 p->aa_loop[alphax][axi];
      Fy[i][j] = LOG_ADD(tran_match[alphax][3] + Fm[i][j - 1],
                         tran_y[alphay_1][2] + Fy[i][j - 1]) +
                 p->aa_loop[alphay][ayj];

      // cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " " <<
      // Fy[i][j] << endl; cout << "Fy[i][j]: " << i << " " << j << " " <<
      // tran_match[alphax_1][3] << " " << Fm[i][j-1] << " " <<
      // tran_y[alphay_1][2] << " " << Fy[i][j-1] << " " <<
      // p->aa_loop[alphay][ayj] << endl;
    }
  }

  // the ending state
  // FE = LOG_ADD(Fm[lenx][leny]+tran_match[x_alphabet[lenx]][4],
  // Fx[lenx][leny]+tran_x[x_alphabet[lenx]][3],
  // Fy[lenx][leny]+tran_y[y_alphabet[leny]][3]);
  FE = LOG_ADD(Fm[lenx][leny] + tran_match[x_alphabet[lenx]][4],
               Fx[lenx][leny] + tran_x[x_alphabet[lenx]][3],
               Fy[lenx][leny] + tran_y[y_alphabet[leny]][3]);
  // cout << Fm[lenx][leny]+tran_match[x_alphabet[lenx]][4] << " " <<
  // Fx[lenx][leny]+tran_x[x_alphabet[lenx]][3] << " " <<
  // Fy[lenx][leny]+tran_y[y_alphabet[leny]][3] << endl;

  // cout << "tran_y[y_alphabet[leny]][3]: " << tran_y[y_alphabet[leny]][3] <<
  // endl;

  if (debug > 1) fprintf(stdout, "FE: %f\t", FE);
  fflush(stdout);  // cout << "FE: " << FE << endl;
}

void hmm_psipred::backward() {
  int i, j, k, m;
  int *ax = x->alignment[1];
  int *ay = y->alignment[1];

  int alphax, alphay;
  int alphax_1, alphay_1;
  int axi_1, ayj_1;

  Bm = gmatrix<float>(lenx, leny);
  Bx = gmatrix<float>(lenx, leny);
  By = gmatrix<float>(lenx, leny);

  for (i = lenx; i >= 0; i--) {
    for (j = leny; j >= 0; j--) {
      // match states and insertion states together
      // alphax = x_alphabet[i];
      // alphay = y_alphabet[j];
      alphax = x_alphabet[i];
      alphay = y_alphabet[j];
      if ((i == lenx) && (j == leny)) {
        Bm[i][j] = tran_match[alphax][4];
        Bx[i][j] = tran_x[alphax][3];
        By[i][j] = tran_y[alphay][3];
        // cout << "By[i][j]: " << i << " " << j << " " << By[i][j] << endl;
        continue;
      }
      if (i == lenx) {
        // By[i][j] =
        // By[i][j+1]+tran_y[alphay][2]+aa_loop[y_alphabet[j+1]][ay[j+1]];
        By[i][j] = By[i][j + 1] + tran_y[alphay][2] +
                   aa_loop[y_alphabet[j + 1]][ay[j + 1]];
        // Bm[i][j] =
        // By[i][j+1]+tran_match[alphax][3]+aa_loop[y_alphabet[j+1]][ay[j+1]];
        Bm[i][j] = By[i][j + 1] + tran_match[alphax][3] +
                   aa_loop[y_alphabet[j + 1]][ay[j + 1]];
        // cout << "By[i][j]: " << i << " " << j << " " << By[i][j] << endl;
        // cout << "Bm[i][j]: " << i << " " << j << " " << Bm[i][j] << endl;
        Bx[i][j] = LOG_ZERO;
        continue;
      }
      if (j == leny) {
        // Bx[i][j] =
        // Bx[i+1][j]+tran_x[alphax][2]+aa_loop[x_alphabet[i+1]][ax[i+1]];
        // Bm[i][j] =
        // Bx[i+1][j]+tran_match[alphax][2]+aa_loop[x_alphabet[i+1]][ax[i+1]];
        Bx[i][j] = Bx[i + 1][j] + tran_x[alphax][2] +
                   aa_loop[x_alphabet[i + 1]][ax[i + 1]];
        Bm[i][j] = Bx[i + 1][j] + tran_match[alphax][2] +
                   aa_loop[x_alphabet[i + 1]][ax[i + 1]];
        By[i][j] = LOG_ZERO;
        // cout << "Bx[i][j]: " << i << " " << j << " " << Bx[i][j] << endl;
        // cout << "Bm[i][j]: " << i << " " << j << " " << Bm[i][j] << endl;
        continue;
      }
      // alphax_1 = x_alphabet[i+1];
      // alphay_1 = y_alphabet[j+1];
      alphax_1 = x_alphabet[i + 1];
      alphay_1 = y_alphabet[j + 1];
      axi_1 = ax[i + 1];
      ayj_1 = ay[j + 1];

      Bm[i][j] = LOG_ADD(
          Bm[i + 1][j + 1] + tran_match[alphax][1] +
              aa_pair[alphax_1][axi_1][ayj_1] +
              ss_w * ss_pair[alphax_1][alphay_1],
          Bx[i + 1][j] + tran_match[alphax][2] + aa_loop[alphax_1][axi_1],
          By[i][j + 1] + tran_match[alphax][3] + aa_loop[alphay_1][ayj_1]);
      // cout << "Bm[i][j]: " << i << " " << j << " " << Bm[i][j] << endl;
      // cout << setprecision(10) <<
      // Bm[i+1][j+1]+tran_x[alphax][1]+aa_pair[alphax_1][axi_1][ayj_1]+ ss_w *
      // ss_pair[alphax_1][alphay_1] << " " <<
      // Bx[i+1][j]+tran_x[alphax][2]+aa_loop[alphax_1][axi_1] << endl;
      // if(Bm[i+1][j+1]+tran_x[alphax][1]+aa_pair[alphax_1][axi_1][ayj_1]+ss_w
      // * ss_pair[alphax_1][alphay_1]==Bx[i+1][j]+tran_x[alphax][2]+aa_loop[alphax_1][axi_1]) { cout << "x equal y" << endl; }
      Bx[i][j] =
          LOG_ADD(Bm[i + 1][j + 1] + tran_x[alphax][1] +
                      aa_pair[alphax_1][axi_1][ayj_1] +
                      ss_w * ss_pair[alphax_1][alphay_1],
                  Bx[i + 1][j] + tran_x[alphax][2] + aa_loop[alphax_1][axi_1]);
      // cout << "Bx[i][j]: " << i << " " << j << " " << Bx[i][j] << endl;
      By[i][j] =
          LOG_ADD(Bm[i + 1][j + 1] + tran_y[alphay][1] +
                      aa_pair[alphax_1][axi_1][ayj_1] +
                      ss_w * ss_pair[alphax_1][alphay_1],
                  By[i][j + 1] + tran_y[alphay][2] + aa_loop[alphay_1][ayj_1]);
      // cout << "By[i][j]: " << i << " " << j << " " << By[i][j] << endl;
    }
  }

  // cout << tran_begin[1] << " " << tran_match[x_alphabet[1]][3] << " " <<
  // tran_y[y_alphabet[2]][3] << " " << aa_pair[x_alphabet[1]][ax[1]][ay[1]] << "
  // " << ss_pair[x_alphabet[1]][y_alphabet[1]] << " " <<
  // aa_loop[y_alphabet[2]][ay[2]] << endl;

  // Bm[0][0] =
  // Bm[1][1]+tran_begin[1]+aa_pair[x_alphabet[1]][ax[1]][ay[1]]+ss_pair[x_alphabet[1]][y_alphabet[1]];
  // Bx[0][0] = Bx[1][0]+tran_begin[2]+aa_loop[x_alphabet[1]][ax[1]];
  // By[0][0] = By[0][1]+tran_begin[3]+aa_loop[y_alphabet[1]][ay[1]];
  Bm[0][0] = Bm[1][1] + tran_begin[1] + aa_pair[x_alphabet[1]][ax[1]][ay[1]] +
             ss_pair[x_alphabet[1]][y_alphabet[1]];
  Bx[0][0] = Bx[1][0] + tran_begin[2] + aa_loop[x_alphabet[1]][ax[1]];
  By[0][0] = By[0][1] + tran_begin[3] + aa_loop[y_alphabet[1]][ay[1]];

  // BE =
  // LOG_ADD(Bm[1][1]+tran_begin[1]+aa_pair[x_alphabet[1]][ax[1]][ay[1]]+ss_pair[x_alphabet[1]][y_alphabet[1]],
  // Bx[1][0]+tran_begin[2]+aa_loop[x_alphabet[1]][ax[1]],
  // By[0][1]+tran_begin[3]+aa_loop[y_alphabet[1]][ay[1]]);
  BE = LOG_ADD(Bm[1][1] + tran_begin[1] + aa_pair[x_alphabet[1]][ax[1]][ay[1]] +
                   ss_w * ss_pair[x_alphabet[1]][y_alphabet[1]],
               Bx[1][0] + tran_begin[2] + aa_loop[x_alphabet[1]][ax[1]],
               By[0][1] + tran_begin[3] + aa_loop[y_alphabet[1]][ay[1]]);

  // cout <<
  // Bm[1][1]+tran_begin[1]+aa_pair[x_alphabet[1]][ax[1]][ay[1]]+ss_pair[x_alphabet[1]][y_alphabet[1]]
  // << " " << Bx[1][0]+tran_begin[2]+aa_loop[x_alphabet[1]][ax[1]] << " " <<
  // By[0][1]+tran_begin[3]+aa_loop[y_alphabet[1]][ay[1]] << endl;

  if (debug > 1) fprintf(stdout, "BE: %f\n", BE);

  // ScoreType full_prob;
  if (debug > 1)
    cout << "Obtain posterior probabilities (>0.01): " << x->aname[0] << "\t"
         << y->aname[0] << endl;
  probMat = gmatrix<float>(lenx, leny);
  for (i = 1; i <= lenx; i++) {
    for (j = 1; j <= leny; j++) {
      probMat[i][j] = exp(Fm[i][j] + Bm[i][j] - FE);
      if (debug > 1)
        cout << "probMat: " << i << " " << j << " " << setprecision(7)
             << probMat[i][j] << endl;
    }
  }

  // fprintf(stdout, "================\n");
}

void hmm_psipred::forward1() {
  int i, j, k, m;
  int *ax = x->alignment[1];
  int *ay = y->alignment[1];

  int alphax, alphay;
  int alphax_1, alphay_1;
  int axi, ayj;

  Fm = gmatrix<float>(lenx, leny);
  Fx = gmatrix<float>(lenx, leny);
  Fy = gmatrix<float>(lenx, leny);

  if (debug > 1) cout << "Forward algorithm:" << endl;

  cout << setprecision(10);

  for (i = 0; i <= lenx; i++) {
    for (j = 0; j <= leny; j++) {
      if ((i == 0) && (j == 0)) {
        Fm[i][j] = Fx[i][j] = Fy[i][j] = LOG_ZERO;
        if (debug > 1)
          cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
               << Fy[i][j] << endl;
        continue;
      }

      if ((i == 0) && (j == 1)) {
        // cout << tran_begin[3]<<" "<<y_alphabet[j]<<" "<< ay[j]<<" " <<
        // aa_loop[y_alphabet[j]][ay[j]]<<endl; Fy[i][j] = tran_begin[3] +
        // aa_loop[y_alphabet[j]][ay[j]];
        Fy[i][j] = tran_begin[3] + score_bg_y[1];
        // cout << i << " " << j << " " << Fy[i][j] << endl;
        Fx[i][j] = Fm[i][j] = LOG_ZERO;
        if (debug > 1)
          cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
               << Fy[i][j] << endl;
        continue;
      }
      if ((i == 0) && (j > 1)) {
        // cout << Fy[i][j-1] << " " << tran_y[y_alphabet[j-1]][2] << " " <<
        // aa_loop[y_alphabet[j]][ay[j]] << endl; Fy[i][j] = Fy[i][j-1] +
        // tran_y[y_alphabet[j-1]][2] + aa_loop[y_alphabet[j]][ay[j]];
        Fy[i][j] = Fy[i][j - 1] + tran_y[0][2] + score_bg_y[j];
        Fx[i][j] = Fm[i][j] = LOG_ZERO;
        // cout << i << " " << j << " " << Fy[i][j] << endl;
        if (debug > 1)
          cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
               << Fy[i][j] << endl;
        continue;
      }
      if ((i == 1) && (j == 0)) {
        // cout << tran_begin[2] << " " << aa_loop[x_alphabet[i]][ax[i]] <<
        // endl; Fx[i][j] = tran_begin[2] + aa_loop[x_alphabet[i]][ax[i]];
        Fx[i][j] = tran_begin[2] + score_bg_x[1];
        Fy[i][j] = Fm[i][j] = LOG_ZERO;
        // cout << i << " " << j << " " << Fx[i][j] << endl;
        if (debug > 1)
          cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
               << Fy[i][j] << endl;
        continue;
      }
      if ((i > 1) && (j == 0)) {
        // Fx[i][j] = Fx[i-1][j] + tran_x[x_alphabet[i-1]][2] +
        // aa_loop[x_alphabet[i]][ax[i]];
        Fx[i][j] = Fx[i - 1][j] + tran_x[0][2] + score_bg_x[i];
        // cout << i << " " << j << " " << Fx[i][j] << endl;
        Fy[i][j] = Fm[i][j] = LOG_ZERO;
        if (debug > 1)
          cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
               << Fy[i][j] << endl;
        continue;
      }
      alphax = x_alphabet[i];
      alphay = y_alphabet[j];
      axi = ax[i];
      ayj = ay[j];
      if ((i == 1) && (j == 1)) {
        // Fm[i][j] = tran_begin[1] + aa_pair[alphax][axi][ayj] +
        // ss_pair[alphax][alphay];
        Fm[i][j] = tran_begin[1] + score_matrix[1][1];
        // cout << i << " " << j << " " << Fm[i][j] << endl;
        Fx[i][j] = Fy[i][j] = LOG_ZERO;
        if (debug > 1)
          cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
               << Fy[i][j] << endl;
        continue;
      }

      alphax_1 = x_alphabet[i - 1];
      alphay_1 = y_alphabet[j - 1];
      if (i == 1) alphax_1 = 1;
      if (j == 1) alphay_1 = 3;
      // Fm[i][j] = LOG_ADD(tran_match[alphax_1][1]+Fm[i-1][j-1],
      // tran_x[0][1]+Fx[i-1][j-1], tran_y[0][1]+Fy[i-1][j-1]) +
      // p->aa_pair[alphax][axi][ayj] + p->ss_pair[alphax][alphay]; Fx[i][j] =
      // LOG_ADD(tran_match[alphax_1][2]+Fm[i-1][j], tran_x[0][2]+Fx[i-1][j]) +
      // p->aa_loop[alphax][axi]; Fy[i][j] =
      // LOG_ADD(tran_match[alphax][3]+Fm[i][j-1],
      // tran_y[alphay_1][2]+Fy[i][j-1]) + p->aa_loop[alphay][ayj];
      Fm[i][j] = LOG_ADD(tran_match[alphax_1][1] + Fm[i - 1][j - 1],
                         tran_x[0][1] + Fx[i - 1][j - 1],
                         tran_y[0][1] + Fy[i - 1][j - 1]) +
                 score_matrix[i][j];
      Fx[i][j] = LOG_ADD(tran_match[alphax_1][2] + Fm[i - 1][j],
                         tran_x[0][2] + Fx[i - 1][j]) +
                 score_bg_x[i];
      Fy[i][j] = LOG_ADD(tran_match[alphax][3] + Fm[i][j - 1],
                         tran_y[0][2] + Fy[i][j - 1]) +
                 score_bg_y[j];

      if (debug > 1)
        cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
             << Fy[i][j] << endl;
      // cout << "Fy[i][j]: " << i << " " << j << " " << tran_match[alphax_1][3]
      // << " " << Fm[i][j-1] << " " << tran_y[alphay_1][2] << " " << Fy[i][j-1]
      // << " " << p->aa_loop[alphay][ayj] << endl;
    }
  }

  // the ending state
  // FE = LOG_ADD(Fm[lenx][leny]+tran_match[x_alphabet[lenx]][4],
  // Fx[lenx][leny]+tran_x[x_alphabet[lenx]][3],
  // Fy[lenx][leny]+tran_y[y_alphabet[leny]][3]);
  FE = LOG_ADD(Fm[lenx][leny] + tran_match[x_alphabet[lenx]][4],
               Fx[lenx][leny] + tran_x[0][3], Fy[lenx][leny] + tran_y[0][3]);
  // cout << Fm[lenx][leny]+tran_match[x_alphabet[lenx]][4] << " " <<
  // Fx[lenx][leny]+tran_x[x_alphabet[lenx]][3] << " " <<
  // Fy[lenx][leny]+tran_y[y_alphabet[leny]][3] << endl;

  // cout << "tran_y[y_alphabet[leny]][3]: " << tran_y[y_alphabet[leny]][3] <<
  // endl;

  if (debug > 1)
    fprintf(stdout, "FE: %f\n", FE);  // cout << "FE: " << FE << endl;
}

void hmm_psipred::backward1() {
  int i, j, k, m;
  int *ax = x->alignment[1];
  int *ay = y->alignment[1];

  int alphax, alphay;
  int alphax_1, alphay_1;
  int axi_1, ayj_1;

  Bm = gmatrix<float>(lenx, leny);
  Bx = gmatrix<float>(lenx, leny);
  By = gmatrix<float>(lenx, leny);

  for (i = lenx; i >= 0; i--) {
    for (j = leny; j >= 0; j--) {
      // match states and insertion states together
      alphax = x_alphabet[i];
      alphay = y_alphabet[j];
      if ((i == lenx) && (j == leny)) {
        Bm[i][j] = tran_match[alphax][4];
        Bx[i][j] = tran_x[0][3];
        By[i][j] = tran_y[0][3];
        // cout << "By[i][j]: " << i << " " << j << " " << By[i][j] << endl;
        if (debug > 1)
          cout << i << " " << j << " " << Bm[i][j] << " " << Bx[i][j] << " "
               << By[i][j] << endl;
        continue;
      }
      if (i == lenx) {
        // By[i][j] = By[i][j+1]+tran_y[0][2]+aa_loop[y_alphabet[j+1]][ay[j+1]];
        // Bm[i][j] =
        // By[i][j+1]+tran_match[alphax][3]+aa_loop[y_alphabet[j+1]][ay[j+1]];
        By[i][j] = By[i][j + 1] + tran_y[0][2] + score_bg_y[j + 1];
        Bm[i][j] = By[i][j + 1] + tran_match[alphax][3] + score_bg_y[j + 1];
        // cout << "By[i][j]: " << i << " " << j << " " << By[i][j] << endl;
        // cout << "Bm[i][j]: " << i << " " << j << " " << Bm[i][j] << endl;
        Bx[i][j] = LOG_ZERO;
        if (debug > 1)
          cout << i << " " << j << " " << Bm[i][j] << " " << Bx[i][j] << " "
               << By[i][j] << endl;
        continue;
      }
      if (j == leny) {
        // Bx[i][j] =
        // Bx[i+1][j]+tran_x[alphax][2]+aa_loop[x_alphabet[i+1]][ax[i+1]];
        // Bm[i][j] =
        // Bx[i+1][j]+tran_match[alphax][2]+aa_loop[x_alphabet[i+1]][ax[i+1]];
        Bx[i][j] = Bx[i + 1][j] + tran_x[0][2] + score_bg_x[i + 1];
        Bm[i][j] = Bx[i + 1][j] + tran_match[alphax][2] + score_bg_x[i + 1];
        By[i][j] = LOG_ZERO;
        // cout << "Bx[i][j]: " << i << " " << j << " " << Bx[i][j] << endl;
        // cout << "Bm[i][j]: " << i << " " << j << " " << Bm[i][j] << endl;
        if (debug > 1)
          cout << i << " " << j << " " << Bm[i][j] << " " << Bx[i][j] << " "
               << By[i][j] << endl;
        continue;
      }
      alphax_1 = x_alphabet[i + 1];
      alphay_1 = y_alphabet[j + 1];
      axi_1 = ax[i + 1];
      ayj_1 = ay[j + 1];

      // Bm[i][j] =
      // LOG_ADD(Bm[i+1][j+1]+tran_match[alphax][1]+aa_pair[alphax_1][axi_1][ayj_1]+ss_pair[alphax_1][alphay_1],
      // Bx[i+1][j]+tran_match[alphax][2]+aa_loop[alphax_1][axi_1],
      // By[i][j+1]+tran_match[alphax][3]+aa_loop[alphay_1][ayj_1]); Bx[i][j] =
      // LOG_ADD(Bm[i+1][j+1]+tran_x[alphax][1]+aa_pair[alphax_1][axi_1][ayj_1]+ss_pair[alphax_1][alphay_1],
      // Bx[i+1][j]+tran_x[alphax][2]+aa_loop[alphax_1][axi_1]); By[i][j] =
      // LOG_ADD(Bm[i+1][j+1]+tran_y[alphay][1]+aa_pair[alphax_1][axi_1][ayj_1]+ss_pair[alphax_1][alphay_1],
      // By[i][j+1]+tran_y[alphay][2]+aa_loop[alphay_1][ayj_1]);
      Bm[i][j] = LOG_ADD(
          Bm[i + 1][j + 1] + tran_match[alphax][1] + score_matrix[i + 1][j + 1],
          Bx[i + 1][j] + tran_match[alphax][2] + score_bg_x[i + 1],
          By[i][j + 1] + tran_match[alphax][3] + score_bg_y[j + 1]);
      Bx[i][j] =
          LOG_ADD(Bm[i + 1][j + 1] + tran_x[0][1] + score_matrix[i + 1][j + 1],
                  Bx[i + 1][j] + tran_x[0][2] + score_bg_x[i + 1]);
      By[i][j] =
          LOG_ADD(Bm[i + 1][j + 1] + tran_y[0][1] + score_matrix[i + 1][j + 1],
                  By[i][j + 1] + tran_y[0][2] + score_bg_y[j + 1]);
      if (debug > 1)
        cout << i << " " << j << " " << Bm[i][j] << " " << Bx[i][j] << " "
             << By[i][j] << endl;
      // cout << "Bx[i][j]: " << i << " " << j << " " << Bx[i][j] << endl;
      // cout << "By[i][j]: " << i << " " << j << " " << By[i][j] << endl;
      // cout << "Bm[i][j]: " << i << " " << j << " " << Bm[i][j] << endl;
    }
  }

  // cout << tran_begin[1] << " " << tran_match[x_alphabet[1]][3] << " " <<
  // tran_y[y_alphabet[2]][3] << " " << aa_pair[x_alphabet[1]][ax[1]][ay[1]] << "
  // " << ss_pair[x_alphabet[1]][y_alphabet[1]] << " " <<
  // aa_loop[y_alphabet[2]][ay[2]] << endl;

  // Bm[0][0] =
  // Bm[1][1]+tran_begin[1]+aa_pair[x_alphabet[1]][ax[1]][ay[1]]+ss_pair[x_alphabet[1]][y_alphabet[1]];
  // Bx[0][0] = Bx[1][0]+tran_begin[2]+aa_loop[x_alphabet[1]][ax[1]];
  // By[0][0] = By[0][1]+tran_begin[3]+aa_loop[y_alphabet[1]][ay[1]];

  // BE =
  // LOG_ADD(Bm[1][1]+tran_begin[1]+aa_pair[x_alphabet[1]][ax[1]][ay[1]]+ss_pair[x_alphabet[1]][y_alphabet[1]],
  // Bx[1][0]+tran_begin[2]+aa_loop[x_alphabet[1]][ax[1]],
  // By[0][1]+tran_begin[3]+aa_loop[y_alphabet[1]][ay[1]]);
  BE = LOG_ADD(Bm[1][1] + tran_begin[1] + score_matrix[1][1],
               Bx[1][0] + tran_begin[2] + score_bg_x[1],
               By[0][1] + tran_begin[3] + score_bg_y[1]);

  // cout <<
  // Bm[1][1]+tran_begin[1]+aa_pair[x_alphabet[1]][ax[1]][ay[1]]+ss_pair[x_alphabet[1]][y_alphabet[1]]
  // << " " << Bx[1][0]+tran_begin[2]+aa_loop[x_alphabet[1]][ax[1]] << " " <<
  // By[0][1]+tran_begin[3]+aa_loop[y_alphabet[1]][ay[1]] << endl;

  if (debug > 1) fprintf(stdout, "BE: %f\n", BE);

  // ScoreType full_prob;
  if (debug > 1)
    cout << "Obtain posterior probabilities (>0.01): " << x->aname[0] << "\t"
         << y->aname[0] << endl;
  probMat = gmatrix<float>(lenx, leny);
  for (i = 1; i <= lenx; i++) {
    for (j = 1; j <= leny; j++) {
      probMat[i][j] = exp(Fm[i][j] + Bm[i][j] - FE);
      if (debug > 1)
        cout << "probMat: " << i << " " << j << " " << setprecision(7)
             << probMat[i][j] << endl;
    }
  }

  // fprintf(stdout, "================\n");
}

void hmm_psipred::forward_no_end_penalty() {
  int i, j, k, m;
  int *ax = x->alignment[1];
  int *ay = y->alignment[1];

  int alphax, alphay;
  int alphax_1, alphay_1;
  int axi, ayj;

  Fm = gmatrix<float>(lenx, leny);
  Fx = gmatrix<float>(lenx, leny);
  Fy = gmatrix<float>(lenx, leny);

  if (debug > 1) cout << "Forward algorithm:" << endl;

  cout << setprecision(10);

  for (i = 0; i <= lenx; i++) {
    for (j = 0; j <= leny; j++) {
      if ((i == 0) && (j == 0)) {
        Fm[i][j] = Fx[i][j] = Fy[i][j] = LOG_ZERO;
        if (debug > 1)
          cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
               << Fy[i][j] << endl;
        continue;
      }

      if ((i == 0) && (j == 1)) {
        // cout << tran_begin[3]<<" "<<y_alphabet[j]<<" "<< ay[j]<<" " <<
        // aa_loop[y_alphabet[j]][ay[j]]<<endl; Fy[i][j] = tran_begin[3] +
        // aa_loop[y_alphabet[j]][ay[j]];
        Fy[i][j] = tran_begin[3] + score_bg_y[1];
        // cout << i << " " << j << " " << Fy[i][j] << endl;
        Fx[i][j] = Fm[i][j] = LOG_ZERO;
        if (debug > 1)
          cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
               << Fy[i][j] << endl;
        continue;
      }
      if ((i == 0) && (j > 1)) {
        // cout << Fy[i][j-1] << " " << tran_y[y_alphabet[j-1]][2] << " " <<
        // aa_loop[y_alphabet[j]][ay[j]] << endl; Fy[i][j] = Fy[i][j-1] +
        // tran_y[y_alphabet[j-1]][2] + aa_loop[y_alphabet[j]][ay[j]];
        ////Fy[i][j] = Fy[i][j-1] + tran_y[0][2] + score_bg_y[j];
        Fy[i][j] =
            Fy[i][j - 1] + tran_y[0][2] * weight_end_penalty + score_bg_y[j];
        Fx[i][j] = Fm[i][j] = LOG_ZERO;
        // cout << i << " " << j << " " << Fy[i][j] << endl;
        if (debug > 1)
          cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
               << Fy[i][j] << endl;
        continue;
      }
      if ((i == 1) && (j == 0)) {
        // cout << tran_begin[2] << " " << aa_loop[x_alphabet[i]][ax[i]] <<
        // endl; Fx[i][j] = tran_begin[2] + aa_loop[x_alphabet[i]][ax[i]];
        Fx[i][j] = tran_begin[2] + score_bg_x[1];
        Fy[i][j] = Fm[i][j] = LOG_ZERO;
        // cout << i << " " << j << " " << Fx[i][j] << endl;
        if (debug > 1)
          cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
               << Fy[i][j] << endl;
        continue;
      }
      if ((i > 1) && (j == 0)) {
        // Fx[i][j] = Fx[i-1][j] + tran_x[x_alphabet[i-1]][2] +
        // aa_loop[x_alphabet[i]][ax[i]];
        ////Fx[i][j] = Fx[i-1][j] + tran_x[0][2] + score_bg_x[i];
        Fx[i][j] =
            Fx[i - 1][j] + tran_x[0][2] * weight_end_penalty + score_bg_x[i];
        // cout << i << " " << j << " " << Fx[i][j] << endl;
        Fy[i][j] = Fm[i][j] = LOG_ZERO;
        if (debug > 1)
          cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
               << Fy[i][j] << endl;
        continue;
      }
      alphax = x_alphabet[i];
      alphay = y_alphabet[j];
      axi = ax[i];
      ayj = ay[j];
      if ((i == 1) && (j == 1)) {
        // Fm[i][j] = tran_begin[1] + aa_pair[alphax][axi][ayj] +
        // ss_pair[alphax][alphay];
        Fm[i][j] = tran_begin[1] + score_matrix[1][1];
        // cout << i << " " << j << " " << Fm[i][j] << endl;
        Fx[i][j] = Fy[i][j] = LOG_ZERO;
        if (debug > 1)
          cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
               << Fy[i][j] << endl;
        continue;
      }

      alphax_1 = x_alphabet[i - 1];
      alphay_1 = y_alphabet[j - 1];
      if (i == 1) alphax_1 = 1;
      if (j == 1) alphay_1 = 3;
      // Fm[i][j] = LOG_ADD(tran_match[alphax_1][1]+Fm[i-1][j-1],
      // tran_x[0][1]+Fx[i-1][j-1], tran_y[0][1]+Fy[i-1][j-1]) +
      // p->aa_pair[alphax][axi][ayj] + p->ss_pair[alphax][alphay]; Fx[i][j] =
      // LOG_ADD(tran_match[alphax_1][2]+Fm[i-1][j], tran_x[0][2]+Fx[i-1][j]) +
      // p->aa_loop[alphax][axi]; Fy[i][j] =
      // LOG_ADD(tran_match[alphax][3]+Fm[i][j-1],
      // tran_y[alphay_1][2]+Fy[i][j-1]) + p->aa_loop[alphay][ayj];
      if (i == 1) {
        Fm[i][j] =
            LOG_ADD(tran_match[alphax_1][1] + Fm[i - 1][j - 1],
                    tran_x[0][1] * weight_end_penalty + Fx[i - 1][j - 1],
                    tran_y[0][1] * weight_end_penalty + Fy[i - 1][j - 1]) +
            score_matrix[i][j];
      } else if (j == 1) {
        Fm[i][j] =
            LOG_ADD(tran_match[alphax_1][1] + Fm[i - 1][j - 1],
                    tran_x[0][1] * weight_end_penalty + Fx[i - 1][j - 1],
                    tran_y[0][1] * weight_end_penalty + Fy[i - 1][j - 1]) +
            score_matrix[i][j];
      } else {
        Fm[i][j] = LOG_ADD(tran_match[alphax_1][1] + Fm[i - 1][j - 1],
                           tran_x[0][1] + Fx[i - 1][j - 1],
                           tran_y[0][1] + Fy[i - 1][j - 1]) +
                   score_matrix[i][j];
      }

      if (j != leny) {
        Fx[i][j] = LOG_ADD(tran_match[alphax_1][2] + Fm[i - 1][j],
                           tran_x[0][2] + Fx[i - 1][j]) +
                   score_bg_x[i];
      } else {
        Fx[i][j] =
            LOG_ADD(tran_match[alphax_1][2] * weight_end_penalty + Fm[i - 1][j],
                    tran_x[0][2] * weight_end_penalty + Fx[i - 1][j]) +
            score_bg_x[i];
      }
      if (i != lenx) {
        Fy[i][j] = LOG_ADD(tran_match[alphax][3] + Fm[i][j - 1],
                           tran_y[0][2] + Fy[i][j - 1]) +
                   score_bg_y[j];
      } else {
        Fy[i][j] =
            LOG_ADD(tran_match[alphax][3] * weight_end_penalty + Fm[i][j - 1],
                    tran_y[0][2] * weight_end_penalty + Fy[i][j - 1]) +
            score_bg_y[j];
      }

      if (debug > 1)
        cout << i << " " << j << " " << Fm[i][j] << " " << Fx[i][j] << " "
             << Fy[i][j] << endl;
      // cout << "Fy[i][j]: " << i << " " << j << " " << tran_match[alphax_1][3]
      // << " " << Fm[i][j-1] << " " << tran_y[alphay_1][2] << " " << Fy[i][j-1]
      // << " " << p->aa_loop[alphay][ayj] << endl;
    }
  }

  // the ending state
  // FE = LOG_ADD(Fm[lenx][leny]+tran_match[x_alphabet[lenx]][4],
  // Fx[lenx][leny]+tran_x[x_alphabet[lenx]][3],
  // Fy[lenx][leny]+tran_y[y_alphabet[leny]][3]);
  FE = LOG_ADD(Fm[lenx][leny] + tran_match[x_alphabet[lenx]][4],
               Fx[lenx][leny] + tran_x[0][3], Fy[lenx][leny] + tran_y[0][3]);
  // cout << Fm[lenx][leny]+tran_match[x_alphabet[lenx]][4] << " " <<
  // Fx[lenx][leny]+tran_x[x_alphabet[lenx]][3] << " " <<
  // Fy[lenx][leny]+tran_y[y_alphabet[leny]][3] << endl;

  // cout << "tran_y[y_alphabet[leny]][3]: " << tran_y[y_alphabet[leny]][3] <<
  // endl;

  if (debug > -1) fprintf(stdout, "FE: %12.5f\t", FE);
  fflush(stdout);  // cout << "FE: " << FE << endl;
}

void hmm_psipred::backward_no_end_penalty() {
  int i, j, k, m;
  int *ax = x->alignment[1];
  int *ay = y->alignment[1];

  int alphax, alphay;
  int alphax_1, alphay_1;
  int axi_1, ayj_1;

  Bm = gmatrix<float>(lenx, leny);
  Bx = gmatrix<float>(lenx, leny);
  By = gmatrix<float>(lenx, leny);

  for (i = lenx; i >= 0; i--) {
    for (j = leny; j >= 0; j--) {
      // match states and insertion states together
      alphax = x_alphabet[i];
      alphay = y_alphabet[j];
      if ((i == lenx) && (j == leny)) {
        Bm[i][j] = tran_match[alphax][4];
        Bx[i][j] = tran_x[0][3];
        By[i][j] = tran_y[0][3];
        // cout << "By[i][j]: " << i << " " << j << " " << By[i][j] << endl;
        if (debug > 1)
          cout << i << " " << j << " " << Bm[i][j] << " " << Bx[i][j] << " "
               << By[i][j] << endl;
        continue;
      }
      if (i == lenx) {
        // By[i][j] = By[i][j+1]+tran_y[0][2]+aa_loop[y_alphabet[j+1]][ay[j+1]];
        // Bm[i][j] =
        // By[i][j+1]+tran_match[alphax][3]+aa_loop[y_alphabet[j+1]][ay[j+1]];
        ////By[i][j] = By[i][j+1]+tran_y[0][2]+score_bg_y[j+1];
        ////Bm[i][j] = By[i][j+1]+tran_match[alphax][3]+score_bg_y[j+1];
        By[i][j] = By[i][j + 1] + tran_y[0][2] * weight_end_penalty +
                   score_bg_y[j + 1];
        Bm[i][j] = By[i][j + 1] + tran_match[alphax][3] * weight_end_penalty +
                   score_bg_y[j + 1];
        // cout << "By[i][j]: " << i << " " << j << " " << By[i][j] << endl;
        // cout << "Bm[i][j]: " << i << " " << j << " " << Bm[i][j] << endl;
        Bx[i][j] = LOG_ZERO;
        if (debug > 1)
          cout << i << " " << j << " " << Bm[i][j] << " " << Bx[i][j] << " "
               << By[i][j] << endl;
        continue;
      }
      if (j == leny) {
        // Bx[i][j] =
        // Bx[i+1][j]+tran_x[alphax][2]+aa_loop[x_alphabet[i+1]][ax[i+1]];
        // Bm[i][j] =
        // Bx[i+1][j]+tran_match[alphax][2]+aa_loop[x_alphabet[i+1]][ax[i+1]];
        ////Bx[i][j] = Bx[i+1][j]+tran_x[0][2]+score_bg_x[i+1];
        ////Bm[i][j] = Bx[i+1][j]+tran_match[alphax][2]+score_bg_x[i+1];
        Bx[i][j] = Bx[i + 1][j] + tran_x[0][2] * weight_end_penalty +
                   score_bg_x[i + 1];
        Bm[i][j] = Bx[i + 1][j] + tran_match[alphax][2] * weight_end_penalty +
                   score_bg_x[i + 1];
        By[i][j] = LOG_ZERO;
        // cout << "Bx[i][j]: " << i << " " << j << " " << Bx[i][j] << endl;
        // cout << "Bm[i][j]: " << i << " " << j << " " << Bm[i][j] << endl;
        if (debug > 1)
          cout << i << " " << j << " " << Bm[i][j] << " " << Bx[i][j] << " "
               << By[i][j] << endl;
        continue;
      }
      alphax_1 = x_alphabet[i + 1];
      alphay_1 = y_alphabet[j + 1];
      axi_1 = ax[i + 1];
      ayj_1 = ay[j + 1];

      // Bm[i][j] =
      // LOG_ADD(Bm[i+1][j+1]+tran_match[alphax][1]+aa_pair[alphax_1][axi_1][ayj_1]+ss_pair[alphax_1][alphay_1],
      // Bx[i+1][j]+tran_match[alphax][2]+aa_loop[alphax_1][axi_1],
      // By[i][j+1]+tran_match[alphax][3]+aa_loop[alphay_1][ayj_1]); Bx[i][j] =
      // LOG_ADD(Bm[i+1][j+1]+tran_x[alphax][1]+aa_pair[alphax_1][axi_1][ayj_1]+ss_pair[alphax_1][alphay_1],
      // Bx[i+1][j]+tran_x[alphax][2]+aa_loop[alphax_1][axi_1]); By[i][j] =
      // LOG_ADD(Bm[i+1][j+1]+tran_y[alphay][1]+aa_pair[alphax_1][axi_1][ayj_1]+ss_pair[alphax_1][alphay_1],
      // By[i][j+1]+tran_y[alphay][2]+aa_loop[alphay_1][ayj_1]);
      Bm[i][j] = LOG_ADD(
          Bm[i + 1][j + 1] + tran_match[alphax][1] + score_matrix[i + 1][j + 1],
          Bx[i + 1][j] + tran_match[alphax][2] + score_bg_x[i + 1],
          By[i][j + 1] + tran_match[alphax][3] + score_bg_y[j + 1]);
      if (j != 0) {
        Bx[i][j] = LOG_ADD(
            Bm[i + 1][j + 1] + tran_x[0][1] + score_matrix[i + 1][j + 1],
            Bx[i + 1][j] + tran_x[0][2] + score_bg_x[i + 1]);
      } else {
        Bx[i][j] =
            LOG_ADD(Bm[i + 1][j + 1] + tran_x[0][1] * weight_end_penalty +
                        score_matrix[i + 1][j + 1],
                    Bx[i + 1][j] + tran_x[0][2] * weight_end_penalty +
                        score_bg_x[i + 1]);
      }
      if (i != 0) {
        By[i][j] = LOG_ADD(
            Bm[i + 1][j + 1] + tran_y[0][1] + score_matrix[i + 1][j + 1],
            By[i][j + 1] + tran_y[0][2] + score_bg_y[j + 1]);
      } else {
        By[i][j] =
            LOG_ADD(Bm[i + 1][j + 1] + tran_y[0][1] * weight_end_penalty +
                        score_matrix[i + 1][j + 1],
                    By[i][j + 1] + tran_y[0][2] * weight_end_penalty +
                        score_bg_y[j + 1]);
      }
      if (debug > 1)
        cout << i << " " << j << " " << Bm[i][j] << " " << Bx[i][j] << " "
             << By[i][j] << endl;
      // cout << "Bx[i][j]: " << i << " " << j << " " << Bx[i][j] << endl;
      // cout << "By[i][j]: " << i << " " << j << " " << By[i][j] << endl;
      // cout << "Bm[i][j]: " << i << " " << j << " " << Bm[i][j] << endl;
    }
  }

  // cout << tran_begin[1] << " " << tran_match[x_alphabet[1]][3] << " " <<
  // tran_y[y_alphabet[2]][3] << " " << aa_pair[x_alphabet[1]][ax[1]][ay[1]] << "
  // " << ss_pair[x_alphabet[1]][y_alphabet[1]] << " " <<
  // aa_loop[y_alphabet[2]][ay[2]] << endl;

  // Bm[0][0] =
  // Bm[1][1]+tran_begin[1]+aa_pair[x_alphabet[1]][ax[1]][ay[1]]+ss_pair[x_alphabet[1]][y_alphabet[1]];
  // Bx[0][0] = Bx[1][0]+tran_begin[2]+aa_loop[x_alphabet[1]][ax[1]];
  // By[0][0] = By[0][1]+tran_begin[3]+aa_loop[y_alphabet[1]][ay[1]];

  // BE =
  // LOG_ADD(Bm[1][1]+tran_begin[1]+aa_pair[x_alphabet[1]][ax[1]][ay[1]]+ss_pair[x_alphabet[1]][y_alphabet[1]],
  // Bx[1][0]+tran_begin[2]+aa_loop[x_alphabet[1]][ax[1]],
  // By[0][1]+tran_begin[3]+aa_loop[y_alphabet[1]][ay[1]]);
  BE = LOG_ADD(Bm[1][1] + tran_begin[1] + score_matrix[1][1],
               Bx[1][0] + tran_begin[2] + score_bg_x[1],
               By[0][1] + tran_begin[3] + score_bg_y[1]);

  // cout <<
  // Bm[1][1]+tran_begin[1]+aa_pair[x_alphabet[1]][ax[1]][ay[1]]+ss_pair[x_alphabet[1]][y_alphabet[1]]
  // << " " << Bx[1][0]+tran_begin[2]+aa_loop[x_alphabet[1]][ax[1]] << " " <<
  // By[0][1]+tran_begin[3]+aa_loop[y_alphabet[1]][ay[1]] << endl;

  if (debug > -1) fprintf(stdout, "BE: %f\n", BE);

  // ScoreType full_prob;
  if (debug > 1)
    cout << "Obtain posterior probabilities (>0.01): " << x->aname[0] << "\t"
         << y->aname[0] << endl;
  probMat = gmatrix<float>(lenx, leny);
  for (i = 1; i <= lenx; i++) {
    for (j = 1; j <= leny; j++) {
      probMat[i][j] = exp(Fm[i][j] + Bm[i][j] - FE);
      // double tmp = exp(Fm[i][j]+Bm[i][j]-FE) + exp(Fx[i][j]+Bx[i][j]-FE) +
      // exp(Fy[i][j]+By[i][j]-FE); cout << "all three - i: " << i << " j: " << j
      // << " "  << tmp << endl;
      if (debug > 1)
        cout << "probMat: " << i << " " << j << " " << setprecision(7)
             << probMat[i][j] << endl;
    }
  }

  // fprintf(stdout, "================\n");
}

int *hmm_psipred::posteriorAlignment() {
  int i, j, k;
  int *path = new int[lenx + leny + 1];

  int **scoreMat = imatrix(lenx, leny);

  for (i = 1; i <= lenx; i++) {
    for (j = 1; j <= leny; j++) {
      scoreMat[i][j] = (int)(probMat[i][j] * 1000);
    }
  }

  MM galign;
  galign.setM(lenx);
  galign.setN(leny);
  galign.set_g(0);  // gap open penalty
  galign.set_h(0);  // gap extension penalty

  cout << "Compute pairwise consistency alignment:" << endl;
  galign.dp(scoreMat);
  for (i = 1; i < galign.print_ptr; i++) {
    path[i] = galign.displ[i];
    // cout << "i: " << i << " " << path[i] << endl;
  }
  // len = galign.print_ptr - 1;
  // cout << "len: " << len << endl;
  // cout << "Pairwise consistency alignment ends here" << endl;

  free_imatrix(scoreMat, lenx, leny);

  char *a, *b;
  a = new char[3000];
  b = new char[3000];

  int i1, i2;
  i1 = i2 = 0;
  int x1 = 0;
  for (i = 1; i < galign.print_ptr; i++) {
    if (path[i] < 0) {
      for (k = 1; k <= 0 - path[i]; k++) {
        a[x1] = x->aseq[0][i1];
        b[x1] = '-';
        i1++;
        x1++;
      }
      continue;
    }
    if (path[i] > 0) {
      for (k = 1; k <= path[i]; k++) {
        a[x1] = '-';
        b[x1] = y->aseq[0][i2];
        i2++;
        x1++;
      }
      continue;
    }
    if (path[i] == 0) {
      a[x1] = x->aseq[0][i1];
      b[x1] = y->aseq[0][i2];
      i1++;
      i2++;
      x1++;
    }
  }
  a[x1] = '\0';
  b[x1] = '\0';

  cout << left << setw(40) << x->aname[0] << a << endl;
  cout << left << setw(40) << y->aname[0] << b << endl;

  delete[] a;
  delete[] b;

  return path;
}
