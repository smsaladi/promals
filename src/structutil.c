#include "btree_template.h"
#include "constraint.h"
#include "multiple.h"
#include "param.h"
#include "progressiveAlignHMM.h"
#include "sequences.h"
#include "subalign.h"

static int debug_here = 11;
void get_added_gap_number_1(char *seq1, char *seq2, int *added_gap_number_1,
                            int dim);

hmm_psipred_parameters *params1;

#define MIN_HERE(x, y) (((x) < (y)) ? (x) : (y))
#define MAX_HERE(x, y) (((x) > (y)) ? (x) : (y))

/*
// combine structure alignment only for dali
void combine_structure_alignments(btree<tnode> &alltree, vector<seq_str_aln *>
&ssaln) {

        int i, j, k, l;

        float **structmat, **tmpMat;
        sparseMatrix ***smat = alltree.smat;

        cout << "sequence weight: " << sequence_weight << endl;

        for(i=0;i<(int)alltree.preAligned.size();i++) {
                for(j=i+1;j<(int)alltree.preAligned.size();j++) {
                        //cout << "+++++++++++" << i << " " << j << endl;
                        structmat = map_structure(ssaln[i], ssaln[j]);
                        //cout << "This place" << endl;
                        if(!structmat) continue;
                        //cout << "This place" << endl;
                        tmpMat = smat[i][j]->sparseCrs2Regular();
                        for(k=1;k<=alltree.preAligned[i]->aln->alilen;k++) {
                                for(l=1;l<=alltree.preAligned[j]->aln->alilen;l++)
{
                                        //tmpMat[k][l] +=
structmat[k][l]*struct_weight; tmpMat[k][l] = tmpMat[k][l]*sequence_weight +
structmat[k][l]*struct_weight;
                                }
                        }
                        //cout << "This place" << endl;
                        delete smat[i][j];
                        delete smat[j][i];
                        smat[i][j] = new
sparseMatrix(tmpMat,alltree.preAligned[i]->aln->alilen,
alltree.preAligned[j]->aln->alilen); smat[j][i] = smat[i][j]->transpose();
                        free_gmatrix<float>(tmpMat,
alltree.preAligned[i]->aln->alilen, alltree.preAligned[j]->aln->alilen);
                        free_gmatrix<float>(structmat,
alltree.preAligned[i]->aln->alilen, alltree.preAligned[j]->aln->alilen);
                        //cout << "This place" << endl;
                }
        }
}
*/

// combine structure alignment for dali, fast and tmalign
void combine_structure_alignments3(btree<tnode> &alltree,
                                   vector<seq_str_aln *> &ssaln) {
  int i, j, k, l;

  print_time_diff("consistency_scoring");
  print_section_info("Below combine structural constraints, old");
  cout << "realign_psiblast: " << realign_psiblast << endl << endl;
  float **structmat_dali = NULL, **structmat_fast = NULL,
        **structmat_tmalign = NULL, **tmpMat, **structmat;
  sparseMatrix ***smat = alltree.smat;

  // cout << "sequence weight: " << sequence_weight << endl;

  for (i = 0; i < (int)alltree.preAligned.size(); i++) {
    for (j = i + 1; j < (int)alltree.preAligned.size(); j++) {
      // cout << "+++++++++++" << i << " " << j << endl;
      structmat_dali = NULL;
      structmat_fast = NULL;
      structmat_tmalign = NULL;
      if (realign_psiblast) {
        if (use_dali)
          structmat_dali = map_structure_promals(
              alltree.preAligned[i]->aux_align,
              alltree.preAligned[j]->aux_align, ssaln[i], ssaln[j], "dali");
        if (use_fast)
          structmat_fast = map_structure_promals(
              alltree.preAligned[i]->aux_align,
              alltree.preAligned[j]->aux_align, ssaln[i], ssaln[j], "fast");
        cout << "Here fast: " << endl;
        if (use_tmalign)
          structmat_tmalign = map_structure_promals(
              alltree.preAligned[i]->aux_align,
              alltree.preAligned[j]->aux_align, ssaln[i], ssaln[j], "tmalign");
      } else {
        if (use_dali)
          structmat_dali = map_structure_nopromals(ssaln[i], ssaln[j], "dali");
        if (use_fast)
          structmat_fast = map_structure_nopromals(ssaln[i], ssaln[j], "fast");
        if (use_tmalign)
          structmat_tmalign =
              map_structure_nopromals(ssaln[i], ssaln[j], "tmalign");
      }
      if ((!structmat_dali) && (!structmat_fast) && (!structmat_tmalign))
        continue;
      structmat = gmatrix<float>(alltree.preAligned[i]->aln->alilen,
                                 alltree.preAligned[j]->aln->alilen);
      for (k = 1; k <= alltree.preAligned[i]->aln->alilen; k++) {
        for (l = 1; l <= alltree.preAligned[j]->aln->alilen; l++) {
          structmat[k][l] = 0;
        }
      }
      if (use_dali)
        combine_structmat(structmat, structmat_dali,
                          alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen, 1.0);
      if (use_fast)
        combine_structmat(structmat, structmat_fast,
                          alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen, 1.0);
      // cout << "Here fast2: " << endl;
      if (use_tmalign)
        combine_structmat(structmat, structmat_tmalign,
                          alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen, 1.0);
      // debugging
      /*cout << alltree.preAligned[i]->aln->aname[0] << " " <<
      alltree.preAligned[j]->aln->aname[0] << endl;
      for(k=1;k<=alltree.preAligned[i]->aln->alilen;k++) {
              for(l=1;l<=alltree.preAligned[j]->aln->alilen;l++) {
                      if(structmat_dali[k][l]) {
                              cout << k << " " << l << " " <<
      alltree.preAligned[i]->aln->aseq[0][k-1] << " " <<
      alltree.preAligned[j]->aln->aseq[0][l-1] << " " << structmat_dali[k][l] <<
      endl;
                      }
              }
      }*/

      // cout << "This place" << endl;
      tmpMat = smat[i][j]->sparseCrs2Regular();
      for (k = 1; k <= alltree.preAligned[i]->aln->alilen; k++) {
        for (l = 1; l <= alltree.preAligned[j]->aln->alilen; l++) {
          // tmpMat[k][l] += structmat[k][l]*struct_weight;
          tmpMat[k][l] =
              tmpMat[k][l] * sequence_weight + structmat[k][l] * struct_weight;
        }
      }
      // cout << "This place" << endl;
      delete smat[i][j];
      delete smat[j][i];
      smat[i][j] = new sparseMatrix(tmpMat, alltree.preAligned[i]->aln->alilen,
                                    alltree.preAligned[j]->aln->alilen);
      smat[j][i] = smat[i][j]->transpose();
      free_gmatrix<float>(tmpMat, alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen);
      free_gmatrix<float>(structmat, alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen);
      free_gmatrix<float>(structmat_dali, alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen);
      free_gmatrix<float>(structmat_fast, alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen);
      free_gmatrix<float>(structmat_tmalign, alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen);
      // cout << "This place" << endl;
    }
  }
}

// combine structure alignment for dali, fast and tmalign
// use updated database
void combine_structure_alignments3_updated(btree<tnode> &alltree,
                                           vector<seq_str_aln *> &ssaln) {
  int i, j, k, l;

  print_section_info("Below combine structural constraints, updated");
  cout << "realign_psiblast: " << realign_psiblast << endl << endl;
  float **structmat_dali = NULL, **structmat_fast = NULL,
        **structmat_tmalign = NULL, **tmpMat, **structmat;
  sparseMatrix ***smat = alltree.smat;

  // cout << "sequence weight: " << sequence_weight << endl;

  for (i = 0; i < (int)alltree.preAligned.size(); i++) {
    for (j = i + 1; j < (int)alltree.preAligned.size(); j++) {
      // cout << "+++++++++++" << i << " " << j << endl;
      structmat_dali = NULL;
      structmat_fast = NULL;
      structmat_tmalign = NULL;
      if (realign_psiblast) {
        structmat_dali = map_structure_promals_updated(
            alltree.preAligned[i]->aux_align, alltree.preAligned[j]->aux_align,
            ssaln[i], ssaln[j], use_dali, use_fast, use_tmalign);
        // if(use_dali) structmat_dali = map_structure(ssaln[i], ssaln[j]);
        /*
        if(use_dali) structmat_dali =
        map_structure_promals_updated(alltree.preAligned[i]->aux_align,
        alltree.preAligned[j]->aux_align, ssaln[i], ssaln[j], "dali");
        //if(use_fast) structmat_fast = map_structure_fast(ssaln[i], ssaln[j]);
        if(use_fast) structmat_fast =
        map_structure_promals_updated(alltree.preAligned[i]->aux_align,
        alltree.preAligned[j]->aux_align, ssaln[i], ssaln[j], "fast");
        //cout << "Here fast updated: " << endl;
        //if(use_tmalign) structmat_tmalign = map_structure_tmalign(ssaln[i],
        ssaln[j]); if(use_tmalign) structmat_tmalign =
        map_structure_promals_updated(alltree.preAligned[i]->aux_align,
        alltree.preAligned[j]->aux_align, ssaln[i], ssaln[j], "tmalign");
        */
      } else {
        if (use_dali)
          structmat_dali =
              map_structure_nopromals_updated(ssaln[i], ssaln[j], "dali");
        if (use_fast)
          structmat_fast =
              map_structure_nopromals_updated(ssaln[i], ssaln[j], "fast");
        if (use_tmalign)
          structmat_tmalign =
              map_structure_nopromals_updated(ssaln[i], ssaln[j], "tmalign");
      }
      if ((!structmat_dali) && (!structmat_fast) && (!structmat_tmalign))
        continue;
      structmat = gmatrix<float>(alltree.preAligned[i]->aln->alilen,
                                 alltree.preAligned[j]->aln->alilen);
      for (k = 1; k <= alltree.preAligned[i]->aln->alilen; k++) {
        for (l = 1; l <= alltree.preAligned[j]->aln->alilen; l++) {
          structmat[k][l] = 0;
        }
      }
      if (use_dali)
        combine_structmat(structmat, structmat_dali,
                          alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen, 1.0);
      if (use_fast)
        combine_structmat(structmat, structmat_fast,
                          alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen, 1.0);
      // cout << "Here fast2: " << endl;
      if (use_tmalign)
        combine_structmat(structmat, structmat_tmalign,
                          alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen, 1.0);
      // debugging
      /*cout << alltree.preAligned[i]->aln->aname[0] << " " <<
      alltree.preAligned[j]->aln->aname[0] << endl;
      for(k=1;k<=alltree.preAligned[i]->aln->alilen;k++) {
              for(l=1;l<=alltree.preAligned[j]->aln->alilen;l++) {
                      if(structmat_dali[k][l]) {
                              cout << k << " " << l << " " <<
      alltree.preAligned[i]->aln->aseq[0][k-1] << " " <<
      alltree.preAligned[j]->aln->aseq[0][l-1] << " " << structmat_dali[k][l] <<
      endl;
                      }
              }
      }*/

      // cout << "This place" << endl;
      tmpMat = smat[i][j]->sparseCrs2Regular();
      for (k = 1; k <= alltree.preAligned[i]->aln->alilen; k++) {
        for (l = 1; l <= alltree.preAligned[j]->aln->alilen; l++) {
          // tmpMat[k][l] += structmat[k][l]*struct_weight;
          tmpMat[k][l] =
              tmpMat[k][l] * sequence_weight + structmat[k][l] * struct_weight;
        }
      }
      // cout << "This place" << endl;
      delete smat[i][j];
      delete smat[j][i];
      smat[i][j] = new sparseMatrix(tmpMat, alltree.preAligned[i]->aln->alilen,
                                    alltree.preAligned[j]->aln->alilen);
      smat[j][i] = smat[i][j]->transpose();
      free_gmatrix<float>(tmpMat, alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen);
      free_gmatrix<float>(structmat, alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen);
      free_gmatrix<float>(structmat_dali, alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen);
      free_gmatrix<float>(structmat_fast, alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen);
      free_gmatrix<float>(structmat_tmalign, alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen);
      // cout << "This place" << endl;
    }
  }
  print_time_diff("calculate and combine structural_aligments");
}

//
void combine_structmat(float **mat1, float **mat2, int d1, int d2,
                       float weight) {
  int i, j;
  if (!mat2) return;

  for (i = 1; i <= d1; i++) {
    for (j = 1; j <= d2; j++) {
      if (mat2[i][j]) {
        if (mat1[i][j]) {
          mat1[i][j] += weight;
        } else {
          mat1[i][j] = mat2[i][j];
        }
      }
    }
  }
  return;
}

// combine structure alignment for dali, fast and tmalign
void combine_constraint(btree<tnode> &alltree, constraint &cons,
                        float constraint_w) {
  int i, j, k, l;
  int ci, cj;  // index for sequences in cons

  // float **structmat_dali=NULL, **structmat_fast=NULL,
  // **structmat_tmalign=NULL, **tmpMat, **structmat;
  float **structmat_constraint = NULL, **tmpMat;
  sparseMatrix ***smat = alltree.smat;

  // cout << "sequence weight: " << sequence_weight << endl;

  for (i = 0; i < (int)alltree.preAligned.size(); i++) {
    // check if index i can be found in the constraint rep_index
    ci = -1;
    for (k = 1; k <= cons.nseqs; k++) {
      if (cons.rep_index[k] == i) ci = k;
    }
    if (ci == -1)
      continue;  // if not found, there is no constraint on prealigned group i
    for (j = i + 1; j < (int)alltree.preAligned.size(); j++) {
      // check if index j can be found in the constraint rep_index
      cj = -1;
      for (k = 1; k <= cons.nseqs; k++) {
        if (cons.rep_index[k] == j) cj = k;
      }
      if (cj == -1)
        continue;  // if not found, there is no constraint on prealigned group i

      // cout << "+++++++++++" << i << " " << j << endl;
      structmat_constraint = NULL;
      structmat_constraint = cons.get_constraint_matrix(ci, cj);
      if (!structmat_constraint) continue;
      // cout << "This place" << endl;
      tmpMat = smat[i][j]->sparseCrs2Regular();
      for (k = 1; k <= alltree.preAligned[i]->aln->alilen; k++) {
        for (l = 1; l <= alltree.preAligned[j]->aln->alilen; l++) {
          // tmpMat[k][l] += structmat[k][l]*struct_weight;
          tmpMat[k][l] =
              tmpMat[k][l] + structmat_constraint[k][l] * constraint_w;
          cout << k << " " << l << " "
               << alltree.preAligned[i]->aln->aseq[0][k - 1] << " "
               << alltree.preAligned[j]->aln->aseq[0][l - 1] << " "
               << tmpMat[k][l] << endl;
        }
      }
      // cout << "This place" << endl;
      delete smat[i][j];
      delete smat[j][i];
      smat[i][j] = new sparseMatrix(tmpMat, alltree.preAligned[i]->aln->alilen,
                                    alltree.preAligned[j]->aln->alilen);
      smat[j][i] = smat[i][j]->transpose();
      free_gmatrix<float>(tmpMat, alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen);
      free_gmatrix<float>(structmat_constraint,
                          alltree.preAligned[i]->aln->alilen,
                          alltree.preAligned[j]->aln->alilen);
      // cout << "This place" << endl;
    }
  }
}

// get aligned residue positions for two dali id names
int **get_salnpos(char *name1, char *name2, int &numalignedp) {
  int i, j, k;
  char dali_dir[100];
  strcpy(dali_dir, "/usr2/db/dalilite_aln_1.69/");

  char result_file[100];
  char index_file[100];
  sprintf(result_file, "%s%s.result", dali_dir, name1);
  sprintf(index_file, "%s%s.result.idx", dali_dir, name1);

  // cout << result_file << endl;
  // cout << index_file << endl;

  int name2_index;
  sscanf(name2, "%d", &name2_index);
  // cout << name2_index << endl;
  int address_daliscore[2];
  FILE *idxfp = fopen(index_file, "r");
  fseek(idxfp, (name2_index - 1) * 8, 0);
  fread(address_daliscore, 4, 2, idxfp);
  // cout << "names: " << name1 << " " << name2 << " " << address_daliscore[0]
  // << " " << address_daliscore[1] << endl;
  fclose(idxfp);
  if (address_daliscore[1] < minDaliZ * 10) return NULL;
  FILE *resultfp = fopen(result_file, "r");
  fseek(resultfp, address_daliscore[0], 0);

  char alnseq1[3000], alnseq2[3000];

  fgets(alnseq1, 3000, resultfp);
  fgets(alnseq1, 3000, resultfp);
  fgets(alnseq1, 3000, resultfp);
  fgets(alnseq2, 3000, resultfp);

  numalignedp = 0;
  for (i = 0; i < strlen(alnseq1); i++) {
    if ((alnseq1[i] >= 65) && (alnseq1[i] <= 90)) {
      numalignedp++;
    }
  }
  // cout << "numalignedp: " << numalignedp << endl;
  // cout << "alnseq1: " << alnseq1 << endl;
  // cout << "alnseq2: " << alnseq2 << endl;
  int **salpos = imatrix(1, numalignedp);
  int count1 = 0;  // residue number in sequence 1
  int count2 = 0;  // residue number in sequence 2
  int countp = 0;  // number of aligned positions
  for (i = 0; i < strlen(alnseq1); i++) {
    if (alnseq1[i] != '-') count1++;
    if (alnseq2[i] != '-') count2++;
    if ((alnseq1[i] >= 65) && (alnseq1[i] <= 90)) {
      salpos[0][countp] = count1;
      salpos[1][countp] = count2;
      // cout << countp << " " <<  salpos[0][countp] << " " << salpos[1][countp]
      // << endl;
      countp++;
    }
  }
  // cout << "HEERElllll" <<endl;
  fclose(resultfp);
  return salpos;
}

// get aligned residue positions for two dali id names
int **get_salnpos_fast(char *name1, char *name2, int &numalignedp) {
  int i, j, k;
  char dali_dir[100], index_dir[100];
  strcpy(dali_dir, "/usr2/db/fast_1.69/");
  // strcpy(dali_dir, "/home/jpei/db/fast_1.69/");
  strcpy(index_dir, "/home/jpei/db/fast_1.69/");

  char result_file[100];
  char index_file[100];
  sprintf(result_file, "%s%s.fast.result", dali_dir, name1);
  sprintf(index_file, "%s%s.fast.result.index", index_dir, name1);
  cout << result_file << endl;
  cout << index_file << endl;

  // cout << result_file << endl;
  // cout << index_file << endl;

  int name2_index;
  sscanf(name2, "%d", &name2_index);
  // cout << name2_index << endl;
  int address_daliscore[2];
  FILE *idxfp = fopen(index_file, "r");
  fseek(idxfp, (name2_index)*4, 0);
  fread(address_daliscore, 4, 1, idxfp);
  cout << name1 << " " << name2 << " " << address_daliscore[0] << endl;
  fclose(idxfp);
  // if(address_daliscore[1]<minDaliZ * 10) return NULL;
  FILE *resultfp = fopen(result_file, "r");
  fseek(resultfp, address_daliscore[0], 0);

  char alnseq1[3000], alnseq2[3000];

  fgets(alnseq1, 3000, resultfp);
  // check names
  char combinednames[100];
  sprintf(combinednames, "%s %s\n", name1, name2);
  if (strcmp(alnseq1, combinednames) != 0) {
    cout << "combined names: " << combinednames << endl;
    cout << "combined names in result file: " << alnseq1 << endl;
    cout << "not the same name" << endl;
    exit(0);
  }
  cout << "combined names: " << combinednames << endl;
  cout << "combined names in result file: " << alnseq1 << endl;
  // fgets( alnseq1, 3000, resultfp);
  fgets(alnseq1, 3000, resultfp);
  fgets(alnseq2, 3000, resultfp);

  numalignedp = 0;
  for (i = 0; i < strlen(alnseq1); i++) {
    if ((alnseq1[i] >= 65) && (alnseq1[i] <= 90)) {
      numalignedp++;
    }
  }
  // cout << "numalignedp: " << numalignedp << endl;
  // cout << alnseq1 << endl;
  // cout << alnseq2 << endl;
  int **salpos = imatrix(1, numalignedp);
  int count1 = 0;  // residue number in sequence 1
  int count2 = 0;  // residue number in sequence 2
  int countp = 0;  // number of aligned positions
  for (i = 0; i < strlen(alnseq1); i++) {
    if (alnseq1[i] != '-') count1++;
    if (alnseq2[i] != '-') count2++;
    if ((alnseq1[i] >= 65) && (alnseq1[i] <= 90)) {
      salpos[0][countp] = count1;
      salpos[1][countp] = count2;
      // cout << countp << " " <<  salpos[0][countp] << " " << salpos[1][countp]
      // << endl;
      countp++;
    }
  }
  // cout << "HEERElllll" <<endl;
  fclose(resultfp);
  return salpos;
}

// get aligned residue positions for two dali id names
int **get_salnpos_tmalign(char *name1, char *name2, int &numalignedp) {
  int i, j, k;
  char dali_dir[100], index_dir[100];
  strcpy(dali_dir, "/usr2/db/tmalign_1.69/");
  // strcpy(dali_dir, "/home/jpei/db/tmalign_1.69/");
  strcpy(index_dir, "/home/jpei/db/tmalign_1.69/");

  char result_file[100];
  char index_file[100];
  sprintf(result_file, "%s%s.tmalign.result", dali_dir, name1);
  sprintf(index_file, "%s%s.tmalign.result.index", index_dir, name1);

  cout << result_file << endl;
  cout << index_file << endl;

  int name2_index;
  sscanf(name2, "%d", &name2_index);
  // cout << name2_index << endl;
  int address_daliscore[2];
  FILE *idxfp = fopen(index_file, "r");
  // cout << name2_index << endl;
  if (!idxfp) cout << "not readable" << endl;
  fseek(idxfp, (name2_index)*4, 0);
  // cout << name2_index << endl;
  fread(address_daliscore, 4, 1, idxfp);
  // cout << name2_index << endl;
  cout << name1 << " " << name2 << " " << address_daliscore[0] << endl;
  fclose(idxfp);
  // if(address_daliscore[1]<minDaliZ * 10) return NULL;
  FILE *resultfp = fopen(result_file, "r");
  fseek(resultfp, address_daliscore[0], 0);

  char alnseq1[3000], alnseq2[3000];

  fgets(alnseq1, 3000, resultfp);
  // check names
  char combinednames[100];
  sprintf(combinednames, "%s %s\n", name2, name1);
  if (strcmp(alnseq1, combinednames) != 0) {
    cout << "combined names: " << combinednames << endl;
    cout << "combined names in result file: " << alnseq1 << endl;
  }
  // fgets( alnseq1, 3000, resultfp);
  fgets(alnseq2, 3000, resultfp);
  fgets(alnseq1, 3000, resultfp);

  numalignedp = 0;
  for (i = 0; i < strlen(alnseq1); i++) {
    if ((alnseq1[i] >= 65) && (alnseq1[i] <= 90)) {
      numalignedp++;
    }
  }
  // cout << "numalignedp: " << numalignedp << endl;
  // cout << alnseq1 << endl;
  // cout << alnseq2 << endl;
  int **salpos = imatrix(1, numalignedp);
  int count1 = 0;  // residue number in sequence 1
  int count2 = 0;  // residue number in sequence 2
  int countp = 0;  // number of aligned positions
  for (i = 0; i < strlen(alnseq1); i++) {
    if (alnseq1[i] != '-') count1++;
    if (alnseq2[i] != '-') count2++;
    if ((alnseq1[i] >= 65) && (alnseq1[i] <= 90)) {
      salpos[0][countp] = count1;
      salpos[1][countp] = count2;
      // cout << countp << " " <<  salpos[0][countp] << " " << salpos[1][countp]
      // << endl;
      countp++;
    }
  }
  // cout << "HEERElllll" <<endl;
  fclose(resultfp);
  return salpos;
}

// get aligned residue positions for two dali id names
int **get_salnpos_updated(seq_str_aln *a1, seq_str_aln *a2, int ind1, int ind2,
                          int &numalignedp, const char *prog_name) {
  int i, j, k;
  char dali_dir[100], index_dir[100];
  strcpy(dali_dir, "/usr2/db/tmalign_1.69/");
  // strcpy(dali_dir, "/home/jpei/db/tmalign_1.69/");
  strcpy(index_dir, "/home/jpei/db/tmalign_1.69/");
  int **salpos = NULL;

  if (debug_here > 11) {
    cout << "a1 id:    " << a1->id[ind1] << endl;
    cout << "a1 start: " << a1->str_start[ind1] << endl;
    cout << "a1 end:   " << a1->str_end[ind1] << endl;
    cout << "a2 id:    " << a2->id[ind2] << endl;
    cout << "a2 start: " << a2->str_start[ind2] << endl;
    cout << "a2 end:   " << a2->str_end[ind2] << endl;
  }

  // cout << a1->id[ind1] << " " << a2->id[ind2] << " " << prog_name << endl;
  fprintf(stdout, "# %s %s by %s\n", a1->id[ind1], a2->id[ind2], prog_name);
  // check if the pdb code ids are the same; if the same, no need to run
  // structural comparison program
  if (strcmp(a1->id[ind1], a2->id[ind2]) == 0) {
    int maxstart = MAX_HERE(a1->str_start[ind1], a2->str_start[ind2]);
    int minend = MIN_HERE(a1->str_end[ind1], a2->str_end[ind2]);
    numalignedp = minend - maxstart + 1;
    if (numalignedp <= 0) return NULL;
    salpos = imatrix(1, numalignedp);
    if (debug_here > 11) cout << "same id" << endl;
    if (debug_here > 11)
      cout << "maxstart: " << maxstart << " minend: " << minend << endl;
    for (i = maxstart; i <= minend; i++) {
      salpos[0][i - maxstart] = i;
      salpos[1][i - maxstart] = i;
      // Debug
      if (debug_here > 11)
        cout << salpos[0][i - maxstart] << " " << salpos[1][i - maxstart]
             << endl;
    }
    return salpos;
  }

  // pdb code ids are not the same, run structural comparison programs by
  // calling a python wrapper
  char command[500];
  char tmpstr[10000];
  char *ss;
  int start_recording = 0;
  // sprintf(command,
  // "/home/jpei/evolvable_database/bin/get_struct_constraint.py %s %s -%s
  // -dirname %s -start1 %d -start2 %d -end1 %d -end2 %d -clear 1", a1->id[ind1],
  // a2->id[ind2], prog_name, blast_dir, a1->str_start[ind1],
  // a2->str_start[ind2], a1->str_end[ind1], a2->str_end[ind2]);
  sprintf(command,
          "python %s/prog/saln/get_struct_constraint.py %s %s -%s -dirname %s "
          "-start1 %d -start2 %d -end1 %d -end2 %d -clear 1 -pdbdir "
          "%s/db/structure_db -progdir %s/bin",
          program_dir, a1->id[ind1], a2->id[ind2], prog_name, blast_dir,
          a1->str_start[ind1], a2->str_start[ind2], a1->str_end[ind1],
          a2->str_end[ind2], program_dir, program_dir);
  cout << "Command: " << endl;
  cout << command << endl;
  if (debug_here > 11) cout << command << endl;
  FILE *pfp = popen(command, "r");
  if (!pfp) return NULL;
  int printline = 1;
  while (fgets(tmpstr, 10000, pfp)) {
    // Debug
    if (printline) {
      cout << tmpstr;
    }
    if (strncmp(tmpstr, "aligned positions", 17) == 0) {
      // Debug
      printline = 0;
      if (debug_here > 11) cout << tmpstr << endl;
      for (ss = tmpstr; !isdigit(*ss); ss++)
        ;
      numalignedp = atoi(ss);
      salpos = imatrix(1, numalignedp);
      start_recording = 1;
      continue;
    }
    if (!start_recording) continue;
    salpos[0][start_recording - 1] = atoi(tmpstr);
    for (ss = tmpstr; isdigit(*ss); ss++)
      ;
    for (ss = ss; !isdigit(*ss); ss++)
      ;
    salpos[1][start_recording - 1] = atoi(ss);
    // Debug
    // cout << salpos[0][start_recording-1] << " " <<
    // salpos[1][start_recording-1] << endl;
    start_recording++;
  }
  pclose(pfp);
  cout << endl;

  return salpos;
}

float **map_structure_nopromals(seq_str_aln *a1, seq_str_aln *a2,
                                const char *prog_name) {
  int i, j, k;

  int **salnpos;  // structure alignment positions; two arrays of integer values
  // [0..1][0..numalignedp-1]
  // Example:
  // aCKE-EKaa
  // -CEEKEKa-
  // array 1: 2 3 4 5 6
  // array 2: 1 2 3 5 6
  // if no structural alignment with z-score > 5, salnpos = 0
  int numalignedp;

  float **tmpMat = gmatrix<float>(a1->len, a2->len);
  for (i = 1; i <= a1->len; i++)
    for (j = 1; j <= a2->len; j++) tmpMat[i][j] = 0;
  // p1 and p2 are the correponding position numbers of the queries
  int p1, p2;
  // cout << "here " << endl;
  // cout << a1->query << endl;
  // cout << a2->query << endl;
  // cout << a1->nhits << " " << a2->nhits << endl;
  ////cout << a1->id[1] << "  ----  " <<  a2->id[1] << endl;
  int tmpMat_count = 0;
  for (i = 1; i <= a1->nhits; i++) {
    for (j = 1; j <= a2->nhits; j++) {
      // cout << "HHHHH" << endl;
      // cout << a1->id[i] << "  ----  " <<  a2->id[j] << endl;
      // salnpos = get_salnpos_tmalign(a1->id[i], a2->id[j], numalignedp);
      if (strcmp(prog_name, "dali") == 0) {
        salnpos = get_salnpos(a1->id[i], a2->id[j], numalignedp);
        // cout << "here" << endl;
      } else if (strcmp(prog_name, "fast") == 0) {
        salnpos = get_salnpos_fast(a1->id[i], a2->id[j], numalignedp);
      } else if (strcmp(prog_name, "tmalign") == 0) {
        salnpos = get_salnpos_tmalign(a1->id[i], a2->id[j], numalignedp);
      }
      if (!salnpos) continue;
      // cout << numalignedp << endl;
      for (k = 0; k < numalignedp; k++) {
        // cout << "1: " << salnpos[0][k] <<  " " << a1->aln[i][salnpos[0][k]]
        // << endl;
        p1 = a1->aln[i][salnpos[0][k]];
        if (!p1) {
          // cout << "not p1" << endl;
          continue;
        }
        // cout << "2: " << endl; cout << salnpos[1][k] <<  endl; cout <<  " "
        // << a2->aln[j][salnpos[1][k]] << endl;
        p2 = a2->aln[j][salnpos[1][k]];
        if (!p2) continue;
        if (tmpMat[p1][p2] == 0) {
          tmpMat[p1][p2] = 1;
          tmpMat_count++;
        }
      }
      free_imatrix(salnpos, 1, numalignedp);
    }
  }
  // cout << "tmpMat_count: " << tmpMat_count << endl;
  if (tmpMat_count)
    return tmpMat;
  else {
    free_gmatrix<float>(tmpMat, a1->len, a2->len);
    return NULL;
  }
  // cout << "This place" << endl;
}

// prog_name: the name of the structural alignment program
float **map_structure_promals(subalign *auxa, subalign *auxb, seq_str_aln *a1,
                              seq_str_aln *a2, const char *prog_name) {
  int i, j, k, l;

  int **salnpos =
      NULL;  // structure alignment positions; two arrays of integer values
  // [0..1][0..numalignedp-1]
  // Example:
  // aCKE-EKaa
  // -CEEKEKa-
  // array 1: 2 3 4 5 6
  // array 2: 1 2 3 5 6
  // if no structural alignment with z-score > 5, salnpos = 0
  int numalignedp;

  float **tmpMat = gmatrix<float>(a1->len, a2->len);
  for (i = 1; i <= a1->len; i++)
    for (j = 1; j <= a2->len; j++) tmpMat[i][j] = 0;
  // p1 and p2 are the correponding position numbers of the queries
  int p1, p2;
  // cout << "here " << endl;
  // cout << a1->query << endl;
  // cout << a2->query << endl;
  // cout << a1->nhits << " " << a2->nhits << endl;
  ////cout << a1->id[1] << "  ----  " <<  a2->id[1] << endl;
  int tmpMat_count = 0;

  float **realmat;
  int lenx, leny;

  for (i = 1; i <= a1->nhits; i++) {
    for (j = 1; j <= a2->nhits; j++) {
      // 1. get the structural alignment positions
      // cout << "here" << endl;
      if (strcmp(prog_name, "dali") == 0) {
        salnpos = get_salnpos(a1->id[i], a2->id[j], numalignedp);
        // cout << "here" << endl;
      } else if (strcmp(prog_name, "fast") == 0) {
        salnpos = get_salnpos_fast(a1->id[i], a2->id[j], numalignedp);
      } else if (strcmp(prog_name, "tmalign") == 0) {
        salnpos = get_salnpos_tmalign(a1->id[i], a2->id[j], numalignedp);
      }
      if (!salnpos) continue;
      // cout << "here" << endl;

      // 2. transform the salnpos into a sparse matrix
      // cout << "slen size: " << a1->slen.size() << endl;
      // for(k=1;k<a1->slen.size();k++) {cout << a1->slen[k] << endl;}
      float **m_bc = gmatrix<float>(a1->slen[i], a2->slen[j]);
      // cout << "here1" << endl;
      for (k = 1; k <= a1->slen[i]; k++)
        for (l = 1; l <= a2->slen[j]; l++) m_bc[k][l] = 0;
      // cout << "here" << endl;
      // cout << a1->query << " " << a1->id[i] << endl;
      // cout << a2->query << " " << a2->id[j] << endl;
      for (k = 0; k < numalignedp; k++) {
        // cout << a1->slen[i] << " " << a2->slen[j] << " " << salnpos[0][k] <<
        // " " << salnpos[1][k] << endl;
        m_bc[salnpos[0][k]][salnpos[1][k]] = 1.0;
      }
      // cout << "here" << endl;
      sparseMatrix *m_bc_S = new sparseMatrix(m_bc, a1->slen[i], a2->slen[j]);
      // cout << "here" << endl;
      free_gmatrix<float>(m_bc, a1->slen[i], a2->slen[j]);
      // cout << "here" << endl;

      // 3. get the residue match probability matrix between target a and
      // template b
      // hmm_psipred *hmmProfPair = new hmm_psipred(auxa, auxb);
      hmm_psipred *hmmProfPair = new hmm_psipred(auxa, a1->prof[i]);
      // cout << "here" << endl;
      hmmProfPair->set_parameters(params1);
      // cout << "here" << endl;
      hmmProfPair->get_scores(ss_w, score_w);
      // cout << "here" << endl;
      hmmProfPair->forward1();
      hmmProfPair->backward1();
      realmat = hmmProfPair->probMat;
      lenx = hmmProfPair->lenx;
      leny = hmmProfPair->leny;

      for (k = 1; k <= lenx; k++) {
        for (l = 1; l <= leny; l++) {
          if (realmat[k][l] < minProb) realmat[k][l] = 0;
        }
      }
      sparseMatrix *m_ab_S = new sparseMatrix(realmat, lenx, leny);
      delete hmmProfPair;
      // cout << "here" << endl;

      // 4. get the residue match probability matrix between template c and
      // target d
      hmmProfPair = new hmm_psipred(a2->prof[j], auxb);
      hmmProfPair->set_parameters(params1);
      hmmProfPair->get_scores(ss_w, score_w);
      hmmProfPair->forward1();
      hmmProfPair->backward1();
      realmat = hmmProfPair->probMat;
      lenx = hmmProfPair->lenx;
      leny = hmmProfPair->leny;

      for (k = 1; k <= lenx; k++) {
        for (l = 1; l <= leny; l++) {
          if (realmat[k][l] < minProb) realmat[k][l] = 0;
        }
      }
      sparseMatrix *m_cd_S = new sparseMatrix(realmat, lenx, leny);
      delete hmmProfPair;
      // cout << "here" << endl;

      // 5. multiply all three matrices
      lenx = auxa->alilen;
      leny = a2->prof[j]->alilen;
      // cout << "lens: " << lenx << " " << leny << endl;
      realmat = gmatrix<float>(lenx, leny);
      for (k = 1; k <= lenx; k++)
        for (l = 1; l <= leny; l++) realmat[k][l] = 0;
      relaxTwoSparse(m_ab_S, m_bc_S, realmat);
      // cout << "here___" << endl;
      sparseMatrix *m_abc_S =
          new sparseMatrix(realmat, lenx, a2->prof[j]->alilen);
      // cout << "here___" << endl;
      free_gmatrix(realmat, lenx, leny);
      leny = auxb->alilen;
      realmat = gmatrix<float>(lenx, leny);
      for (k = 1; k <= lenx; k++)
        for (l = 1; l <= leny; l++) realmat[k][l] = 0;
      relaxTwoSparse(m_abc_S, m_cd_S, realmat);
      // for(k=1;k<=lenx;k++) for(l=1;l<=leny;l++) if(realmat[k][l]<minProb)
      // realmat[k][l]=0; cout << "here___" << endl;

      // 6. add realmat to tmpMat
      for (k = 1; k <= lenx; k++)
        for (l = 1; l <= leny; l++) tmpMat[k][l] += realmat[k][l];
      free_gmatrix(realmat, lenx, leny);
      // cout << "here___=======" << endl;

      // 7. clean up the memories
      delete m_ab_S;
      // cout << "here___*******" << endl;
      delete m_bc_S;
      // cout << "here___*******" << endl;
      delete m_cd_S;
      // cout << "here___*******" << endl;
      delete m_abc_S;
      // cout << "here___*******" << endl;
      free_imatrix(salnpos, 1, numalignedp);
      // cout << "here___*******" << endl;
    }
  }
  // filter the matrix with minProb
  for (i = 1; i <= a1->len; i++) {
    for (j = 1; j <= a2->len; j++) {
      if (tmpMat[i][j] < minProb)
        tmpMat[i][j] = 0;
      else {
        // cout << "tmpMat: " << i << " " << j << " " << auxa->aseq[0][i-1] << "
        // " << auxb->aseq[0][j-1]<< tmpMat[i][j] << endl;
        tmpMat_count++;
      }
    }
  }

  // cout << "tmpMat_count: " << tmpMat_count << endl;
  if (tmpMat_count)
    return tmpMat;
  else {
    free_gmatrix<float>(tmpMat, a1->len, a2->len);
    return NULL;
  }
  // cout << "This place" << endl;
}

float **map_structure_nopromals_updated(seq_str_aln *a1, seq_str_aln *a2,
                                        const char *prog_name) {
  int i, j, k, l;

  int **salnpos;  // structure alignment positions; two arrays of integer values
  // [0..1][0..numalignedp-1]
  // Example:
  // aCKE-EKaa
  // -CEEKEKa-
  // array 1: 2 3 4 5 6
  // array 2: 1 2 3 5 6
  // if no structural alignment with z-score > 5, salnpos = 0
  int numalignedp;

  float **tmpMat = gmatrix<float>(a1->len, a2->len);
  for (i = 1; i <= a1->len; i++)
    for (j = 1; j <= a2->len; j++) tmpMat[i][j] = 0;
  // p1 and p2 are the correponding position numbers of the queries
  int p1, p2;
  // cout << "here " << endl;
  // cout << a1->query << endl;
  // cout << a2->query << endl;
  // cout << a1->nhits << " " << a2->nhits << endl;
  ////cout << a1->id[1] << "  ----  " <<  a2->id[1] << endl;
  int tmpMat_count = 0;
  for (i = 1; i <= a1->nhits; i++) {
    for (j = 1; j <= a2->nhits; j++) {
      // Debug here
      // cout << "HHHHH" << endl;
      // cout << a1->id[i] << "  ----  " <<  a2->id[j] << endl;
      // salnpos = get_salnpos_tmalign(a1->id[i], a2->id[j], numalignedp);
      salnpos = get_salnpos_updated(a1, a2, i, j, numalignedp, prog_name);
      /*
      if(strcmp(prog_name, "dali")==0) {
              salnpos = get_salnpos(a1->id[i], a2->id[j], numalignedp);
      //cout << "here" << endl;
      }
      else if(strcmp(prog_name, "fast")==0) {
              salnpos = get_salnpos_fast(a1->id[i], a2->id[j], numalignedp);
      }
      else if(strcmp(prog_name, "tmalign")==0) {
              salnpos = get_salnpos_tmalign(a1->id[i], a2->id[j], numalignedp);
      }
      */
      if (!salnpos) {
        if (strcmp("dali", prog_name) == 0) {
          cout << "Using_fast_instead" << endl;
          salnpos = get_salnpos_updated(a1, a2, i, j, numalignedp, "fast");
          if (!salnpos) continue;
        } else
          continue;
      }
      // derived alignment
      string s1, s2, tmpa1, tmpa2;
      int prev_p1 = 0, prev_p2 = 0;
      int k1 = 0;
      // Debug here
      if (debug_here > 11) cout << "numalignedp: " << numalignedp << endl;
      for (k = 0; k < numalignedp; k++) {
        if (debug_here > 11)
          cout << "1: " << salnpos[0][k] << " " << a1->aln[i][salnpos[0][k]]
               << endl;
        p1 = a1->aln[i][salnpos[0][k]];
        if (!p1) {
          // cout << "not p1" << endl;
          continue;
        }
        // cout << a2->id[j] << " " <<  a2->slen[j] << " " <<
        // a2->aln[j][a2->aln[j][a2->slen[j]-1]] << endl;
        // for(l=1;l<=a2->slen[j];l++) { cout << "l: " << l << " " <<
        // a2->aln[j][l] << endl; } cout << "2: " << endl; cout << salnpos[1][k]
        // <<  endl; cout <<  " " << a2->aln[j][salnpos[1][k]] << endl;
        p2 = a2->aln[j][salnpos[1][k]];
        if (!p2) continue;
        // add here to print out the derived sequence
        tmpa1 = "";
        tmpa2 = "";
        if (k1 == 0) {
          cout << p1 << " " << p2 << endl;
        }
        k1++;
        for (l = prev_p1; l < p1 - 1; l++) {
          tmpa1 += tolower(a1->query[l]);
        }
        for (l = prev_p2; l < p2 - 1; l++) {
          tmpa2 += tolower(a2->query[l]);
        }
        for (l = 1; l <= p1 - prev_p1 + prev_p2 - p2; l++) {
          tmpa2 += '-';
        }
        for (l = 1; l <= p2 - prev_p2 + prev_p1 - p1; l++) {
          tmpa1 += '-';
        }
        s1 += tmpa1;
        s2 += tmpa2;
        s1 += a1->query[p1 - 1];
        s2 += a2->query[p2 - 1];
        if (tmpMat[p1][p2] == 0) {
          tmpMat[p1][p2] = 1;
          tmpMat_count++;
        }
        prev_p1 = p1;
        prev_p2 = p2;
      }
      cout << s1 << endl;
      cout << s2 << endl;
      cout << endl;
      // Debug here
      // cout << "Here" << endl;
      free_imatrix(salnpos, 1, numalignedp);
    }
  }
  // cout << "tmpMat_count: " << tmpMat_count << endl;
  if (tmpMat_count)
    return tmpMat;
  else {
    free_gmatrix<float>(tmpMat, a1->len, a2->len);
    return NULL;
  }
  // cout << "This place" << endl;
}

// prog_name: the name of the structural alignment program
float **map_structure_promals_updated(subalign *auxa, subalign *auxb,
                                      seq_str_aln *a1, seq_str_aln *a2,
                                      int dali, int fast, int tmalign) {
  int i, j, k, l;

  int **salnpos =
      NULL;  // structure alignment positions; two arrays of integer values
  int **salnpos1 = NULL;
  int **salnpos2 = NULL;
  // [0..1][0..numalignedp-1]
  // Example:
  // aCKE-EKaa
  // -CEEKEKa-
  // array 1: 2 3 4 5 6
  // array 2: 1 2 3 5 6
  // if no structural alignment with z-score > 5, salnpos = 0
  int numalignedp, numalignedp1, numalignedp2;

  float **tmpMat = gmatrix<float>(a1->len, a2->len);
  for (i = 1; i <= a1->len; i++)
    for (j = 1; j <= a2->len; j++) tmpMat[i][j] = 0;
  // p1 and p2 are the correponding position numbers of the queries
  int p1, p2;
  // cout << "here " << endl;
  // cout << a1->query << endl;
  // cout << a2->query << endl;
  // cout << a1->nhits << " " << a2->nhits << endl;
  ////cout << a1->id[1] << "  ----  " <<  a2->id[1] << endl;
  int tmpMat_count = 0;

  float **realmat;
  int lenx, leny;
  int offset_b, offset_c;
  // cout << "nhits1: " << a1->nhits << " nhits2: " << a2->nhits << endl;
  fprintf(stdout, "### %s %s\n", a1->query_name, a2->query_name);

  for (i = 1; i <= a1->nhits; i++) {
    for (j = 1; j <= a2->nhits; j++) {
      // 1. get the structural alignment positions
      // cout << "here" << endl;
      if (dali)
        salnpos = get_salnpos_updated(a1, a2, i, j, numalignedp, "dali");
      if (fast)
        salnpos1 = get_salnpos_updated(a1, a2, i, j, numalignedp1, "fast");
      if (tmalign)
        salnpos2 = get_salnpos_updated(a1, a2, i, j, numalignedp2, "tmalign");
      if ((!salnpos) && (!salnpos1) && (!salnpos2)) continue;

      // 1'. check dimensions
      offset_b = 0;
      offset_c = 0;
      lenx = a1->prof[i]->alilen;
      leny = a2->prof[j]->alilen;
      if (a1->slen[i] != lenx) {
        offset_b = a1->str_start[i] - 1;
      }
      if (a2->slen[j] != leny) {
        offset_c = a2->str_start[j] - 1;
      }

      // 2. transform the salnpos into a sparse matrix
      // cout << "slen size: " << a1->slen.size() << endl;
      // for(k=1;k<a1->slen.size();k++) {cout << a1->slen[k] << endl;}
      // float **m_bc = gmatrix<float>(a1->slen[i], a2->slen[j]);
      float **m_bc = gmatrix<float>(lenx, leny);
      // cout << "here1" << endl;
      for (k = 1; k <= lenx; k++)
        for (l = 1; l <= leny; l++) m_bc[k][l] = 0;
      // cout << "here" << endl;
      // cout << a1->query << " " << a1->id[i] << endl;
      // cout << a2->query << " " << a2->id[j] << endl;
      // cout << "numalignedp: " << numalignedp << endl;
      if (salnpos) {
        for (k = 0; k < numalignedp; k++) {
          // cout << k << " " << a1->slen[i] << " " << a2->slen[j] << " " <<
          // salnpos[0][k] << " " << salnpos[1][k] << endl;
          m_bc[salnpos[0][k] - offset_b][salnpos[1][k] - offset_c] += 1.0;
        }
      }
      if (salnpos1) {
        for (k = 0; k < numalignedp1; k++) {
          // cout << k << " " << a1->slen[i] << " " << a2->slen[j] << " " <<
          // salnpos[0][k] << " " << salnpos[1][k] << endl;
          cout << k << " " << salnpos1[0][k] - offset_b << " "
               << salnpos1[1][k] - offset_c << endl;
          m_bc[salnpos1[0][k] - offset_b][salnpos1[1][k] - offset_c] += 1.0;
        }
      }
      if (salnpos2) {
        for (k = 0; k < numalignedp2; k++) {
          // cout << k << " " << a1->slen[i] << " " << a2->slen[j] << " " <<
          // salnpos[0][k] << " " << salnpos[1][k] << endl;
          cout << k << " " << salnpos2[0][k] - offset_b << " "
               << salnpos2[1][k] - offset_c << endl;
          m_bc[salnpos2[0][k] - offset_b][salnpos2[1][k] - offset_c] += 1.0;
        }
      }
      // check m_bc
      if (debug > 1) {
        fprintf(stdout, "m_bc matrix:\n");
        for (k = 1; k <= lenx; k++) {
          for (l = 1; l <= leny; l++) {
            if (m_bc[k][l]) fprintf(stdout, "%d %d %3.2f\n", k, l, m_bc[k][l]);
          }
        }
      }

      // cout << "b dimension: " << a1->slen[i] << endl;
      // cout << "c dimension: " << a1->slen[i] << endl;
      sparseMatrix *m_bc_S = new sparseMatrix(m_bc, lenx, leny);
      free_gmatrix<float>(m_bc, lenx, leny);

      // 3. get the residue match probability matrix between target a and
      // template b
      // hmm_psipred *hmmProfPair = new hmm_psipred(auxa, auxb);
      // cout << "a1 prof: " << endl; a1->prof[i]->printali(70);
      hmm_psipred *hmmProfPair = new hmm_psipred(auxa, a1->prof[i]);
      hmmProfPair->set_parameters(params1);
      hmmProfPair->get_scores(ss_w, score_w);
      hmmProfPair->forward1();
      hmmProfPair->backward1();
      // weight_end_penalty = 0.001;
      // hmmProfPair->forward_no_end_penalty();
      // hmmProfPair->backward_no_end_penalty();
      realmat = hmmProfPair->probMat;
      lenx = hmmProfPair->lenx;
      leny = hmmProfPair->leny;

      for (k = 1; k <= lenx; k++) {
        for (l = 1; l <= leny; l++) {
          if (realmat[k][l] < minProb) realmat[k][l] = 0;
        }
      }
      sparseMatrix *m_ab_S = new sparseMatrix(realmat, lenx, leny);
      if (debug > 1) m_ab_S->printCrs();
      delete hmmProfPair;

      // cout << "b dimension 2: " << leny << endl;
      // cout << "b start: " << a1->str_start[i] << endl;
      // cout << "a dimension: " << lenx << endl;

      // 4. get the residue match probability matrix between template c and
      // target d
      hmmProfPair = new hmm_psipred(a2->prof[j], auxb);
      hmmProfPair->set_parameters(params1);
      hmmProfPair->get_scores(ss_w, score_w);
      hmmProfPair->forward1();
      hmmProfPair->backward1();
      // weight_end_penalty = 0.001;
      // hmmProfPair->forward_no_end_penalty();
      // hmmProfPair->backward_no_end_penalty();
      realmat = hmmProfPair->probMat;
      lenx = hmmProfPair->lenx;
      leny = hmmProfPair->leny;

      for (k = 1; k <= lenx; k++) {
        for (l = 1; l <= leny; l++) {
          if (realmat[k][l] < minProb) realmat[k][l] = 0;
        }
      }
      sparseMatrix *m_cd_S = new sparseMatrix(realmat, lenx, leny);
      if (debug > 1) m_cd_S->printCrs();
      delete hmmProfPair;

      // 5. multiply all three matrices
      lenx = auxa->alilen;
      leny = a2->prof[j]->alilen;
      realmat = gmatrix<float>(lenx, leny);
      for (k = 1; k <= lenx; k++)
        for (l = 1; l <= leny; l++) realmat[k][l] = 0;
      relaxTwoSparse(m_ab_S, m_bc_S, realmat);
      sparseMatrix *m_abc_S = new sparseMatrix(realmat, lenx, leny);
      free_gmatrix(realmat, lenx, leny);
      leny = auxb->alilen;
      realmat = gmatrix<float>(lenx, leny);
      for (k = 1; k <= lenx; k++)
        for (l = 1; l <= leny; l++) realmat[k][l] = 0;
      relaxTwoSparse(m_abc_S, m_cd_S, realmat);
      // for(k=1;k<=lenx;k++) for(l=1;l<=leny;l++) if(realmat[k][l]<minProb)
      // realmat[k][l]=0;

      // 6. add realmat to tmpMat
      for (k = 1; k <= lenx; k++)
        for (l = 1; l <= leny; l++) tmpMat[k][l] += realmat[k][l];
      free_gmatrix(realmat, lenx, leny);

      // 7. clean up the memories
      delete m_ab_S;
      // cout << "here___*******" << endl;
      delete m_bc_S;
      // cout << "here___*******" << endl;
      delete m_cd_S;
      // cout << "here___*******" << endl;
      delete m_abc_S;
      // cout << "here___*******" << endl;
      free_imatrix(salnpos, 1, numalignedp);
    }
  }
  // filter the matrix with minProb
  for (i = 1; i <= a1->len; i++) {
    for (j = 1; j <= a2->len; j++) {
      if (tmpMat[i][j] < minProb)
        tmpMat[i][j] = 0;
      else {
        // cout << "tmpMat: " << i << " " << j << " " << auxa->aseq[0][i-1] << "
        // " << auxb->aseq[0][j-1]<< tmpMat[i][j] << endl;
        tmpMat_count++;
      }
    }
  }
  cout << endl;

  // cout << "tmpMat_count: " << tmpMat_count << endl;
  if (tmpMat_count)
    return tmpMat;
  else {
    free_gmatrix<float>(tmpMat, a1->len, a2->len);
    return NULL;
  }
  // cout << "This place" << endl;
}

// prog_name: the name of the structural alignment program
/*
float **map_structure_promals_updated(subalign *auxa, subalign *auxb,
seq_str_aln *a1, seq_str_aln *a2, char *prog_name) {

        int i, j, k, l;

        int **salnpos=NULL; // structure alignment positions; two arrays of
integer values
        // [0..1][0..numalignedp-1]
        // Example:
        // aCKE-EKaa
        // -CEEKEKa-
        // array 1: 2 3 4 5 6
        // array 2: 1 2 3 5 6
        // if no structural alignment with z-score > 5, salnpos = 0
        int numalignedp;

        float **tmpMat = gmatrix<float>(a1->len, a2->len);
        for(i=1;i<=a1->len;i++) for(j=1;j<=a2->len;j++) tmpMat[i][j] = 0;
        // p1 and p2 are the correponding position numbers of the queries
        int p1, p2;
        //cout << "here " << endl;
        //cout << a1->query << endl;
        //cout << a2->query << endl;
        //cout << a1->nhits << " " << a2->nhits << endl;
        ////cout << a1->id[1] << "  ----  " <<  a2->id[1] << endl;
        int tmpMat_count = 0;

        float **realmat;
        int lenx, leny;
        int offset_b, offset_c;
        //cout << "nhits1: " << a1->nhits << " nhits2: " << a2->nhits << endl;
        fprintf(stdout, "%s %s\n", a1->query_name, a2->query_name);

        for(i=1;i<=a1->nhits;i++) {
                for(j=1;j<=a2->nhits;j++) {

                        // 1. get the structural alignment positions
                        //cout << "here" << endl;
                        salnpos = get_salnpos_updated(a1, a2, i, j, numalignedp,
prog_name); if(!salnpos) continue;

                        // 1'. check dimensions
                        offset_b = 0; offset_c=0;
                        lenx = a1->prof[i]->alilen;
                        leny = a2->prof[j]->alilen;
                        if(a1->slen[i]!=lenx) { offset_b = a1->str_start[i]-1; }
                        if(a2->slen[j]!=leny) { offset_c = a2->str_start[j]-1; }

                        // 2. transform the salnpos into a sparse matrix
                        //cout << "slen size: " << a1->slen.size() << endl;
                        //for(k=1;k<a1->slen.size();k++) {cout << a1->slen[k] <<
endl;}
                        //float **m_bc = gmatrix<float>(a1->slen[i],
a2->slen[j]); float **m_bc = gmatrix<float>(lenx, leny);
                        //cout << "here1" << endl;
                        for(k=1;k<=lenx;k++) for(l=1;l<=leny;l++) m_bc[k][l] =
0;
                        //cout << "here" << endl;
                        //cout << a1->query << " " << a1->id[i] << endl;
                        //cout << a2->query << " " << a2->id[j] << endl;
                        //cout << "numalignedp: " << numalignedp << endl;
                        for(k=0;k<numalignedp;k++) {
                                //cout << k << " " << a1->slen[i] << " " <<
a2->slen[j] << " " << salnpos[0][k] << " " << salnpos[1][k] << endl;
                                m_bc[salnpos[0][k]-offset_b][salnpos[1][k]-offset_c]
= 1.0;
                        }
                        //cout << "b dimension: " << a1->slen[i] << endl;
                        //cout << "c dimension: " << a1->slen[i] << endl;
                        sparseMatrix *m_bc_S = new sparseMatrix(m_bc, lenx,
leny); free_gmatrix<float>(m_bc, lenx, leny);

                        // 3. get the residue match probability matrix between
target a and template b
                        //hmm_psipred *hmmProfPair = new hmm_psipred(auxa,
auxb);
                        //cout << "a1 prof: " << endl;
a1->prof[i]->printali(70); hmm_psipred *hmmProfPair = new hmm_psipred(auxa,
a1->prof[i]); hmmProfPair->set_parameters(params1);
                        hmmProfPair->get_scores(ss_w, score_w);
                        hmmProfPair->forward1();
                        hmmProfPair->backward1();
                        //weight_end_penalty = 0.001;
                        //hmmProfPair->forward_no_end_penalty();
                        //hmmProfPair->backward_no_end_penalty();
                        realmat = hmmProfPair->probMat;
                        lenx = hmmProfPair->lenx; leny = hmmProfPair->leny;

                        for(k=1;k<=lenx;k++){for(l=1;l<=leny;l++){if(realmat[k][l]<minProb)realmat[k][l]
= 0;}} sparseMatrix *m_ab_S = new sparseMatrix(realmat, lenx, leny); if(debug>1)
m_ab_S->printCrs(); delete hmmProfPair;

                        //cout << "b dimension 2: " << leny << endl;
                        //cout << "b start: " << a1->str_start[i] << endl;
                        //cout << "a dimension: " << lenx << endl;

                        // 4. get the residue match probability matrix between
template c and target d hmmProfPair = new hmm_psipred(a2->prof[j], auxb);
                        hmmProfPair->set_parameters(params1);
                        hmmProfPair->get_scores(ss_w, score_w);
                        hmmProfPair->forward1();
                        hmmProfPair->backward1();
                        //weight_end_penalty = 0.001;
                        //hmmProfPair->forward_no_end_penalty();
                        //hmmProfPair->backward_no_end_penalty();
                        realmat = hmmProfPair->probMat;
                        lenx = hmmProfPair->lenx; leny = hmmProfPair->leny;

                        for(k=1;k<=lenx;k++){for(l=1;l<=leny;l++){if(realmat[k][l]<minProb)realmat[k][l]
= 0;}} sparseMatrix *m_cd_S = new sparseMatrix(realmat, lenx, leny); if(debug>1)
m_cd_S->printCrs(); delete hmmProfPair;

                        // 5. multiply all three matrices
                        lenx = auxa->alilen;
                        leny = a2->prof[j]->alilen;
                        realmat = gmatrix<float>(lenx, leny);
                        for(k=1;k<=lenx;k++) for(l=1;l<=leny;l++) realmat[k][l]
= 0; relaxTwoSparse(m_ab_S, m_bc_S, realmat); sparseMatrix *m_abc_S = new
sparseMatrix(realmat, lenx, leny); free_gmatrix(realmat, lenx, leny); leny =
auxb->alilen; realmat = gmatrix<float>(lenx, leny); for(k=1;k<=lenx;k++)
for(l=1;l<=leny;l++) realmat[k][l] = 0; relaxTwoSparse(m_abc_S, m_cd_S,
realmat);
                        //for(k=1;k<=lenx;k++) for(l=1;l<=leny;l++)
if(realmat[k][l]<minProb) realmat[k][l]=0;

                        // 6. add realmat to tmpMat
                        for(k=1;k<=lenx;k++) for(l=1;l<=leny;l++) tmpMat[k][l]
+= realmat[k][l]; free_gmatrix(realmat, lenx, leny);

                        // 7. clean up the memories
                        delete m_ab_S;
                        //cout << "here___*******" << endl;
                        delete m_bc_S;
                        //cout << "here___*******" << endl;
                        delete m_cd_S;
                        //cout << "here___*******" << endl;
                        delete m_abc_S;
                        //cout << "here___*******" << endl;
                        free_imatrix(salnpos, 1, numalignedp);
                }
        }
        // filter the matrix with minProb
        for(i=1;i<=a1->len;i++) {
                for(j=1;j<=a2->len;j++) {
                        if (tmpMat[i][j] < minProb) tmpMat[i][j]=0;
                        else {
                                //cout << "tmpMat: " << i << " " << j << " " <<
auxa->aseq[0][i-1] << " " << auxb->aseq[0][j-1]<< tmpMat[i][j] << endl;
                                tmpMat_count++;
                        }
                }
        }

        //cout << "tmpMat_count: " << tmpMat_count << endl;
        if(tmpMat_count) return tmpMat;
        else {
                free_gmatrix<float>(tmpMat, a1->len, a2->len);
                return NULL;
        }
        //cout << "This place" << endl;
}
*/

// merge a and b according to the sequence with seq_name
// b is added to a
subalign *merge_align_by_one_sequence(subalign *a, subalign *b,
                                      char *seq_name) {
  int i, j, k, l;

  int ia = -1, ib = -1;

  for (i = 0; i < a->nal; i++) {
    if (strcmp(seq_name, a->aname[i]) == 0) {
      ia = i;
      break;
    }
  }
  if (ia < 0) {
    cout << "Name " << seq_name << " is not present in a" << endl;
    exit(0);
  }
  for (i = 0; i < b->nal; i++) {
    if (strcmp(seq_name, b->aname[i]) == 0) {
      ib = i;
      break;
    }
  }
  if (ib < 0) {
    cout << "Name " << seq_name << " is not present in b" << endl;
    exit(0);
  }

  // determine the number of gaps of the linking sequence in a, and in b
  int ngapa = 0, ngapb = 0;
  char *no_gap_seqa = cvector(a->alilen);
  char *no_gap_seqb = cvector(b->alilen);
  int tmp_index = 0;
  for (i = 0; i < a->alilen; i++) {
    if (a->aseq[ia][i] != '-') {
      no_gap_seqa[tmp_index] = a->aseq[ia][i];
      tmp_index++;
    } else {
      ngapa += 1;
    }
  }
  no_gap_seqa[tmp_index] = '\0';

  if (debug_here > 11) cout << "ngapa: " << ngapa << endl;

  tmp_index = 0;
  for (i = 0; i < b->alilen; i++) {
    if (b->aseq[ib][i] != '-') {
      no_gap_seqb[tmp_index] = b->aseq[ib][i];
      tmp_index++;
    } else {
      ngapb += 1;
    }
  }
  no_gap_seqb[tmp_index] = '\0';

  if (debug_here > 11) cout << "ngapb: " << ngapb << endl;
  if (debug_here > 11) cout << "no_gap_seqb: " << no_gap_seqb << endl;

  int none_gap_len = strlen(no_gap_seqb);

  if (strcmp(no_gap_seqa, no_gap_seqb) != 0) {
    cout << "None-gapped sequences are different in subalign a and subalign b"
         << endl;
    cout << no_gap_seqa << endl;
    cout << no_gap_seqb << endl;
    exit(0);
  }

  // find the gap pattern arrays for a and b
  int *gap_pattern_a = ivector(none_gap_len);
  int *gap_pattern_b = ivector(none_gap_len);

  int tmp_gap_count = 0;
  tmp_index = 0;
  for (i = 0; i < a->alilen; i++) {
    if (a->aseq[ia][i] != '-') {
      gap_pattern_a[tmp_index] = tmp_gap_count;
      if (debug > 1)
        cout << tmp_index << " " << gap_pattern_a[tmp_index] << endl;
      tmp_index += 1;
      tmp_gap_count = 0;
    } else
      tmp_gap_count += 1;
  }
  gap_pattern_a[tmp_index] = tmp_gap_count;
  if (debug > 1) cout << tmp_index << " " << gap_pattern_a[tmp_index] << endl;

  tmp_gap_count = 0;
  tmp_index = 0;
  for (i = 0; i < b->alilen; i++) {
    if (b->aseq[ib][i] != '-') {
      gap_pattern_b[tmp_index] = tmp_gap_count;
      if (debug > 1)
        cout << tmp_index << " " << gap_pattern_b[tmp_index] << endl;
      tmp_index += 1;
      tmp_gap_count = 0;
    } else
      tmp_gap_count += 1;
  }
  gap_pattern_b[tmp_index] = tmp_gap_count;
  if (debug > 1) cout << tmp_index << " " << gap_pattern_b[tmp_index] << endl;

  // set up a new align
  subalign *new_aln = new subalign();
  new_aln->nal = a->nal + b->nal - 1;
  new_aln->alilen = strlen(no_gap_seqb) + ngapa + ngapb;
  new_aln->mnamelen = 0;
  for (i = 0; i < a->nal; i++) {
    if (strlen(a->aname[i]) > new_aln->mnamelen)
      new_aln->mnamelen = strlen(a->aname[i]);
  }
  for (i = 0; i < b->nal; i++) {
    if (strlen(b->aname[i]) > new_aln->mnamelen)
      new_aln->mnamelen = strlen(b->aname[i]);
  }
  if (debug_here > 11) {
    cout << "new_aln: " << new_aln->nal << " " << new_aln->alilen << " "
         << new_aln->mnamelen << endl;
  }

  new_aln->aseq = cmatrix(new_aln->nal, new_aln->alilen + 1);
  new_aln->aname = cmatrix(new_aln->nal, new_aln->mnamelen + 1);

  for (i = 0; i < a->nal; i++) {
    strcpy(new_aln->aname[i], a->aname[i]);
  }
  for (i = 0; i < b->nal; i++) {
    if (i == ib) continue;
    if (i < ib) strcpy(new_aln->aname[i + a->nal], b->aname[i]);
    if (i > ib) strcpy(new_aln->aname[i + a->nal - 1], b->aname[i]);
  }

  // generate new sequences
  // atmp_len: temporary length of the new sequences
  // tmp_index: the index of non-gapped residues in the linking sequence
  int atmp_len = 0;
  tmp_index = 0;
  for (i = 0; i < a->nal; i++) {
    if (i != ia) continue;
    for (j = 0; j < a->alilen; j++) {
      if (a->aseq[i][j] != '-') {
        for (l = 1; l <= gap_pattern_b[tmp_index]; l++) {
          for (k = 0; k < a->nal; k++) {
            new_aln->aseq[k][atmp_len] = '-';
          }
          atmp_len++;
        }
        for (k = 0; k < a->nal; k++) {
          new_aln->aseq[k][atmp_len] = a->aseq[k][j];
        }
        atmp_len++;
        tmp_index += 1;
      } else {
        for (k = 0; k < a->nal; k++) {
          new_aln->aseq[k][atmp_len] = a->aseq[k][j];
        }
        atmp_len++;
      }
    }
    // gaps at the C-terminal
    for (l = 1; l <= gap_pattern_b[tmp_index]; l++) {
      for (k = 0; k < a->nal; k++) {
        new_aln->aseq[k][atmp_len] = '-';
      }
      atmp_len++;
    }
    // closing the sequences
    for (k = 0; k < a->nal; k++) {
      new_aln->aseq[k][atmp_len] = '\0';
    }
  }

  atmp_len = 0;
  tmp_index = 0;
  for (i = 0; i < b->nal; i++) {
    if (i != ib) continue;
    for (j = 0; j < b->alilen; j++) {
      if (b->aseq[i][j] != '-') {
        for (l = 1; l <= gap_pattern_a[tmp_index]; l++) {
          for (k = 0; k < b->nal; k++) {
            // skipping the linking sequence
            if (k == ib) continue;
            if (k < ib) new_aln->aseq[k + a->nal][atmp_len] = '-';
            if (k > ib) new_aln->aseq[k + a->nal - 1][atmp_len] = '-';
          }
          atmp_len++;
        }
        for (k = 0; k < b->nal; k++) {
          if (k == ib) continue;
          if (k < ib) new_aln->aseq[k + a->nal][atmp_len] = b->aseq[k][j];
          if (k > ib) new_aln->aseq[k + a->nal - 1][atmp_len] = b->aseq[k][j];
        }
        atmp_len++;
        tmp_index += 1;
      } else {
        for (k = 0; k < b->nal; k++) {
          if (k == ib) continue;
          if (k < ib) new_aln->aseq[k + a->nal][atmp_len] = b->aseq[k][j];
          if (k > ib) new_aln->aseq[k + a->nal - 1][atmp_len] = b->aseq[k][j];
        }
        atmp_len++;
      }
    }
    for (l = 1; l <= gap_pattern_a[tmp_index]; l++) {
      for (k = 0; k < b->nal; k++) {
        if (k == ib) continue;
        if (k < ib) new_aln->aseq[k + a->nal][atmp_len] = '-';
        if (k > ib) new_aln->aseq[k + a->nal - 1][atmp_len] = '-';
      }
      atmp_len++;
    }
    for (k = 0; k < b->nal; k++) {
      if (k == ib) continue;
      if (k < ib) new_aln->aseq[k + a->nal][atmp_len] = '\0';
      if (k > ib) new_aln->aseq[k + a->nal - 1][atmp_len] = '\0';
    }
  }

  if (debug > 1) {
    a->printali(80);
    b->printali(80);
    new_aln->printali(80);
  }

  delete[] gap_pattern_a;
  delete[] gap_pattern_b;

  return new_aln;
}

// rep_index_in_b: repre
subalign *merge_master_and_slaves(subalign *a, vector<subalign *> &b,
                                  vector<char *> &seq_name) {
  int i, j, k, l;

  int *is_rep_index_in_a, *rep_index_in_a;
  int *rep_index_in_b;

  int ns = b.size();

  if (ns == 0) return NULL;

  // 1. determine rep_index_in_b
  int found;
  rep_index_in_b = ivector(ns);
  for (i = 0; i < ns; i++) {
    found = 0;
    for (j = 0; j < b[i]->nal; j++) {
      if (strcmp(seq_name[i], b[i]->aname[j]) == 0) {
        rep_index_in_b[i] = j;
        found = 1;
        if (debug > 1) cout << "rep_index_in_b: " << rep_index_in_b[i] << endl;
        break;
      }
    }
    if (found == 0) {
      cout << "Error: cannot find the representative name: " << seq_name[i]
           << endl;
      b[i]->printali(80);
      exit(0);
    }
  }

  // 2. determine rep_index_in_a and is_rep_index_in_a
  is_rep_index_in_a = ivector(a->nal);
  rep_index_in_a = ivector(ns);
  found = 0;
  for (i = 0; i < a->nal; i++) {
    is_rep_index_in_a[i] = -1;
    for (j = 0; j < ns; j++) {
      if (strcmp(a->aname[i], seq_name[j]) == 0) {
        is_rep_index_in_a[i] = j;
        rep_index_in_a[j] = i;
        if (debug > 1) cout << "rep_index_in_a: " << rep_index_in_a[j] << endl;
        found++;
        break;
      }
    }
  }
  if (found != ns) {
    cout << "Error: cannot find all the representative names in subalign a"
         << endl;
    exit(0);
  }

  // 3. determine the added_gap_number array [0..a->alilen]
  int lena = a->alilen;
  int *added_gap_number = ivector(a->alilen);
  for (i = 0; i <= lena; i++) added_gap_number[i] = 0;
  int *added_gap_number_1 = ivector(a->alilen);
  for (i = 0; i < ns; i++) {
    get_added_gap_number_1(a->aseq[rep_index_in_a[i]],
                           b[i]->aseq[rep_index_in_b[i]], added_gap_number_1,
                           a->alilen);
    for (j = 0; j <= lena; j++) {
      if (added_gap_number[j] < added_gap_number_1[j]) {
        added_gap_number[j] = added_gap_number_1[j];
      }
    }
  }
  if (debug > 1)
    for (j = 0; j <= lena; j++) {
      cout << "added_gap_number: " << j << " " << added_gap_number[j] << endl;
    }
  delete[] added_gap_number_1;

  // 4. get a new sequence alignment
  subalign *newalign = new subalign();
  newalign->nal = a->nal;
  newalign->alilen = a->alilen;
  newalign->mnamelen = a->mnamelen;
  int lenb;
  // cout << "nal0: " << newalign->nal << endl;
  for (i = 0; i < ns; i++) {
    newalign->nal += (b[i]->nal - 1);
    // cout << b[i]->nal-1 << endl;
    if (newalign->mnamelen < b[i]->mnamelen)
      newalign->mnamelen = b[i]->mnamelen;
  }
  // cout << "newalign: " << endl;
  // cout << "nal: " << newalign->nal << endl;
  // cout << "mnamelen: " << newalign->mnamelen << endl;
  for (i = 0; i <= lena; i++) {
    newalign->alilen += added_gap_number[i];
  }
  // cout << "alilen: " << newalign->alilen << endl;

  newalign->aname = cmatrix(newalign->nal, newalign->mnamelen);
  newalign->aseq = cmatrix(newalign->nal, newalign->alilen);
  char **newname, **newseq;
  newname = newalign->aname;
  newseq = newalign->aseq;
  int ni = 0;
  int nsi = 0;
  int bi;
  int repindex = 0;
  for (i = 0; i < a->nal; i++) {
    strcpy(newname[ni], a->aname[i]);
    nsi = 0;
    for (j = 0; j <= lena; j++) {
      for (k = 0; k < added_gap_number[j]; k++) {
        newseq[ni][nsi] = '-';
        nsi++;
      }
      if (j != lena) {
        newseq[ni][nsi] = a->aseq[i][j];
        nsi++;
      }
    }
    newseq[ni][nsi] = '\0';
    repindex = ni;
    ni++;
    // cout << "ni: " << ni << endl;

    if (is_rep_index_in_a[i] == -1) continue;
    bi = is_rep_index_in_a[i];
    added_gap_number_1 = ivector(b[bi]->alilen);
    // get_added_gap_number_1(b[bi]->aseq[rep_index_in_b[bi]], a->aseq[i],
    // added_gap_number_1, b[bi]->alilen);
    get_added_gap_number_1(b[bi]->aseq[rep_index_in_b[bi]], newseq[repindex],
                           added_gap_number_1, b[bi]->alilen);
    lenb = b[bi]->alilen;
    for (j = 0; j < b[bi]->nal; j++) {
      if (strcmp(seq_name[bi], b[bi]->aname[j]) == 0) continue;
      strcpy(newname[ni], b[bi]->aname[j]);
      nsi = 0;
      for (k = 0; k <= lenb; k++) {
        for (l = 0; l < added_gap_number_1[k]; l++) {
          newseq[ni][nsi] = '-';
          nsi++;
        }
        if (k != lenb) {
          newseq[ni][nsi] = b[bi]->aseq[j][k];
          nsi++;
        }
      }
      newseq[ni][nsi] = '\0';
      ni++;
      // cout << "ni: " << ni << endl;
    }
    delete[] added_gap_number_1;
  }
  delete[] added_gap_number;
  // newalign->printali(80);
  return newalign;
}

void get_added_gap_number_1(char *seq1, char *seq2, int *added_gap_number_1,
                            int dim) {
  int i, j, k;

  if (debug > 1) cout << "seq1: " << seq1 << endl;
  if (debug > 1) cout << "seq2: " << seq2 << endl;
  char *s, *ss;
  int g1, g2;
  int i1 = 0;
  for (i = 0; i <= dim; i++) added_gap_number_1[i] = 0;
  s = seq1;
  ss = seq2;
  while (i1 <= dim) {
    g1 = 0;
    g2 = 0;
    for (s = s; *s == '-'; s++) {
      g1++;
      i1++;
    }
    for (ss = ss; *ss == '-'; ss++) {
      g2++;
    }
    if (g2 > g1) added_gap_number_1[i1] = g2 - g1;
    if (debug > 1)
      cout << *s << " " << *ss << " " << g1 << " " << g2 << " " << i1 << endl;
    i1++;
    if (i1 > dim) break;
    s++;
    ss++;
  }
  // cout << seq1 << endl;
  // cout << seq2 << endl;
  // for(i=0;i<=dim;i++) { cout << added_gap_number_1[i] << " "; } cout << endl;
}

char ssint2ss(int i) {
  if (i == 1) return 'H';
  if (i == 2) return 'E';
  if (i == 3) return 'C';
}
