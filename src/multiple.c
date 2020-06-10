#include "multiple.h"
#include "btree_template.h"
#include "constraint.h"
#include "param.h"
#include "progressiveAlignHMM.h"
#include "sequences.h"
#include "subalign.h"

static int debug_here = 11;

FILE *logfp;
static int arrayindex = 0;

multiple::multiple(char *filename) {
  int i, j, k;

  strcpy(inputfileName, filename);

  char logfile[500];
  strcpy(logfile, inputfileName);
  strcat(logfile, ".promals.logfile");
  logfp = fopen(logfile, "w");
  if (!logfp) {
    cout << "cannot open log file " << logfile << " for write" << endl;
    exit(1);
  }

  distance_cutoff_similar = 0.4;
  N_small = 20;

  if (relax_count > 0) N_small = relax_count;

  max_cluster_elem = 200;

  // 1. read the sequence file
  allseqs.readFasta(filename, 1, 1);

  // 2. filter out highly similar sequences by cdhit, store them in similaraln,
  // and reassign allseqs
  if (filter_similar) {
    // the following function do the following things
    // 1. run cdhit on all the sequences, obtain clstr file
    // 2. read clstr file and run mafft on individual clusters to get similaraln
    // 3. change the allseqs object, so that it only contain cluster
    // representatives
    do_filter_similar(allseqs, filename, cdhit_c_option, similaraln);
    print_time_diff("filter_similar");

    get_rep_names(allseqs, similaraln, repnames);
    // for(i=0;i<repnames.size();i++) {cout << i << " " << repnames[i] << endl;}
    // cout << "This is the place" << endl;
  }

  // 3. calculating the distance matrix
  // a much simpler and much faster way of getting distances
  fprintf(logfp, "Start calculating the distance matrix ...");
  fflush(logfp);
  print_section_info("Below get distance matrix");
  allseqs.get_kmer_array(6);
  allseqs.get_kmer_distance(6);
  fprintf(logfp, " Done.\n\n");
  fflush(logfp);
  print_time_diff("calculate_distance_matrix");
}

// use a kcenter approach to derive clusters of sequences, use mafft to align
// sequences with each cluster and build a tree for cluster representatives
void multiple::build_cluster() {
  max_iteration = 20;
  max_initiation = 500;

  fprintf(logfp, "Start aligning similar sequences ...");
  fflush(logfp);
  // 1. obtain kcenter clusters
  get_kcenter_clusters(allseqs.nseqs, max_group_number, max_cluster_elem,
                       max_iteration, max_initiation);
  print_time_diff("get_kcenter_clusters");

  // 2. make alignments within clusters and build tree for cluster
  // representatives
  cluster2tree_and_similarSet();
  fprintf(logfp, " Done.\n");
  // print_time_diff("cluster2tree_and_similarSet_by_mafft");

  // 3. assign nodes to preAligned vector
  alltree.obtainPreAligned(alltree.root);
  fprintf(logfp, "      Number of input sequences: %d\n", allseqs.nseqs);
  fprintf(logfp, "      Number of pre-aligned groups: %d\n",
          alltree.preAligned.size());
  fprintf(logfp, "\n");
  fflush(logfp);

  // 4. prepare for restricted consistency measure
  map_allseqs_pos_to_tnode();
  get_distance_matrix_for_preAligned(N_small);
}

// derive clusters based on the UPGMA tree of all the sequences after filtering
void multiple::build_tree() {
  // build k-mer dist based tree
  fprintf(logfp, "Start making a guide tree ...");
  fflush(logfp);
  print_section_info("Below build the tree");
  alltree.UPGMA(allseqs.distMat, allseqs.seq, allseqs.name, allseqs.nseqs);
  print_time_diff("build_tree");
  fprintf(logfp, " Done.\n\n");
  fflush(logfp);
  // cout << "finished UPGMA" << endl;
}

void multiple::get_kcenter_clusters(int dim, int nc, int maxelem, int maxiter,
                                    int maxinit) {
  int i, j, k;

  // 1. build a similarity matrix
  int **simmat = imatrix(allseqs.nseqs, allseqs.nseqs);
  for (i = 1; i <= allseqs.nseqs; i++) {
    for (j = 1; j <= allseqs.nseqs; j++) {
      simmat[i][j] = int((1 - allseqs.distMat[i][j]) * 100);
    }
  }

  // 2. get kcenter clusters
  clusterindex = kcenter(simmat, dim, nc, maxelem, maxiter, maxinit);
  // Debug here, each sequence has a cluster index number, which range is
  // [1..nc]
  if (debug_here > 11) cout << "cluster indexes:" << endl;
  if (debug_here > 11)
    for (i = 1; i <= allseqs.nseqs; i++) {
      cout << i << " " << clusterindex[i] << endl;
    }

  // 3. clear up the similarity matrix
  free_imatrix(simmat, allseqs.nseqs, allseqs.nseqs);
}

void multiple::cluster2tree_and_similarSet() {
  int i, j, k;

  // 1. from clusterindex, get the mafft alignments for clusters
  int *index = ivector(max_cluster_elem);
  int clustersize;
  char tmpdir[200];
  char command[200];
  long int a1 = time(NULL);
  srand(a1);
  
  print_section_info("Below align each cluster");
  for (i = 1; i <= max_group_number; i++) {
    // 1.1 get the indexes in allseqs for each cluster
    clustersize = 0;
    for (j = 1; j <= allseqs.nseqs; j++) {
      if (clusterindex[j] == i) {
        clustersize++;
        index[clustersize] = j;
      }
    }
    // 1.2 write the sequences to a temporary file
    // 1.2.1 set up a directory
    int myrand = rand();
    
    sprintf(tmpdir, "/tmp/%d_%d", getpid(), myrand);
    sprintf(command, "mkdir %s", tmpdir);
    system(command);
    
    char tmpfa[300];
    sprintf(tmpfa, "%s/tmp.fa", tmpdir);
    ofstream ofp(tmpfa, ios::out);
    for (j = 1; j <= clustersize; j++) {
      ofp << ">" << allseqs.name[index[j]] << endl;
      ofp << allseqs.seq[index[j]] << endl;
    }
    ofp.close();
    ofp.clear();

    // 1.3 run mafft on the sequence
    if (clustersize > 1)
      sprintf(command,
              "%s --localpair --maxiterate 10 %s/tmp.fa 1>%s/tmp.aln.fa "
              "2> %s/tmp.err",
              mafft, tmpdir, tmpdir, tmpdir);
    else
      sprintf(command, "cp %s/tmp.fa %s/tmp.aln.fa", tmpdir, tmpdir);
    system(command);
    char outalnfile[200];
    sprintf(outalnfile, "%s/tmp.aln.fa", tmpdir);
    subalign *a = new subalign(outalnfile, "fasta", 25);
    prealn.push_back(a);
    
    if (debug_here > -11) cout << "mafft alignment for cluster " << i << endl;
    if (debug_here > -11) a->printali(80);
    
    // 1.4. clear up the tmp directory
    sprintf(command, "rm -rf %s", tmpdir);
    system(command);
  }
  print_time_diff("running mafft to get prealn");

  // 2. select representatives, and get alignments with single sequences
  // obtain a sequences object of representatives
  print_section_info("Below select representatives");
  sequences repseqs;
  repseqs.nseqs = max_group_number;
  repseqs.seq.push_back(string(""));
  repseqs.name.push_back(string(""));
  for (j = 0; j < max_group_number; j++) {
    cout << "cluster number: " << j + 1 << endl;
    int tmp_index = prealn[j]->select_representative_henikoff(0.5);
    char *tmp_name = new char[strlen(prealn[j]->aname[tmp_index]) + 1];
    strcpy(tmp_name, prealn[j]->aname[tmp_index]);
    int count_aa = 0;
    for (i = 1; i <= prealn[j]->alilen; i++) {
      if (prealn[j]->aseq[tmp_index][i] != '-') count_aa++;
    }
    char *tmp_seq = new char[count_aa + 1];
    int tmp_array_index = 0;
    for (i = 0; i < prealn[j]->alilen; i++) {
      if (prealn[j]->aseq[tmp_index][i] != '-') {
        tmp_seq[tmp_array_index] = prealn[j]->aseq[tmp_index][i];
        tmp_array_index++;
      }
    }
    tmp_seq[tmp_array_index] = '\0';
    repseqs.seq.push_back(string(tmp_seq));
    repseqs.name.push_back(string(tmp_name));
    delete[] tmp_name;
    delete[] tmp_seq;
  }

  // 3 build the tree based on repseqs
  repseqs.get_kmer_array(6);
  repseqs.get_kmer_distance(6);
  alltree.UPGMA(repseqs.distMat, repseqs.seq, repseqs.name, repseqs.nseqs);

  // 4. assign prealn to similarset
  for (i = 1; i <= max_group_number; i++) {
    alltree.v[i]->similarSet = prealn[i - 1];
    alltree.v[i]->aligned = 1;
  }
  print_time_diff("select_rep, tree_rep, assign similarset");
}

void multiple::set_distance_cutoff_similar(double dist_cutoff) {
  distance_cutoff_similar = dist_cutoff;
}

// for a large number of sequences, find the cutoff that results in a fixed
// number of pre-aligned groups
void multiple::set_distance_cutoff_similar(double dist_cutoff, int Ngroup) {
  int i;
  if (allseqs.nseqs <= Ngroup) {
    distance_cutoff_similar = dist_cutoff;
    return;
  }

  // find the distance to the leaf node for each node
  double dist2leaf[allseqs.nseqs];
  arrayindex = 1;
  distance2leaf(alltree.root, dist2leaf);
  sort(allseqs.nseqs - 1, dist2leaf);
  double tmp_dist_cutoff = (dist2leaf[allseqs.nseqs - 1 - (Ngroup - 1)] +
                            dist2leaf[allseqs.nseqs - 1 - (Ngroup - 2)]) /
                           2 * 2;
  if (tmp_dist_cutoff > dist_cutoff) {
    distance_cutoff_similar = tmp_dist_cutoff;
    fprintf(logfp,
            "The original identity threshold %3.2f results in the number of "
            "pre-aligned groups > %d.\n",
            1 - dist_cutoff, Ngroup);
    fprintf(logfp,
            "The new identity threshold is adjusted to %3.2f to set the number "
            "of pre-aligned groups to %d.\n\n",
            1 - tmp_dist_cutoff, Ngroup);
    fprintf(stdout,
            "  identity threshold adjusted to: %3.2f (to save time)\n\n",
            1 - distance_cutoff_similar);
  } else
    distance_cutoff_similar = dist_cutoff;
}

void multiple::distance2leaf(tnode *r, double *array) {
  double dist;
  if (r->childL == NULL) {
    return;
  }
  dist = 0;
  tnode *tmp = r;
  while (tmp->childL != NULL) {
    dist += tmp->childL->branchlen;
    tmp = tmp->childL;
  }
  array[arrayindex] = dist;
  arrayindex++;
  distance2leaf(r->childL, array);
  distance2leaf(r->childR, array);
}

void multiple::alignSimilar() {
  int i, j;

  // 1. progressively align similar sequences using general substitution matrix
  print_section_info("Below align similar sequences");
  if (allseqs.nseqs != 1) {
    fprintf(logfp, "Start aligning similar sequences ...");
    fflush(logfp);
    alltree.progressiveAlignHMM_FastStage_mafft(alltree.root,
                                                distance_cutoff_similar / 2);
    fprintf(logfp, " Done.\n");
    fflush(logfp);
  } else {
    alltree.root->aligned = 1;
    alltree.root->aln->printali(80);
  }
  
  // 2. store pre-aligned groups in the stopped nodes and select one
  // representative from each group
  print_section_info("Below select representatives");
  store_similar_henikoff(alltree.root);

  // 3. obtain nodes of for preAligned vector
  alltree.obtainPreAligned(alltree.root);
  fprintf(logfp, "      Number of input sequences: %d\n", allseqs.nseqs);
  fprintf(logfp, "      Number of pre-aligned groups: %d\n",
          alltree.preAligned.size());
  fprintf(logfp, "\n");
  fflush(logfp);

  // 4. prepare for restricted consistency
  map_allseqs_pos_to_tnode();
  get_distance_matrix_for_preAligned(N_small);

  // output distance matrix for representatives
  if (debug > 1) {
    for (i = 0; i < alltree.preAligned.size(); i++) {
      for (j = 0; j < alltree.preAligned.size(); j++) {
        cout << dist_matrix_preAligned[i][j] << " ";
      }
      cout << endl;
    }
  }
  if (debug_here > 1) {
    cout << "NUMBER OF SEQUENCES: " << allseqs.nseqs
         << " NUMBER OF GROUPS: " << alltree.preAligned.size() << endl;
  }

  print_time_diff("align_similar_and_select_reps");
}

extern hmm_psipred_parameters *params1;

void multiple::alignDivergent_psipred(int use_homologs) {
  int i, j;

  if (alltree.preAligned.size() == 1) {
    fprintf(logfp, "Start getting structural information ...\n");
  } else
    fprintf(logfp, "Start aligning divergent groups ...\n");
  fflush(logfp);
  // right now, just the option of multim - for probablistic consistency
  hmm_psipred_parameters params(psipred_env_number);
  params.read_parameters(psipred_parameter_file, psipred_env_number, 1);
  params1 = &params;

  // select representatives
  for (i = 0; i < alltree.preAligned.size(); i++) {
    alltree.preAligned[i]->aln->get_oneSeqAln(0);
  }

  print_section_info("Below run or read psiblast alignment and psipred");

  // use_homologs==1: use profile of pre-aligned group
  if (use_homologs == 1) {
    for (i = 0; i < alltree.preAligned.size(); i++)
      alltree.preAligned[i]->aux_align =
          alltree.preAligned[i]->similarSet->purge_align_one_seq_name(
              alltree.preAligned[i]->aln->aname[0]);
  }
  // use_homologs==2: use database homologs
  else if (use_homologs == 2) {
    fprintf(logfp, "      - Running PSI-BLAST and PSIPRED ...\n");
    fprintf(logfp, "        ");
    fflush(logfp);
    for (i = 0; i < alltree.preAligned.size(); i++) {
      if (clean_blast_before)
        clean_blast_psipred(alltree.preAligned[i]->aln->aname[0]);
      alltree.preAligned[i]->aux_align =
          get_blastpgp_alignment(alltree.preAligned[i]->aln->aname[0],
                                 alltree.preAligned[i]->aln->aname[0],
                                 alltree.preAligned[i]->aln->aseq[0]);
      alltree.preAligned[i]->aln->get_ss_prof1(
          blast_dir, alltree.preAligned[i]->aln->aname[0], runpsipred1_command);
      if (!alltree.preAligned[i]->aux_align) {
        cout << "reading blastpgp alignment failed "
             << alltree.preAligned[i]->aln->aname[0] << endl;
        alltree.preAligned[i]->aux_align = alltree.preAligned[i]->aln;
      }
      if (clean_blast_after)
        clean_blast_psipred(alltree.preAligned[i]->aln->aname[0]);
      fprintf(logfp, "*");
      fflush(logfp);
    }
    fprintf(logfp, "\n");
  }
  // use_homologs==0: use single sequence as profile
  else {
    for (i = 0; i < alltree.preAligned.size(); i++)
      alltree.preAligned[i]->aux_align = alltree.preAligned[i]->aln;
  }

  print_time_diff("psiblast_and_psipred");
  print_section_info("Below calculate profiles");
  subalign *taln;
  fprintf(logfp, "      - Calculating profiles ...\n");
  fprintf(logfp, "        ");
  fflush(logfp);
  for (i = 0; i < alltree.preAligned.size(); i++) {
    alltree.preAligned[i]->aux_align->ss = alltree.preAligned[i]->aln->ss;
    taln = alltree.preAligned[i]->aux_align;
    taln->prof(1);
    // prof_freq allocation
    taln->prof_freq = dmatrix(taln->prof_len, 20);
    // use_ss_freq:
    // 0 - use blosum62 matrix to get pseudocount
    // 1 - use aa_pair1[0] and aa_bg1[0] to get pseudocount
    // 2 - use aa_pair1 and aa_bg1 and ss->alphabet1 to get pseudocount
    if (use_ss_freq == 0) {
      taln->pseudoCounts(taln->prof_effn, taln->prof_nef, taln->prof_len,
                         taln->prof_freq);
    } else if (use_ss_freq == 1) {
      taln->pseudoCounts(taln->prof_effn, taln->prof_nef, taln->prof_len,
                         taln->prof_freq, params.aa_pair1[0], params.aa_bg1[0]);
    } else if (use_ss_freq == 2) {
      if (psipred_env_number == 9)
        taln->pseudoCounts(taln->prof_effn, taln->prof_nef, taln->prof_len,
                           taln->prof_freq, params.aa_pair1, params.aa_bg1,
                           taln->ss->alphabet1);
      if (psipred_env_number == 3)
        taln->pseudoCounts(taln->prof_effn, taln->prof_nef, taln->prof_len,
                           taln->prof_freq, params.aa_pair1, params.aa_bg1,
                           taln->ss->sstype);
    }
    taln->log_pseudoCounts();
    
    fprintf(logfp, "*");
    fflush(logfp);
  }
  fprintf(logfp, "\n");

  print_time_diff("calculate_profiles");
  if (debug_here > 11) {
    cout << "Before consistency" << endl;
  }
  fprintf(logfp, "      - Making consistency scoring fuction ...\n");
  fprintf(logfp, "        ");
  fflush(logfp);
  
  if (struct_weight == 0) {
    print_section_info("Below make consistency scoring");
    alltree.profileConsistency_psipred(&params, dist_matrix_preAligned, 1,
                                       use_homologs);
    alltree.relaxConsistMatrix(dist_matrix_preAligned, 1, minProb);
    alltree.relaxConsistMatrix(dist_matrix_preAligned, 1, minProb * 0.1);
  } else if (!use_updated_database) {
    print_section_info("Below get seq_str_aln, old database");
    ssaln = *(get_seq_str_alns());

    print_section_info("Below make consistency scoring");
    alltree.profileConsistency_psipred(&params, dist_matrix_preAligned, 1,
                                       use_homologs);
    print_time_diff("profile_profile_alignments");
    if (before_relax_combine == 0)
      combine_structure_alignments3(alltree, ssaln);
    alltree.relaxConsistMatrix(dist_matrix_preAligned, 1, minProb);
    print_time_diff("relax_consistency_first_round");
    if (before_relax_combine == 1)
      combine_structure_alignments3(alltree, ssaln);
    alltree.relaxConsistMatrix(dist_matrix_preAligned, 1, minProb * 0.1);
    print_time_diff("relax_consistency_second_round");
    if (before_relax_combine == 2)
      combine_structure_alignments3(alltree, ssaln);
  } else {
    print_section_info("Below get seq_str_aln, updated database");
    ssaln = *(get_seq_str_alns1());

    print_section_info("Below make consistency scoring");
    alltree.profileConsistency_psipred(&params, dist_matrix_preAligned, 1,
                                       use_homologs);
    print_time_diff("profile_profile_alignments");
    if (before_relax_combine == 0)
      combine_structure_alignments3_updated(alltree, ssaln);
    alltree.relaxConsistMatrix(dist_matrix_preAligned, 1, minProb);
    print_time_diff("relax_consistency_first_round");
    if (before_relax_combine == 1)
      combine_structure_alignments3_updated(alltree, ssaln);
    alltree.relaxConsistMatrix(dist_matrix_preAligned, 1, minProb * 0.1);
    print_time_diff("relax_consistency_second_round");
    if (before_relax_combine == 2)
      combine_structure_alignments3_updated(alltree, ssaln);
  }
  print_time_diff("structure comparisons and constraints");

  // read and combine structure constraints
  if (strlen(constraint_file) != 0) {
    print_section_info("Below combine outside structural constraints");
    vector<constraint *> vcons = read_multiple_constraint(constraint_file);
    for (i = 0; i < vcons.size(); i++) {
      if (vcons[i]->nseqs == 0) continue;
      cout << "constraint: " << i << endl;
      vcons[i]->assign_prealigned(&alltree.preAligned);
      vcons[i]->assign_originalseq(&allseqs);
      vcons[i]->printSeqs();
      vcons[i]->allocate_index();
      vcons[i]->checkNamesSequences();
      combine_constraint(alltree, *(vcons[i]), pdb_weight);
    }
  }
  // read and combine user-defined constraints
  if (strlen(user_constraint) != 0) {
    print_section_info("Below combine user constraints");
    vector<constraint *> vcons1 = read_multiple_constraint(user_constraint);
    for (i = 0; i < vcons1.size(); i++) {
      if (vcons1[i]->nseqs == 0) continue;
      cout << "constraint: " << i << endl;
      vcons1[i]->assign_prealigned(&alltree.preAligned);
      vcons1[i]->assign_originalseq(&allseqs);
      vcons1[i]->printSeqs();
      vcons1[i]->allocate_index();
      vcons1[i]->checkNamesSequences();
      combine_constraint(alltree, *(vcons1[i]), user_constraint_weight);
    }
  }

  if (debug_here > 11) {
    cout << "After consistency" << endl;
  }

  print_section_info("Below compute consistency alignment");
  fprintf(logfp, "      - Making progressive alignments ...\n");
  fflush(logfp);
  alltree.computeConsistencyAlignment(alltree.root);
  print_time_diff("compute_consistency_alignment");

  int iterRef_round = 0;
  alltree.iterativeRefinement(iterRef_round);
  cout << "iterative refinement for " << iterRef_round << " rounds." << endl;
  print_time_diff("iterative_refinement");

  print_section_info("Below refine and print alignment");
  output_alignment();
  fprintf(logfp, "\nPROMALS is now finished\n");
  fflush(logfp);
  fclose(logfp);

}

void multiple::output_alignment() {
  int i;
  char outFileName[200];
  strcpy(outFileName, inputfileName);
  int tmpLen = strlen(outFileName);
  if ((outFileName[tmpLen - 1] == 'a') && (outFileName[tmpLen - 2] == 'f') &&
      (outFileName[tmpLen - 3] == '.')) {
    outFileName[tmpLen - 3] = '\0';
  }
  strcat(outFileName, ".promals.aln");
  if (!outFile.empty()) {
    strcpy(outFileName, outFile.c_str());
  }
  alltree.printAlignmentFromAbs(alltree.root, outFileName, similaraln,
                                repnames);
  cout << "  program finished" << endl << endl;

  cout << "original sequence names:" << endl;
  for (i = 1; i <= allseqs.nseqs; i++) {
    cout << allseqs.name[i] << endl;
  }
}

// store pre-aligned groups in the stopped nodes and select one representative
// from each group
void multiple::store_similar_henikoff(tnode *r) {
  int i, j;

  // cout << "here" << endl;
  if (!r->aligned) {
    store_similar_henikoff(r->childL);
    store_similar_henikoff(r->childR);
    return;
  }

  if (r->aln->nal == 1) return;

  // store the original subalign to "similarSet"
  // r->similarSet = new subalign(*(r->aln));
  if (debug_here > 11) r->similarSet->printali(80);

  // find a representative sequence for the group, make it the "aln",
  // modify the representative sequence if necessary (adding small letters)
  int tmp_index = r->aln->select_representative_henikoff(0.5);
  r->similarSet =
      new subalign(*(r->aln));  // after modification of representative, store
                                // similar, instead of store similar before
  char *tmp_name = new char[strlen(r->aln->aname[tmp_index]) + 1];
  int count_aa = 0;
  for (i = 1; i <= r->aln->alilen; i++) {
    if (r->aln->aseq[tmp_index][i] != '-') count_aa++;
  }
  char *tmp_seq = new char[count_aa + 1];
  strcpy(tmp_name, r->aln->aname[tmp_index]);
  int tmp_array_index = 0;
  for (i = 0; i < r->aln->alilen; i++) {
    if (r->aln->aseq[tmp_index][i] != '-') {
      tmp_seq[tmp_array_index] = r->aln->aseq[tmp_index][i];
      tmp_array_index++;
    }
  }
  tmp_seq[tmp_array_index] = '\0';
  r->aln = oneSeq2subalign(tmp_seq, tmp_name);
  if (debug_here > 11) r->aln->printali(60);
}

// determine the position in the allseqs for any tnode in preAligned vector
//                p_seq
void multiple::map_allseqs_pos_to_tnode() {
  int i, j;
  subalign *tmpaln;

  if (debug_here > 11)
    cout << alltree.preAligned.size() << "  " << allseqs.nseqs << endl;

  for (i = 0; i < alltree.preAligned.size(); i++) {
    tmpaln = alltree.preAligned[i]->aln;
    for (j = 1; j <= allseqs.nseqs; j++) {
      if (strcmp(tmpaln->aname[0], allseqs.name[j].c_str()) == 0) {
        alltree.preAligned[i]->p_seq = j;
        if (debug_here > 11)
          cout << i << " " << j << " " << tmpaln->aname[0] << endl;
        continue;
      }
    }
  }
}

// find N sequences with smallest distances to a sequence in the preAligned
// vector store them in a distance matrix
void multiple::get_distance_matrix_for_preAligned(int N_smallest) {
  int i, j;

  double tmp_dist_array[alltree.preAligned.size() + 1];
  int auxilary[alltree.preAligned.size() + 1];

  dist_matrix_preAligned =
      dmatrix(alltree.preAligned.size(), alltree.preAligned.size());

  for (i = 0; i < alltree.preAligned.size(); i++) {
    for (j = 0; j < alltree.preAligned.size(); j++) {
      tmp_dist_array[j + 1] = 0;
      tmp_dist_array[j + 1] = allseqs.distMat[alltree.preAligned[i]->p_seq]
                                             [alltree.preAligned[j]->p_seq];
      auxilary[j + 1] = j;
    }

    sort2(alltree.preAligned.size(), tmp_dist_array, auxilary);

    cout << "Neighbors of " << i << endl;
    for (j = 1; j <= ((N_smallest + 1 < alltree.preAligned.size())
                          ? N_smallest + 1
                          : alltree.preAligned.size());
         j++) {
      dist_matrix_preAligned[i][auxilary[j]] = tmp_dist_array[j];
      cout << alltree.preAligned[i]->aln->aname[0] << "\t"
           << alltree.preAligned[auxilary[j]]->aln->aname[0] << endl;
    }
    cout << endl;
  }
  if (debug_here > 11)
    for (i = 0; i < alltree.preAligned.size(); i++) {
      for (j = 0; j < alltree.preAligned.size(); j++) {
        cout << "distance " << i << " " << j << " "
             << dist_matrix_preAligned[i][j] << endl;
      }
    }
}

vector<seq_str_aln *> *multiple::get_seq_str_alns() {
  int i, j, k;
  char database[500];
  char suffix[20];
  char chkfile[200];
  char options[100];
  char fastafile[500];
  strcpy(suffix, "pdb.br");
  strcpy(options, "-z 1665828471 -e 0.001 -h 0.001");
  strcpy(database, "/home/jpei/promals/src_structure/structure_db/dali.fa");

  k = strlen(blast_dir);
  if (blast_dir[k - 1] != '/') {
    assert(k < 500 - 1);
    blast_dir[k] = '/';
    blast_dir[k + 1] = '\0';
  }

  cout << "\tid_cutoff: " << struct_id_cutoff << endl;
  cout << "\tbelow_id_cutoff: " << below_id_cutoff << endl << endl;
  for (i = 0; i < (int)alltree.preAligned.size(); i++) {
    seq_str_aln *tmpss = new seq_str_aln(alltree.preAligned[i]->aln->aseq[0]);
    tmpss->set_id_cutoff(struct_id_cutoff);
    tmpss->set_below_id_cutoff(below_id_cutoff);
    sprintf(fastafile, "%s%s.fa", blast_dir,
            alltree.preAligned[i]->aln->aname[0]);
    sprintf(chkfile, "%s%s.chk", blast_dir,
            alltree.preAligned[i]->aln->aname[0]);
    tmpss->run_blast(blastpgp_command, fastafile, database, suffix, chkfile,
                     options);
    tmpss->read_blast_results();
    tmpss->get_prof();
    ssaln.push_back(tmpss);
  }
  return &ssaln;
}

// new for all the updated structures
vector<seq_str_aln *> *multiple::get_seq_str_alns1() {
  int i, j, k;
  char database[500];
  char suffix[20];
  char chkfile[200];
  char options[100];
  char fastafile[500];
  strcpy(suffix, "pdb.br");
  strcpy(options, "-z 1665828471 -e 0.001 -h 0.001");
  sprintf(database, "%s/db/structure_db/struct.fasta", program_dir);

  k = strlen(blast_dir);
  if (blast_dir[k - 1] != '/') {
    assert(k < 500 - 1);
    blast_dir[k] = '/';
    blast_dir[k + 1] = '\0';
  }

  cout << "\tid_cutoff: " << struct_id_cutoff << endl;
  cout << "\tbelow_id_cutoff: " << below_id_cutoff << endl << endl;
  for (i = 0; i < (int)alltree.preAligned.size(); i++) {
    seq_str_aln *tmpss = new seq_str_aln(alltree.preAligned[i]->aln->aseq[0]);
    tmpss->set_id_cutoff(struct_id_cutoff);
    tmpss->set_below_id_cutoff(below_id_cutoff);
    tmpss->set_subject_N_C_extension(5);
    sprintf(fastafile, "%s%s.fa", blast_dir,
            alltree.preAligned[i]->aln->aname[0]);
    sprintf(chkfile, "%s%s.chk", blast_dir,
            alltree.preAligned[i]->aln->aname[0]);
    tmpss->run_blast(blastpgp_command, fastafile, database, suffix, chkfile,
                     options);
    tmpss->read_blast_results();
    tmpss->get_prof_update();
    if (debug_here > 11) tmpss->print_result();
    ssaln.push_back(tmpss);
  }
  return &ssaln;
}

// retricted kcenter method
// simmat: simlarity matrix of dimXdim, nc: number of clusters;
// maxelem: maximum number of elements in a cluster
int *kcenter(int **simmat, int dim, int nc, int maxelem, int maxiter,
             int maxinit) {
  int i, j, k, iter, maxclusterindex, maxsim;
  long int sum_maxsumsim, lsum_maxsumsim, ilsum_maxsumsim;
  int maxsumsim, sumsim, tmpvalue, centerindex;
  int centercount, found;
  int init;

  assert(nc < dim);

  // allocate memory
  int **clusters = imatrix(nc, maxelem);
  int *clusterindex = ivector(dim);
  int *tmpclusterindex = ivector(dim);
  int *clustersize = ivector(nc);
  int iscenter;
  long int a1 = time(NULL);
  srand(a1);
  srand(1);
  int randindex;

  print_section_info("Below use kcenter approach to get clusters");
  cout << "\tnumber of clusters: " << nc << endl;
  cout << "\tnumber of elementes: " << dim << endl;
  cout << "\tmaximum number of elements in a cluster: " << maxelem << endl;
  cout << "\tmaximum number of initiations: " << maxinit << endl;
  cout << "\tmaximum number of iterations in each init: " << maxiter << endl;
  cout << endl;

  // Debug here
  // for(i=1;i<=dim;i++) { for(j=1;j<=dim;j++) { cout << simmat[i][j] << "\t"; }
  // cout << endl; }

  // for(random_1

  ilsum_maxsumsim = -10000;
  for (init = 1; init <= maxinit; init++) {
    if (debug_here > 11) cout << "===========================" << endl;
    if (debug_here > 11) cout << "initiation number " << init << endl;

    // initialize the cluster centers
    centercount = 0;
    while (centercount < nc) {
      while (dim <= (randindex = rand() / (RAND_MAX / dim)))
        ;
      randindex++;
      // cout << randindex << endl;
      found = 0;
      for (i = 1; i <= centercount; i++) {
        if (randindex == clusters[i][0]) {
          found = 1;
          break;
        }
      }
      if (found) continue;
      centercount++;
      clusters[centercount][0] = randindex;
    }

    lsum_maxsumsim = -10000;
    for (iter = 1; iter <= maxiter; iter++) {
      if (debug_here > 11) cout << "-----------" << endl;
      if (debug_here > 11) cout << "iteration number " << iter << endl;
      // assign elements to clusters
      for (i = 1; i <= nc; i++) clustersize[i] = 0;
      for (i = 1; i <= dim; i++) {
        maxclusterindex = -1;
        maxsim = -100000;
        for (j = 1; j <= nc; j++) {
          if (clustersize[j] == maxelem) continue;
          if (simmat[i][clusters[j][0]] > maxsim) {
            maxsim = simmat[i][clusters[j][0]];
            maxclusterindex = j;
          }
        }
        clustersize[maxclusterindex]++;
        clusters[maxclusterindex][clustersize[maxclusterindex]] = i;
      }

      // update the cluster centers
      sum_maxsumsim = 0;
      for (i = 1; i <= nc; i++) {
        maxsumsim = -1000000;
        sumsim = 0;
        for (j = 1; j <= clustersize[i]; j++) {
          sumsim = 0;
          for (k = 1; k <= clustersize[i]; k++) {
            sumsim += simmat[clusters[i][j]][clusters[i][k]];
          }
          if (sumsim > maxsumsim) {
            maxsumsim = sumsim;
            centerindex = j;
          }
        }
        // Debug here
        if (debug_here > 11) {
          int minsim = 100000000;
          for (j = 1; j <= clustersize[i]; j++) {
            if (minsim > simmat[clusters[i][centerindex]][clusters[i][j]])
              minsim = simmat[clusters[i][centerindex]][clusters[i][j]];
          }
          cout << "cluster: " << i << " size: " << clustersize[i]
               << " average maxsumsim: " << 1.0 * maxsumsim / clustersize[i]
               << " min: " << minsim << endl;
        }
        sum_maxsumsim += maxsumsim;
        clusters[i][0] = clusters[i][centerindex];
      }

      // print information
      if (debug_here > 11)
        cout << "after round " << iter << " : " << sum_maxsumsim << endl;

      // stop criteria
      if (lsum_maxsumsim < sum_maxsumsim) {
        lsum_maxsumsim = sum_maxsumsim;
        // update the result
        for (i = 1; i <= nc; i++) {
          for (j = 1; j <= clustersize[i]; j++) {
            tmpclusterindex[clusters[i][j]] = i;
          }
        }
      }
    }  // end of the iterations for one initial condition

    if (debug_here > 11) cout << "lsum_maxsumsim: " << lsum_maxsumsim << endl;

    // check if the previous condition gives a better result
    if (lsum_maxsumsim > ilsum_maxsumsim) {
      for (i = 1; i <= dim; i++) clusterindex[i] = tmpclusterindex[i];
      ilsum_maxsumsim = lsum_maxsumsim;
    }
    if (debug_here > 11) cout << "ilsum_maxsumsim: " << ilsum_maxsumsim << endl;

  }  // end of different initiations

  // print out the average identity and min identity of each cluster
  // recreate the clusters
  for (i = 1; i <= nc; i++) clustersize[i] = 0;
  for (i = 1; i <= dim; i++) {
    clustersize[clusterindex[i]] += 1;
    clusters[clusterindex[i]][clustersize[clusterindex[i]]] = i;
  }
  for (i = 1; i <= nc; i++) {
    maxsumsim = -1000000;
    sumsim = 0;
    for (j = 1; j <= clustersize[i]; j++) {
      sumsim = 0;
      for (k = 1; k <= clustersize[i]; k++) {
        sumsim += simmat[clusters[i][j]][clusters[i][k]];
      }
      if (sumsim > maxsumsim) {
        maxsumsim = sumsim;
        centerindex = j;
      }
    }
    int minsim = 100000000;
    for (j = 1; j <= clustersize[i]; j++) {
      if (minsim > simmat[clusters[i][centerindex]][clusters[i][j]])
        minsim = simmat[clusters[i][centerindex]][clusters[i][j]];
    }
    cout << "\tcluster: " << i << " size: " << clustersize[i]
         << " average maxsumsim: " << 1.0 * maxsumsim / clustersize[i]
         << " minsim: " << minsim << endl;
    sum_maxsumsim += maxsumsim;
    clusters[i][0] = clusters[i][centerindex];
  }
  cout << "\tfinal ilsum_maxsumsim: " << ilsum_maxsumsim << endl;

  delete[] tmpclusterindex;
  delete[] clustersize;
  free_imatrix(clusters, nc, maxelem);
  return clusterindex;
}

