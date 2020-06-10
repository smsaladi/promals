#include "all.h"

// int id_cutoff = 30;

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

seq_str_aln::seq_str_aln() {
  query = NULL;
  len = 0;
  id_cutoff = 30;
  nhits = 0;
  subject_N_C_extension = 5;  // default
}

seq_str_aln::seq_str_aln(char *q) {
  len = strlen(q);
  query = new char[len + 1];
  strcpy(query, q);
  aln.push_back(NULL);
  id.push_back(NULL);
  start.push_back(-1);
  end.push_back(-1);
  str_start.push_back(-1);
  str_end.push_back(-1);
  slen.push_back(-1);
  id_cutoff = 30;
  nhits = 0;
  subject_N_C_extension = 5;  // default
}

void seq_str_aln::set_id_cutoff(int a) { id_cutoff = a; }

void seq_str_aln::set_below_id_cutoff(int a) { below_id_cutoff = a; }

seq_str_aln::~seq_str_aln() {
  int i;
  if (query) {
    delete[] query;
  }
  for (i = 0; i < aln.size(); i++) {
    if (aln[i]) delete[] aln[i];
    if (id[i]) delete[] id[i];
  }
  aln.clear();
  id.clear();
  start.clear();
  end.clear();
}

void seq_str_aln::run_blast(char *blastexe, char *fasta_file, char *database,
                            char *suffix, char *chkfile, char *options) {
  int j;
  char blast_cmd[600];
  // sprintf(blast_cmd, "ls %s", fasta_file);
  // system(blast_cmd);
  // check if the checkpoint file exists
  char tmpcmd[600];
  sprintf(tmpcmd, "ls %s", chkfile);
  int a1 = system(tmpcmd);
  if (a1 == 0) {
    sprintf(blast_cmd, "%s -i %s -d %s -o %s.%s -R %s %s 2>/dev/null", blastexe,
            fasta_file, database, fasta_file, suffix, chkfile, options);
  } else {
    char alnfile[600];
    strcpy(alnfile, fasta_file);
    for (j = strlen(alnfile); j >= 0; j--)
      if (alnfile[j] == '.') {
        alnfile[j] = '\0';
        break;
      }
    strcat(alnfile, ".aln");
    sprintf(blast_cmd, "%s -i %s -d %s -o %s.%s -B %s -C %s %s 2>/dev/null",
            blastexe, fasta_file, database, fasta_file, suffix, alnfile,
            chkfile, options);
  }
  sprintf(blastoutput, "%s.%s", fasta_file, suffix);
  cout << blast_cmd << endl;
  // cout << "print working directory" << endl;
  // system("pwd");
  system(blast_cmd);
}

/* format of blast output
>6827
          Length = 388

 Score = 29.3 bits (64), Expect = 0.053,   Method: Composition-based stats.
 Identities = 14/30 (46%), Positives = 16/30 (53%)

Query: 4  YSGLAYTYHLQGNFSAAISYYHKALWLKPD 33
          YS L   Y  +G    AI +Y  AL LKPD
Sbjct: 70 YSNLGNVYKERGQLQEAIEHYRHALRLKPD 99



 Score = 21.6 bits (44), Expect =    11,   Method: Composition-based stats.
 Identities = 11/31 (35%), Positives = 16/31 (51%)

Query: 3   TYSGLAYTYHLQGNFSAAISYYHKALWLKPD 33
           + + LA     QGN   A+  Y KAL + P+
Sbjct: 307 SLNNLANIKREQGNIEEAVRLYRKALEVFPE 337

>0619
          Length = 421

 Score = 15.0 bits (27), Expect =   999,   Method: Composition-based stats.
 Identities = 5/17 (29%), Positives = 8/17 (47%)

Query: 16  NFSAAISYYHKALWLKP 32
           NF    ++     W+KP
Sbjct: 299 NFDINRNFLFAGDWMKP 315


  Database: dali.fa
    Posted date:  Apr 26, 2007  5:34 PM


*/

void seq_str_aln::read_blast_results() {
  int i, j, k;
  char line[200];
  char qr[2000];
  char sj[2000];
  char strfrag[100];
  int qstart, qend;
  int sstart, send;
  int tmpint;
  int identity;
  char name[100];
  int newlinenum = 0;
  int sbjctflag = 0;
  int firstblockflag = 0;
  int qp, sp;
  int subjectlength;
  char tmps[10], tmps1[10];

  char *s, *ss;

  // get ids and ranges
  FILE *brfile = fopen(blastoutput, "r");
  if (brfile == NULL) {
    cout << "Error: cannot open blastoutput file in seq_str_aln" << blastoutput
         << endl;
    exit(1);
  }

  cout << "blastoutput file: " << blastoutput << endl;
  // determine query_name
  for (i = strlen(blastoutput) - 1; i >= 0; i--) {
    if (blastoutput[i] == '/') {
      break;
    }
  }
  strcpy(query_name, blastoutput + i + 1);
  // cout << "query_name: " << query_name << endl;

  while (fgets(line, 200, brfile)) {
    if (line[0] == '>') {
      strcpy(name, &line[1]);
      // cout << name << endl;
      for (i = 0; i < strlen(name); i++) {
        if ((name[i] == ' ') || (name[i] == '\n')) name[i] = '\0';
      }
      sbjctflag = 0;
      continue;
    }
    if (strncmp(line, "          Length = ", 19) == 0) {
      sscanf(line, "%s%s%d", tmps, tmps1, &subjectlength);
      // slen.push_back(subjectlength);
      // cout << "length " << subjectlength << endl;
    }
    if (strncmp(line, " Identities = ", 14) == 0) {
      strcpy(qr, "");
      strcpy(sj, "");
      for (s = line; *s != '('; s++) {
        ;
      }
      sscanf(s + 1, "%d", &identity);
      // cout << identity << endl;
      firstblockflag = 1;
      newlinenum = 0;
      continue;
    }
    if (strncmp(line, "Query:", 6) == 0) {
      for (s = line; *s != ' '; s++) {
        ;
      }
      for (ss = s; *ss == ' '; ss++) {
        ;
      }
      // for(s=ss;*s!=' ';s++) {;}
      // for(ss=s;*ss==' ';ss++) {;}
      if (firstblockflag)
        sscanf(ss, "%d%s%d", &qstart, strfrag, &qend);
      else
        sscanf(ss, "%d%s%d", &tmpint, strfrag, &qend);
      strcat(qr, strfrag);
      // cout << "qstart and qr: " << qstart << " " << qr << " " << qend <<
      // endl;
      sbjctflag = 0;
    }
    if (strncmp(line, "Sbjct:", 6) == 0) {
      for (s = line; *s != ' '; s++) {
        ;
      }
      for (ss = s; *ss == ' '; ss++) {
        ;
      }
      // for(s=ss;*s!=' ';s++) {;}
      // for(ss=s;*ss==' ';ss++) {;}
      if (firstblockflag)
        sscanf(ss, "%d %s %d", &sstart, strfrag, &send);
      else
        sscanf(ss, "%d %s %d", &tmpint, strfrag, &send);
      strcat(sj, strfrag);
      // cout << "sstart and sj: " << sstart << " " << sj << " " << send <<
      // endl;
      sbjctflag = 1;
      newlinenum = 0;
      firstblockflag = 0;
    }
    if (strcmp(line, "\n") == 0) {
      newlinenum++;
    }
    if (sbjctflag && (newlinenum == 2)) {
      // cout << "Identity: " << identity << endl; //cout << id_cutoff << endl;
      // //cout << below_id_cutoff << endl; //cout << qstart << endl; //cout <<
      // qend << endl;
      if (identity < id_cutoff) continue;
      if (identity > below_id_cutoff) continue;
      // if(checkboundary(qstart, qend) == 1) continue;
      if ((checkboundary(qstart, qend) == 1) && (nhits >= max_struct_number))
        continue;

      start.push_back(qstart);
      end.push_back(qend);
      /* old version
      str_start.push_back(sstart);
      str_end.push_back(send);
      */
      // new
      str_start.push_back(MAX(1, sstart - subject_N_C_extension));
      str_end.push_back(MIN(subjectlength, send + subject_N_C_extension));
      slen.push_back(subjectlength);
      char *tmpName = new char[strlen(name) + 1];
      strcpy(tmpName, name);
      id.push_back(tmpName);
      // cout << "name: " << name << endl;
      // cout << "subjectlength: " << subjectlength << endl;

      int *tmpintseq = new int[subjectlength + 1];
      for (i = 0; i <= subjectlength; i++) tmpintseq[i] = 0;
      aln.push_back(tmpintseq);
      qp = qstart;
      sp = sstart;
      // cout << qstart << "  " << sstart << endl;
      // cout << "qr: " << qr << endl;
      // cout << "sj: " << sj << endl;
      for (i = 0; i < strlen(qr); i++) {
        if ((qr[i] != '-') && (sj[i] != '-')) {
          tmpintseq[sp] = qp;
        }
        if (qr[i] != '-') qp++;
        if (sj[i] != '-') sp++;
      }
      nhits++;
      cout << id[nhits];
      cout << " str_start: " << str_start[nhits]
           << " str_end: " << str_end[nhits] << endl;
      fprintf(stdout, "%-5d %s %-5d\n", qstart, qr, qend);
      fprintf(stdout, "%-5d %s %-5d\n", sstart, sj, send);

      // Debug here
      // for(i=1;i<=subjectlength;i++) cout << "aln: " << i << " " <<
      // aln[nhits][i] << endl;
    }
  }
  cout << endl;
  fclose(brfile);
}

int seq_str_aln::checkboundary(int qstart, int qend) {
  int i, j;
  int checkresultflag = 0;
  int extendedlen;
  for (i = 1; i <= nhits; i++) {
    extendedlen = 0;
    if (qstart > end[i]) continue;
    if (qend < start[i]) continue;
    if (start[i] > qstart) extendedlen = start[i] - qstart;
    if (end[i] < qend) extendedlen -= (end[i] - qend);
    if (extendedlen < 30) {
      checkresultflag = 1;
      break;
    }
  }
  return checkresultflag;
}

void seq_str_aln::get_prof() {
  int i, j;
  char domain_sequence[5000];
  char domain_sequence_name[200];
  char blast_dir[300];
  char blast_aln_file[300];
  subalign *tmpaln, *taln;
  hmm_psipred_parameters params(psipred_env_number);
  params.read_parameters(psipred_parameter_file, psipred_env_number, 1);

  // the directory containing all the SCOP domain sequences in fasta format
  // each file contains a single domain sequence
  strcpy(domain_sequence_name,
         "/home/jpei/promals/src_structure/structure_db/promals/blastresult");
  strcpy(blast_dir,
         "/home/jpei/promals/src_structure/structure_db/promals/blastresult");
  // 1. get the subaligns with one sequence: oneseqaln
  oneseqaln.push_back(NULL);
  for (i = 1; i <= nhits; i++) {
    // the name of the fasta file containing the single sequence of the domain
    strcpy(domain_sequence_name,
           "/home/jpei/promals/src_structure/structure_db/promals/blastresult");
    sprintf(domain_sequence_name, "%s/%s.fa", domain_sequence_name, id[i]);
    // cout << "domain sequence name: " << domain_sequence_name << endl;
    ifstream fp(domain_sequence_name, ios::in);
    fp >> domain_sequence;
    fp >> domain_sequence;
    // cout << "domain_sequence: " << domain_sequence << endl;
    // cout << "id: " << id[i] << endl;
    fp.close();
    oneseqaln.push_back(oneSeq2subalign(domain_sequence, id[i]));
  }

  // 2. get the secondary structure profile for oneseqaln
  for (i = 1; i <= nhits; i++) {
    // oneseqaln[i]->repres_name = new char [strlen(oneseqaln[i]->aname[0])+2];
    // strcpy(oneseqaln[i]->repres_name, oneseqaln[i]->aname[0]);
    // cout << oneseqaln[i]->repres_name << endl;
    oneseqaln[i]->select_representative();
    oneseqaln[i]->get_ss_prof1(blast_dir, id[i], runpsipred_command);
  }

  // 3. get the auxilary alignments: prof, read from the blast result directory
  prof.push_back(NULL);
  for (i = 1; i <= nhits; i++) {
    sprintf(blast_aln_file, "%s/%s.aln", blast_dir, id[i]);
    cout << blast_aln_file << endl;
    tmpaln =
        read_blastpgp_alignment(blast_aln_file, id[i], oneseqaln[i]->aseq[0]);
    // cout << "here" << endl;
    prof.push_back(tmpaln);
  }

  // 4. transfer the secondary structure profile from oneseqaln to prof
  for (i = 1; i <= nhits; i++) {
    prof[i]->ss = oneseqaln[i]->ss;
  }
  // cout << "here" << endl;

  // 5. get the numerical profile
  for (i = 1; i <= nhits; i++) {
    prof[i]->prof(1);
    // cout << "here" << endl;
    taln = prof[i];
    taln->prof_freq = dmatrix(taln->prof_len, 20);
    // cout << "here" << endl;
    prof[i]->pseudoCounts(taln->prof_effn, taln->prof_nef, taln->prof_len,
                          taln->prof_freq, params.aa_pair1, params.aa_bg1,
                          taln->ss->sstype);
    // cout << "here" << endl;
    prof[i]->log_pseudoCounts();
    // cout << "here" << endl;
  }
  // cout << "here" << endl;
}

void seq_str_aln::get_prof_update() {
  int i, j;
  char domain_sequence[5000];
  char domain_sequence_name[200];
  char blast_dir[300];
  char blast_aln_file[300];
  subalign *tmpaln, *taln;
  hmm_psipred_parameters params(psipred_env_number);
  params.read_parameters(psipred_parameter_file, psipred_env_number, 1);
  char subdir[10];
  subdir[3] = '\0';

  char blastresults[300], blastresults_[300];
  sprintf(blastresults, "%s/db/structure_db", program_dir);
  sprintf(blastresults_, "%s/db/structure_db/", program_dir);

  // the directory containing all the SCOP domain sequences in fasta format
  // each file contains a single domain sequence
  // strcpy(domain_sequence_name,
  // "/home/jpei/promals/src_structure/structure_db/promals_new/blastresults");
  // strcpy(blast_dir,
  // "/home/jpei/promals/src_structure/structure_db/promals_new/blastresults/");
  strcpy(domain_sequence_name, blastresults);
  strcpy(blast_dir, blastresults_);
  // 1. get the subaligns with one sequence: oneseqaln
  oneseqaln.push_back(NULL);
  for (i = 1; i <= nhits; i++) {
    strncpy(subdir, id[i], 3);
    subdir[3] = '\0';
    // the name of the fasta file containing the single sequence of the domain
    // strcpy(domain_sequence_name,
    // "/home/jpei/promals/src_structure/structure_db/promals_new/blastresults");
    strcpy(domain_sequence_name, blastresults);
    sprintf(domain_sequence_name, "%s/%s/%s.fa", domain_sequence_name, subdir,
            id[i]);
    // cout << "domain sequence name: " << domain_sequence_name << endl;
    ifstream fp(domain_sequence_name, ios::in);
    fp >> domain_sequence;
    fp >> domain_sequence;
    fp >> domain_sequence;  // add one line, since the fasta record defline now
                            // contains two fields
    // cout << "domain_sequence: " << domain_sequence << endl;
    // cout << "id: " << id[i] << endl;
    fp.close();
    oneseqaln.push_back(oneSeq2subalign(domain_sequence, id[i]));
  }

  // 2. get the secondary structure profile for oneseqaln
  for (i = 1; i <= nhits; i++) {
    // oneseqaln[i]->repres_name = new char [strlen(oneseqaln[i]->aname[0])+2];
    // strcpy(oneseqaln[i]->repres_name, oneseqaln[i]->aname[0]);
    // cout << oneseqaln[i]->repres_name << endl;
    oneseqaln[i]->select_representative();
    strncpy(subdir, id[i], 3);
    subdir[3] = '\0';
    // cout << "subdir: " << subdir << endl;
    // strcpy(blast_dir,
    // "/home/jpei/promals/src_structure/structure_db/promals_new/blastresults/");
    strcpy(blast_dir, blastresults_);
    strcat(blast_dir, subdir);
    // cout << "blast_dir: " << blast_dir << endl;
    oneseqaln[i]->get_ss_prof1(blast_dir, id[i], runpsipred_command);
  }

  char tmpid[5000];
  int *is_scopid = ivector(nhits);
  // 3.1. check the identity of the pdb, see if it is a scop id or pdbid
  for (i = 1; i <= nhits; i++) {
    // check the identity of the pdb, see if it is a scop id or pdbid
    strncpy(subdir, id[i], 3);
    subdir[3] = '\0';
    // sprintf(blast_aln_file, "%s/%s/%s.fa",
    // "/home/jpei/promals/src_structure/structure_db/promals_new/blastresults",
    // subdir, id[i]);
    sprintf(blast_aln_file, "%s/%s/%s.fa", blastresults, subdir, id[i]);
    ifstream tmpfa_fp(blast_aln_file, ios::in);
    tmpfa_fp >> tmpid;
    // cout << tmpid << endl;
    tmpfa_fp >> tmpid;
    // cout << tmpid << endl;
    tmpfa_fp.close();
    if (tmpid[0] == 'd')
      is_scopid[i] = 1;
    else
      is_scopid[i] = 0;
    // do not consider scop domains as a whole sequence, also consider only
    // aligned portion
    is_scopid[i] = 0;
    // cout << "nhit: " << i << " is_scopid: " << is_scopid[i] << endl;
  }

  // popen("ls", "r");

  // 3. get the auxilary alignments: prof, read from the blast result directory
  prof.push_back(NULL);
  for (i = 1; i <= nhits; i++) {
    strncpy(subdir, id[i], 3);
    subdir[3] = '\0';
    // sprintf(blast_aln_file, "%s/%s/%s.aln",
    // "/home/jpei/promals/src_structure/structure_db/promals_new/blastresults",
    // subdir, id[i]);
    sprintf(blast_aln_file, "%s/%s/%s.aln", blastresults, subdir, id[i]);
    // cout << "blastalnfile: " << blast_aln_file << endl;
    tmpaln =
        read_blastpgp_alignment(blast_aln_file, id[i], oneseqaln[i]->aseq[0]);
    if (!is_scopid[i]) {
      int *mark = ivector(tmpaln->nal);
      int *mark1 = ivector(tmpaln->alilen);
      for (j = 1; j <= tmpaln->nal; j++) mark[j] = 1;
      for (j = 1; j <= tmpaln->alilen; j++) {
        if (j < str_start[i])
          mark1[j] = 0;
        else if (j > str_end[i])
          mark1[j] = 0;
        else
          mark1[j] = 1;
      }
      taln = tmpaln->sub2align(mark, mark1);
      delete tmpaln;
      delete[] mark;
      delete[] mark1;
      tmpaln = taln;
      // tmpaln->printali(80);
    }
    // tmpaln->printali(80);
    prof.push_back(tmpaln);
  }

  // 4. transfer the secondary structure profile from oneseqaln to prof
  for (i = 1; i <= nhits; i++) {
    // prof[i]->ss = oneseqaln[i]->ss;
    if (is_scopid[i]) {
      prof[i]->ss = oneseqaln[i]->ss;
      // prof[i]->ss->print_result();
    } else {
      // oneseqaln[i]->ss->print_result();
      // cout << "Debug here" << endl;
      prof[i]->ss = oneseqaln[i]->ss->sub_ss_prof(str_start[i], str_end[i]);
      // cout << "Debug here" << endl;
      delete oneseqaln[i]->ss;
      // cout << "Debug here" << endl;
      // prof[i]->ss->print_result();
    }
  }
  // cout << "here" << endl;

  // 5. get the numerical profile
  for (i = 1; i <= nhits; i++) {
    prof[i]->prof(1);
    // cout << "here" << endl;
    taln = prof[i];
    taln->prof_freq = dmatrix(taln->prof_len, 20);
    // cout << "here" << endl;
    prof[i]->pseudoCounts(taln->prof_effn, taln->prof_nef, taln->prof_len,
                          taln->prof_freq, params.aa_pair1, params.aa_bg1,
                          taln->ss->sstype);
    // cout << "here" << endl;
    // prof[i]->log_pseudoCounts();
    // cout << "here" << endl;
  }
  // cout << "here" << endl;
}

void seq_str_aln::print_result() {
  int i, j;

  cout << "Sequence structure aln " << endl;
  for (i = 1; i < id.size(); i++) {
    cout << "id: " << id[i] << endl;
    cout << "start: " << start[i] << " end: " << end[i] << endl;
    cout << "str_start: " << str_start[i] << " str_end: " << str_end[i] << endl;
    cout << "slen: " << slen[i] << endl;
    cout << "aln: ";
    for (j = 1; j <= slen[i]; j++) {
      cout << j << ":" << aln[i][j] << ";  ";
    }
    cout << endl;
    prof[i]->printali(80);
  }
  cout << endl;
  cout << endl;
}
