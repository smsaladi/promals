#include "subalign.h"
#include "amino.h"
#include "header_cpp.h"
#include "regularizer.h"
#include "util.h"
//#include "hmm_psipred.h"

static int mydebug = 11;
static int debug = 11;
static void run_psipred_check_lowercase_letters(char *repres_name,
                                                char *runpsipred1_command);

// Initilize static data member at file scope
int subalign::subaligncount = 0;

// Default constructor
subalign::subalign() {
  alignment = NULL;
  aseq = NULL;
  aname = NULL;
  astart = NULL;
  nal = 0;
  alilen = 0;

  gap_threshold = 0.5;
  gapt_threshold = 1.0;
  done_profile = 0;
  done_prof = 0;
  done_prof_freq = 0;
  beta = 10;
  prof_alphabet1 = 0;
  ss = 0;
  score_bg_aa = 0;
  score_bg_ss = 0;
  bg_type_here = -1;
  done_score_bg = 0;   // sum-of-pair measure
  done_score_bg2 = 0;  // mutual emission measure

  apos_filtr = 0;
  n_effAa = 0;
  pseudoCnt = 0;
  sum_eff_let = 0;
  maskgapRegion = 0;
  gap_content = 0;

  hwt_all = 0;
  distMat = 0;

  prof_pos = 0;
  prof_effn = 0;
  prof_freq = 0;
  prof_sum_eff = 0;
  prof_hwt_all = 0;
  prof_gap_content = 0;

  ss = 0;
  repres_name = 0;
  prof_alphabet1 = 0;
  score_bg_aa = 0;
  score_bg_ss = 0;

  subaligncount++;
}

void subalign::set_gapt(double gapt) { gapt_threshold = gapt; }

// Copy constructor
subalign::subalign(const subalign &init) {
  int i, j;

  nal = init.nal;
  alilen = init.alilen;
  mnamelen = init.mnamelen;
  beta = init.beta;

  // cout << "nal: " << nal << " alilen: " << alilen << " mnamelen " << mnamelen
  // << endl;

  alignment = imatrix(nal, alilen);
  aname = cmatrix(nal, mnamelen + 1);
  aseq = cmatrix(nal, alilen + 1);

  for (i = 1; i <= nal; i++) {
    for (j = 1; j <= alilen; j++) {
      // cout << i << endl;
      alignment[i][j] = init.alignment[i][j];
      aseq[i - 1][j - 1] = init.aseq[i - 1][j - 1];
      // cout << alignment[i][j];
    }
    aseq[i - 1][alilen] = '\0';
    // cout << endl;
  }
  for (i = 0; i < nal; i++) {
    strcpy(aname[i], init.aname[i]);
  }

  gap_threshold = 0.5;
  gapt_threshold = 1.0;
  done_profile = 0;
  done_prof = 0;
  done_prof_freq = 0;
  beta = 10;
  prof_alphabet1 = 0;
  ss = 0;
  score_bg_aa = 0;
  score_bg_ss = 0;
  bg_type_here = -1;
  done_score_bg = 0;
  done_score_bg2 = 0;

  apos_filtr = 0;
  n_effAa = 0;
  pseudoCnt = 0;
  sum_eff_let = 0;
  maskgapRegion = 0;
  gap_content = 0;

  hwt_all = 0;
  distMat = 0;

  prof_pos = 0;
  prof_effn = 0;
  prof_freq = 0;
  prof_sum_eff = 0;
  prof_hwt_all = 0;
  prof_gap_content = 0;

  ss = 0;
  repres_name = 0;
  prof_alphabet1 = 0;
  score_bg_aa = 0;
  score_bg_ss = 0;

  subaligncount++;
  // cout << "alignment count: " << subaligncount << endl;
}

// Constructor reading an input alignment file in fasta format
subalign::subalign(char *filename, const char *format_name, int max_name_len) {
  int i, j, k;
  int maxstr = 100001;
  char line[100001];
  char *s, *ss;

  if (strncmp(format_name, "fasta", 5) != 0) {
    cout << "currently support fasta format" << endl;
    exit(0);
  }

  ifstream fp(filename, ios::in);
  int count = 0;
  vector<string> names;
  vector<string> seqs;
  string tmpseq = "";
  names.push_back(string(""));
  while (fp.good()) {
    fp.getline(line, maxstr);
    for (s = line; isspace(*s); s++)
      ;
    if (strlen(s) == 0) continue;  // empty line
    if (*s == '>') {
      count++;
      for (ss = s + 1; isspace(*ss); ss++)
        ;
      names.push_back(string(ss).substr(0, max_name_len));
      seqs.push_back(tmpseq);
      tmpseq = "";
      continue;
    } else {
      for (ss = s; isspace(*ss); ss++)
        ;
      tmpseq += ss;
    }
  }
  nal = count;
  seqs.push_back(tmpseq);
  int maxlen = 0;
  for (i = 1; i < int(seqs.size()); i++) {
    // cout << seqs[i] << endl;
    if (int(seqs[i].size()) > maxlen) maxlen = seqs[i].size();
  }
  alilen = maxlen;
  mnamelen = 0;
  for (i = 1; i < int(names.size()); i++) {
    // cout << names[i] << endl;
    if (int(names[i].size()) > mnamelen) mnamelen = names[i].size();
  }
  // cout << "mnamelen: " << mnamelen << endl;
  // fill the ends with gap characters
  for (i = 1; i < int(seqs.size()); i++) {
    for (j = 1; j <= alilen - seqs[i].size(); j++) {
      seqs[i] += "-";
    }
  }
  aname = cmatrix(nal, mnamelen + 1);
  aseq = cmatrix(nal, alilen + 1);
  for (i = 0; i < nal; i++) {
    strcpy(aname[i], names[i + 1].c_str());
    strcpy(aseq[i], seqs[i + 1].c_str());
  }
  // new
  // reassign();

  alignment = imatrix(nal, alilen);

  for (i = 1; i <= nal; i++) {
    for (j = 1; j <= alilen; j++) {
      alignment[i][j] = am2num(aseq[i - 1][j - 1]);
    }
  }

  gap_threshold = 0.5;
  gapt_threshold = 1.0;
  done_profile = 0;
  done_prof = 0;
  done_prof_freq = 0;
  beta = 10;
  prof_alphabet1 = 0;

  ss = 0;
  score_bg_aa = 0;
  score_bg_ss = 0;
  bg_type_here = -1;
  done_score_bg = 0;
  done_score_bg2 = 0;

  apos_filtr = 0;
  n_effAa = 0;
  pseudoCnt = 0;
  sum_eff_let = 0;
  maskgapRegion = 0;
  gap_content = 0;

  hwt_all = 0;
  distMat = 0;

  prof_pos = 0;
  prof_effn = 0;
  prof_freq = 0;
  prof_sum_eff = 0;
  prof_hwt_all = 0;
  prof_gap_content = 0;

  ss = 0;
  repres_name = 0;
  prof_alphabet1 = 0;
  score_bg_aa = 0;
  score_bg_ss = 0;
  subaligncount++;
}

// Constructor reading an input alignment file
subalign::subalign(char *filename) {
  int i, j, k;

  readali(filename);

  // new
  reassign();

  alignment = imatrix(nal, alilen);

  for (i = 1; i <= nal; i++) {
    for (j = 1; j <= alilen; j++) {
      alignment[i][j] = am2num(aseq[i - 1][j - 1]);
    }
  }

  gap_threshold = 0.5;
  gapt_threshold = 1.0;
  done_profile = 0;
  done_prof = 0;
  done_prof_freq = 0;
  beta = 10;
  prof_alphabet1 = 0;

  ss = 0;
  score_bg_aa = 0;
  score_bg_ss = 0;
  bg_type_here = -1;
  done_score_bg = 0;
  done_score_bg2 = 0;

  apos_filtr = 0;
  n_effAa = 0;
  pseudoCnt = 0;
  sum_eff_let = 0;
  maskgapRegion = 0;
  gap_content = 0;

  hwt_all = 0;
  distMat = 0;

  prof_pos = 0;
  prof_effn = 0;
  prof_freq = 0;
  prof_sum_eff = 0;
  prof_hwt_all = 0;
  prof_gap_content = 0;

  ss = 0;
  repres_name = 0;
  prof_alphabet1 = 0;
  score_bg_aa = 0;
  score_bg_ss = 0;
  subaligncount++;
}

// Constructor reading an input file name
subalign::subalign(string filename) {
  int i, j, k;

  char filename1[100];
  for (i = 0; i < filename.length(); i++) {
    filename1[i] = filename[i];
  }
  filename1[i] = '\0';
  readali(filename1);

  // new
  reassign();

  alignment = imatrix(nal, alilen);

  for (i = 1; i <= nal; i++) {
    for (j = 1; j <= alilen; j++) {
      alignment[i][j] = am2num(aseq[i - 1][j - 1]);
    }
  }

  gap_threshold = 0.5;
  gapt_threshold = 1.0;
  done_profile = 0;
  done_prof = 0;
  done_prof_freq = 0;
  beta = 10;

  prof_alphabet1 = 0;
  ss = 0;
  score_bg_aa = 0;
  score_bg_ss = 0;
  bg_type_here = -1;
  done_score_bg = 0;
  done_score_bg2 = 0;

  apos_filtr = 0;
  n_effAa = 0;
  pseudoCnt = 0;
  sum_eff_let = 0;
  maskgapRegion = 0;
  gap_content = 0;

  hwt_all = 0;
  distMat = 0;

  prof_pos = 0;
  prof_effn = 0;
  prof_freq = 0;
  prof_sum_eff = 0;
  prof_hwt_all = 0;
  prof_gap_content = 0;

  ss = 0;
  repres_name = 0;
  prof_alphabet1 = 0;
  score_bg_aa = 0;
  score_bg_ss = 0;
  subaligncount++;
}

// Destructor
subalign::~subalign() {
  int i;

  if (alignment) {
    for (i = 0; i <= nal; i++) delete[] alignment[i];
    delete[] alignment;
  }
  if (aseq) {
    for (i = 0; i <= nal; i++) {
      delete[] aseq[i];
    }
    delete[] aseq;
  }
  if (aname) {
    for (i = 0; i <= nal; i++) delete[] aname[i];
    delete[] aname;
  }
  subaligncount--;
  // cout << "alignment count decrease: " << subaligncount << endl;

  if (done_profile) {
    delete[] apos_filtr;
    for (i = 0; i <= alilen; i++) delete[] n_effAa[i];
    delete[] n_effAa;
    for (i = 0; i <= alilen; i++) delete[] pseudoCnt[i];
    delete[] pseudoCnt;
    delete[] sum_eff_let;
    delete[] maskgapRegion;
  }
  if (hwt_all) delete[] hwt_all;
  if (done_prof) {
    delete[] prof_pos;
    delete[] prof_sum_eff;
    delete[] prof_hwt_all;
    delete[] prof_gap_content;
    free_dmatrix(prof_effn, prof_len, 20);
  }
  if (done_prof_freq) {
    free_dmatrix(prof_freq, prof_len, 20);
  }
  if (score_bg_aa) delete[] score_bg_aa;
  if (score_bg_ss) delete[] score_bg_ss;
  if (repres_name) delete[] repres_name;
  if (prof_alphabet1) delete[] prof_alphabet1;
}

// Overload assignment operator
const subalign &subalign::operator=(const subalign &right) {
  int i, j, k;

  if (&right != this) {
    // delete the original arrays
    if (alignment) {
      for (i = 0; i <= nal; i++) {
        delete[] alignment[i];
      }
      delete[] alignment;
    }
    if (aseq) {
      for (i = 0; i <= nal; i++) {
        delete[] aseq[i];
      }
      delete[] aseq;
    }
    if (aname) {
      for (i = 0; i <= nal; i++) {
        delete[] aname[i];
      }
      delete[] aname;
    }

    // assignment of new values
    nal = right.nal;
    alilen = right.alilen;
    mnamelen = right.mnamelen;
    alignment = imatrix(nal, alilen);
    aseq = cmatrix(nal, alilen);
    aname = cmatrix(nal, mnamelen);
    for (i = 1; i <= nal; i++) {
      for (j = 1; j <= alilen; j++) {
        alignment[i][j] = right.alignment[i][j];
        aseq[i - 1][j - 1] = right.aseq[i - 1][j - 1];
      }
      aseq[i - 1][alilen] = '\0';
      strcpy(aname[i - 1], right.aname[i - 1]);
    }
    beta = right.beta;
  }
  return *this;
}

// Return the number of subalign objects intantiated
int subalign::getsubaligncount() { return subaligncount; }

/******************************************************************/
/*adapted from align.c, auxilary routines used for readali...*/
/*memory allocations follow the C language constoms */

void *subalign::mymalloc(int size) {
  void *buf;

  if ((buf = malloc(size)) == NULL) {
    fprintf(stderr, "Not enough memory: %d\n", size);
    exit(1);
  }
  return buf;
}

char *subalign::strsave(char *str) {
  char *buf;
  int l;

  l = strlen(str);
  buf = (char *)mymalloc((l + 1) * sizeof(str[0]));
  strcpy(buf, str);
  return buf;
}

char *subalign::strnsave(char *str, int l) {
  char *buf;

  buf = (char *)mymalloc(l + 1);
  memcpy(buf, str, l);
  buf[l] = '\0';
  return buf;
}

char **subalign::incbuf(int n, char **was) {
  char **buf;

  buf = (char **)mymalloc((n + 1) * sizeof(buf[0]));
  if (n > 0) {
    memcpy(buf, was, n * sizeof(was[0]));
    free(was);
  }
  buf[n] = NULL;
  return buf;
}

int *subalign::incibuf(int n, int *was) {
  int *ibuf;

  ibuf = (int *)mymalloc((n + 1) * sizeof(ibuf[0]));
  if (n > 0) {
    memcpy(ibuf, was, n * sizeof(was[0]));
    free(was);
  }
  ibuf[n] = 0;
  return ibuf;
}

/* end of auxilary routines from align.c ****************************/
/********************************************************************/

int subalign::getNal() { return (nal); }

int subalign::getAlilen() { return (alilen); }

int **subalign::getAlignment() {
  /* int **ali;
  int i, j;

  ali = imatrix(0, nal+1, 0, alilen+1);
  for(j=1;j<=nal;j++)
  for(i=1;i<=alilen;i++) {
          ali[j][i] = alignment[j][i];
  } */

  return alignment;
}

char **subalign::getAseq() {
  /* char **seq;
  int i, j;

  seq = cmatrix(0, nal+1, 0, alilen+1);
  for(j=0;j<nal;j++)
  for(i=0;i<alilen;i++) {
          seq[j][i] = aseq[j][i];
  } */

  return aseq;
}

char **subalign::getAname() {
  /* char **name;
  int i, j;
  int mlen;

  for (i=1, mlen=strlen(aname[0]); i < nal; i++) {
          if (mlen < strlen(aname[i])) {
                  mlen = strlen(aname[i]);
          }
  }

  name = cmatrix(0,nal,0,mlen);
  for(i=0;i<nal;i++) {
          strcpy(name[i], aname[i]);
  } */

  return (aname);
}

void subalign::setNal(int n) { nal = n; }

void subalign::setAlilen(int len) { alilen = len; }

void subalign::setAlignment(int **ali) {
  int i, j;

  if (ali == NULL) {
    fprintf(stderr, "WARNING:: alignment is null\n");
    return;
  }

  alignment = new int *[nal + 1];
  for (i = 1; i <= nal; i++) alignment[i] = new int[alilen + 1];
  // fprintf(stderr, "nal: %d\n", nal);
  for (i = 1; i <= nal; i++) {
    for (j = 1; j <= alilen; j++) {
      alignment[i][j] = ali[i][j];
      //              fprintf(stderr,"%d",ali[i][j]);
    }
  }

  /*alignment = ali;*/
}

void subalign::setAseq(char **seq) {
  int i, j;

  if (seq == NULL) {
    fprintf(stderr, "WARNING:: alignment sequences are null\n");
    return;
  }

  aseq = new char *[nal];
  for (i = 0; i < nal; i++) aseq[i] = new char[alilen + 1];
  for (i = 0; i < nal; i++) {
    for (j = 0; j < alilen; j++) aseq[i][j] = seq[i][j];
    aseq[i][alilen] = '\0';
  }
}

void subalign::setAname(char **name) {
  int i, j;
  int mlen;

  if (name == NULL) {
    fprintf(stderr, "ALign name empty\n");
    return;
  }

  // fprintf(stderr, "nal: %d\n", nal);
  for (i = 0, mlen = strlen(name[0]); i < nal; i++) {
    if (mlen < strlen(name[i])) {
      mlen = strlen(name[i]);
    }
  }

  aname = new char *[nal];
  for (i = 0; i < nal; i++) aname[i] = new char[mlen + 1];
  for (i = 0; i < nal; i++) strcpy(aname[i], name[i]);
  mnamelen = mlen;
}

// Read an input alignment in ClustalW (with or without header) format
void subalign::readali(char *filename) {
  FILE *fp;
  char *s, *ss, *seqbuf;
  int n, l, len, len0;
  char str[100001];
  int MAXSTR = 100001;

  if ((fp = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "No such file: \"%s\"\n", filename);
    exit(1);
  }
  alilen = 0;
  nal = 0;
  n = 0;
  while (fgets(str, MAXSTR, fp) != NULL) {
    if (strncmp(str, "CLUSTAL ", 8) == 0) continue;
    for (ss = str; isspace(*ss); ss++)
      ;
    if ((*ss == '.') || (*ss == '*') || (*ss == ':')) continue;
    if (*ss == '\0') {
      if (n == 0) {
        continue;
      }
      if (nal == 0) {
        if (n == 0) {
          fprintf(stderr, "No alignments read: %s\n", filename);
          exit(1);
        }
        nal = n;
      } else if (n != nal) {
        fprintf(stderr, "Wrong nal, was: %d, now: %d\n", nal, n);
        exit(1);
      }
      n = 0;
      continue;
    }
    for (s = ss; *s != '\0' && !isspace(*s); s++)
      ;
    *s++ = '\0';
    if (nal == 0) {
      astart = incibuf(n, astart);
      alen = incibuf(n, alen);
      aseq = incbuf(n, aseq);
      aname = incbuf(n, aname);
      aname[n] = strsave(ss);
    } else {
      if (n < 0 || n >= nal) {
        fprintf(stderr, "Bad sequence number: %d of %d\n", n, nal);
        exit(1);
      }
      if (strcmp(ss, aname[n]) != 0) {
        fprintf(stderr, "Names do not match");
        fprintf(stderr, ", was: %s, now: %s\n", aname[n], ss);
        exit(1);
      }
    }
    for (ss = s; isspace(*ss); ss++)
      ;
    for (s = ss; isdigit(*s); s++)
      ;
    if (isspace(*s)) {
      if (nal == 0) {
        astart[n] = atoi(ss);
        /*fprintf(stderr, "n: %d  astart: %d\n", n, astart[n]);*/
      }
      for (ss = s; isspace(*ss); ss++)
        ;
    }
    for (s = ss, l = 0; *s != '\0' && !isspace(*s); s++) {
      if (isalpha(*s)) {
        l++;
      }
    }
    len = s - ss;
    if (n == 0) {
      len0 = len;
      alilen += len;
    } else if (len != len0) {
      fprintf(stderr, "wrong len for %s", aname[n]);
      fprintf(stderr, ", was: %d, now: %d\n", len0, len);
      exit(1);
    }
    alen[n] += l;
    if (aseq[n] == NULL) {
      aseq[n] = strnsave(ss, len);
    } else {
      seqbuf = (char *)mymalloc(alilen + 1);
      memcpy(seqbuf, aseq[n], alilen - len);
      free(aseq[n]);
      aseq[n] = seqbuf;
      memcpy(seqbuf + alilen - len, ss, len);
      seqbuf[alilen] = '\0';
    }
    n++;
    // cout << n << endl;
  }
  if (nal == 0) {
    if (n == 0) {
      fprintf(stderr, "No alignments read: %s\n", filename);
      exit(1);
    }
    nal = n;
  } else if (n != 0 && n != nal) {
    fprintf(stderr, "Wrong nal, was: %d, now: %d\n", nal, n);
    exit(1);
  }
  // free(astart);
  free(alen);
  fclose(fp);
}

// convert the C-style malloc allocations to C++ style new allocations
void subalign::reassign() {
  int i, j, k;
  int mlen;
  char **tmpaseq;
  char **tmpaname;

  tmpaseq = cmatrix(nal, alilen);

  for (i = 1, mlen = strlen(aname[0]); i < nal; i++) {
    if (mlen < strlen(aname[i])) {
      mlen = strlen(aname[i]);
    }
  }
  mnamelen = mlen;
  tmpaname = cmatrix(nal, mlen + 1);

  for (i = 0; i < nal; i++) {
    strcpy(tmpaseq[i], aseq[i]);
    strcpy(tmpaname[i], aname[i]);
    free(aseq[i]);
    free(aname[i]);
  }
  free(aseq);
  free(aname);
  aseq = cmatrix(nal, alilen);
  aname = cmatrix(nal, mlen + 1);
  for (i = 0; i < nal; i++) {
    strcpy(aseq[i], tmpaseq[i]);
    strcpy(aname[i], tmpaname[i]);
    aseq[i][alilen] = '\0';
    delete[] tmpaseq[i];
    delete[] tmpaname[i];
  }
  delete[] tmpaseq;
  delete[] tmpaname;
}

// Print alignment with a given block size
void subalign::printali(int blocksize, int header) {
  int i, j, k;

  if (header) {
    cout << endl << endl;
    cout << "blocksize: " << blocksize << endl;
    cout << "alilen: " << alilen << endl;
    cout << "nal: " << nal << endl;
  }
  int nblocks = (alilen - 1) / blocksize + 1;
  if (header) {
    cout << "nblocks: " << nblocks << endl;
  }
  cout << endl;

  for (i = 1; i <= nblocks; i++) {
    for (j = 0; j < nal; j++) {
      cout << setiosflags(ios::left) << setw(mnamelen + 3) << aname[j];
      for (k = 1; k <= blocksize; k++) {
        if ((i - 1) * blocksize + k - 1 >= alilen) break;
        cout << aseq[j][(i - 1) * blocksize + k - 1];
      }
      cout << endl;
    }
    cout << endl << endl;
  }
}

// Print alignment with a given block size
void subalign::printali(int blocksize) {
  int i, j, k;

  cout << endl << endl;
  cout << "blocksize: " << blocksize << endl;
  cout << "alilen: " << alilen << endl;
  cout << "nal: " << nal << endl;
  int nblocks = (alilen - 1) / blocksize + 1;
  cout << "nblocks: " << nblocks << endl;
  cout << endl;

  for (i = 1; i <= nblocks; i++) {
    for (j = 0; j < nal; j++) {
      cout << setiosflags(ios::left) << setw(mnamelen + 3) << aname[j];
      for (k = 1; k <= blocksize; k++) {
        if ((i - 1) * blocksize + k - 1 >= alilen) break;
        cout << aseq[j][(i - 1) * blocksize + k - 1];
      }
      cout << endl;
    }
    cout << endl << endl;
  }
}

// Print alignment to a file
void subalign::printali(const char *filename, int blocksize) {
  int i, j, k;
  ofstream outfile(filename, ios::out);
  if (!outfile) {
    cout << "cannot write the alignment to " << filename << endl;
    exit(1);
  }
  outfile << "CLUSTAL format multiple sequence alignment by PROMALS" << endl;
  outfile << endl << endl;
  int nblocks = (alilen - 1) / blocksize + 1;

  for (i = 1; i <= nblocks; i++) {
    for (j = 0; j < nal; j++) {
      outfile << setiosflags(ios::left) << setw(mnamelen + 3) << aname[j];
      for (k = 1; k <= blocksize; k++) {
        if ((i - 1) * blocksize + k - 1 >= alilen) break;
        outfile << aseq[j][(i - 1) * blocksize + k - 1];
      }
      outfile << endl;
    }
    outfile << endl << endl;
  }
  outfile.close();
}

/* from a given alignment with aa as numbers, computes effective aa counts
(PSIC->our formula) and marks the columns with EFFECTIVE content of gaps >
threshold (effgapmax) */
void subalign::neffsForEachCol_maskGapReg(int **ali, int n, int len,
                                          double effgapmax,
                                          double effgapRegionMin,
                                          double **n_effAa, double *sum_eff_let,
                                          int *maskgapRegion, int *apos_filtr,
                                          int *len_lowgaps, double *nef) {
  int i, j, k, l;
  int alilen_mat, nsymbols_col, nsymbols;
  double nef_loc;
  int ele;
  double effnu[21];
  double sum_let;
  int *mark;
  int flagmark;

  gap_content = new double[alilen + 1];

  mark = new int[n + 11];
  alilen_mat = 0;
  nsymbols = 0;
  for (j = 1; j <= len; j++) {
    nsymbols_col = 0;
    sum_let = 0;

    for (k = 0; k <= 20; ++k) {
      /* Mark sequences that have amino acid  k (or gap, k=0) in this jth
       * position */
      flagmark = 0;
      for (i = 1; i <= n; ++i) {
        mark[i] = 0;

        ele = ali[i][j];
        if (ele == k) {
          mark[i] = 1;
          flagmark = 1;
        }
        ele = ali[i][j] - 25;
        if (ele == k) {
          mark[i] = 1;
          flagmark = 1;
        }
      }

      /* If aa k (or gap) is present in this position call compute k-th
       * effective count */
      if (flagmark == 1) {
        effnu[k] = effective_number_nogaps(ali, mark, n, 1, len);
        nsymbols_col++;

      } else {
        effnu[k] = 0.0;
      }

      if (k > 0) sum_let += effnu[k];
    }

    // calculate raw gap fraction
    if (mydebug > 100) {
      int gap_count = 0;
      for (i = 1; i <= n; i++) {
        if (ali[i][j] == 0) gap_count++;
      }
      double gapfract = 1.0 * gap_count / n;
      // cout << j << "  " << 1.0*effnu[0]/(sum_let + effnu[0]) << "  " <<
      // gapfract << endl;
    }

    if (mydebug > 100)
      cout << j << "  " << aseq[0][j - 1] << " " << sum_let << "  " << effnu[0]
           << endl;

    if (sum_let > 0 && 1.0 * effnu[0] / (sum_let + effnu[0]) < effgapmax) {
      alilen_mat++;
      gap_content[alilen_mat] = 1.0 * effnu[0] / (sum_let + effnu[0]);
      for (k = 0; k <= 20; k++) {
        n_effAa[alilen_mat][k] = effnu[k];
      }
      sum_eff_let[alilen_mat] = sum_let;
      apos_filtr[alilen_mat] = j;
      nsymbols += nsymbols_col;

      if (1.0 * effnu[0] / (sum_let + effnu[0]) < effgapRegionMin) {
        maskgapRegion[alilen_mat] = 0;
      } else {
        maskgapRegion[alilen_mat] = 1;
      }
    }
  }

  // cout << nsymbols << "\t" << alilen_mat << endl;

  nef_loc = 1.0 * nsymbols / alilen_mat;
  *nef = nef_loc;
  *len_lowgaps = alilen_mat;

  average_sum_eff_let = 0;
  for (i = 1; i <= alilen_mat; i++) {
    average_sum_eff_let += sum_eff_let[i];
  }
  average_sum_eff_let /= alilen_mat;

  maskgapRegion[0] = maskgapRegion[alilen_mat + 1] = 0;

  delete[] mark;

  /* for(i=1;i<=alilen_mat;i++) {
      for(j=1;j<=20;j++) {
          cout << n_effAa[i][j] << " ";
      }
      cout << endl;
  } */
}

double subalign::effective_number_nogaps(int **ali, int *marks, int n,
                                         int start, int end) {
  /* from the alignment of n sequences ali[1..n][1..l]
  calculates effective number of sequences that are marked by 1 in mark[1..n]
  for the segment of positions ali[][start..end]
  Neff=ln(1-0.05*N-of-different-letters-per-site)/ln(0.95)
  */

  int i, k, a, flag;
  int amco[21], lettercount = 0, sitecount = 0;
  double letpersite = 0, neff;

  for (k = start; k <= end; ++k) {
    /******************DUMP the condition "consider only positions without gaps
    *in the marked seqs" ***
    ********/
    /*****  flag=0;for(i=1;i<=n;++i)if(marks[i]==1 && ali[i][k]==0)flag=1;
            if(flag==1)continue;
    *****/
    for (a = 0; a <= 20; ++a) amco[a] = 0;
    for (i = 1; i <= n; ++i)
      if (marks[i] == 1) amco[ali[i][k]]++;
    flag = 0;
    for (a = 1; a <= 20; ++a)
      if (amco[a] > 0) {
        flag = 1;
        lettercount++;
      }
    if (flag == 1) sitecount++;
  }
  if (sitecount == 0)
    letpersite = 0;
  else
    letpersite = 1.0 * lettercount / sitecount;

  neff = -log(1.0 - 0.05 * letpersite) / 0.05129329438755;

  return neff;
}

void subalign::log_pseudoCounts() {
  int i, j;

  for (i = 1; i <= prof_len; i++) {
    for (j = 1; j <= 20; j++) {
      prof_freq[i][j] = log(prof_freq[i][j]);
    }
  }
}

void subalign::pseudoCounts(double **matrix, double n_eff, int len,
                            double **pseudoCnt0) {
  int i, j, k;
  double f[21], g[21];
  double sumN;
  double alpha;

  alpha = n_eff - 1;
  if (mydebug > 1) cout << "n_eff: " << n_eff << endl;
  if (mydebug > 1) cout << "alpha: " << alpha << endl;
  if (mydebug > 1) cout << "beta: " << beta << endl;

  if (beta == 0) alpha = 1;

  for (i = 1; i <= len; i++) {
    sumN = 0;
    for (j = 1; j <= 20; j++) sumN += matrix[i][j];
    for (j = 1; j <= 20; j++) {
      f[j] = 1.0 * matrix[i][j] / sumN;
    }
    for (j = 1; j <= 20; j++) {
      g[j] = 0;
      for (k = 1; k <= 20; k++) {
        g[j] += q_blosum62[j][k] * f[k] / robinson_freq[k];
        if (mydebug > 1)
          cout << "k: " << k << " " << q_blosum62[j][k] << " " << f[k] << endl;
      }
      pseudoCnt0[i][j] = (alpha * f[j] + beta * g[j]) / (alpha + beta);
      if (mydebug > 1)
        cout << j << " " << g[j] << " " << pseudoCnt0[i][j] << endl;
    }
  }

  // if(mydebug>-1) cout << "+++++++++++" << endl;
}

void subalign::pseudoCounts(double **matrix, double n_eff, int len,
                            double **pseudoCnt0, float **input_matrix,
                            float *input_bfreq) {
  int i, j, k;
  double f[21], g[21];
  double sumN;
  double alpha;

  alpha = n_eff - 1;
  // cout << n_eff << endl;
  // cout << alpha << endl;
  // cout << beta <<endl;

  if (beta == 0) alpha = 1;

  for (i = 1; i <= len; i++) {
    sumN = 0;
    for (j = 1; j <= 20; j++) sumN += matrix[i][j];
    for (j = 1; j <= 20; j++) {
      f[j] = 1.0 * matrix[i][j] / sumN;
    }
    for (j = 1; j <= 20; j++) {
      g[j] = 0;
      for (k = 1; k <= 20; k++) {
        g[j] += input_matrix[j][k] * f[k] / input_bfreq[k];
        // cout << i << " " << j << " " << k << " " << input_matrix[j][k] << " "
        // << input_bfreq[k] << " " << g[j] << endl;
      }
      pseudoCnt0[i][j] = (alpha * f[j] + beta * g[j]) / (alpha + beta);
    }
  }
  // if(mydebug>-1)cout << "=======++++" << endl;
}

void subalign::pseudoCounts(double **matrix, double n_eff, int len,
                            double **pseudoCnt0, float ***input_matrices,
                            float **input_bfreqs, int *mat_selections) {
  int i, j, k;
  double f[21], g[21];
  double sumN;
  double alpha;

  alpha = n_eff - 1;
  // cout << "neff: " << n_eff << endl;
  // cout << alpha << endl;
  // cout << beta <<endl;

  if (beta == 0) alpha = 1;

  for (i = 1; i <= len; i++) {
    sumN = 0;
    for (j = 1; j <= 20; j++) sumN += matrix[i][j];
    for (j = 1; j <= 20; j++) {
      f[j] = 1.0 * matrix[i][j] / sumN;
      // cout << "j " << j << " " << f[j] << endl;
    }
    for (j = 1; j <= 20; j++) {
      g[j] = 0;
      for (k = 1; k <= 20; k++) {
        // cout << i << " " << j << " " << k << " " << mat_selections[i] <<
        // endl; //" " << input_matrices[mat_selections[i]][j][k] << endl; //" "
        // << input_bfreqs[mat_selections[i]][k] << " " << g[j] << endl;
        g[j] += input_matrices[mat_selections[i]][j][k] * f[k] /
                input_bfreqs[mat_selections[i]][k];
      }
      pseudoCnt0[i][j] = (alpha * f[j] + beta * g[j]) / (alpha + beta);
      // cout << pseudoCnt0[i][j] << endl;
    }
  }
}

void subalign::profile() {
  int i, j, k, m, n;

  apos_filtr = new int[alilen + 1];
  n_effAa = new double *[alilen + 1];
  pseudoCnt = new double *[alilen + 1];
  double *tmp = new double[15];
  fflush(stdout);
  for (i = 0; i <= alilen; i++) {
    n_effAa[i] = new double[20 + 1];
    pseudoCnt[i] = new double[20 + 1];
  }
  sum_eff_let = new double[alilen + 1];
  maskgapRegion = new int[alilen + 2];

  // for(i=0;i<=20;i++) for(j=0;j<=20;j++) qmatrix[i][j] = q_blosum62[i][j];

  neffsForEachCol_maskGapReg(alignment, nal, alilen, gap_threshold,
                             gapt_threshold, n_effAa, sum_eff_let,
                             maskgapRegion, apos_filtr, &alilen_mat, &n_eff);

  pseudoCounts(n_effAa, n_eff, alilen_mat, pseudoCnt);

  done_profile = 1;

  /* for(i=1;i<=alilen_mat;i++) {
      for(j=1;j<=20;j++) {
          cout << n_effAa[i][j] << " ";
      }
      cout << endl;
  }

  cout << endl;
  for(i=1;i<=alilen_mat;i++) {
      for(j=1;j<=20;j++) {
          cout << pseudoCnt[i][j] << " ";
      }
      cout << endl;
  } */
}

void subalign::h_weight_all() {
  int i, j, k;
  int position_num_letters;
  int tmp_mark[21];
  hwt_all = dvector(nal);
  double sum_hwt_all = 0;

  for (j = 1; j <= alilen_mat; j++) {
    position_num_letters = 0;
    for (i = 1; i <= 20; i++) tmp_mark[i] = 0;
    for (i = 1; i <= nal; i++) {
      // cout << j << " " << i << " " << alignment[i][apos_filtr[j]] << endl;
      if (alignment[i][apos_filtr[j]] == 0) continue;  // ignore gaps
      if (tmp_mark[alignment[i][apos_filtr[j]]]) {
        tmp_mark[alignment[i][apos_filtr[j]]]++;
      } else {
        position_num_letters += 1;
        tmp_mark[alignment[i][apos_filtr[j]]]++;
      }
    }
    // cout << "position_num_letters: " << position_num_letters << endl;
    // for(i=1;i<=20;i++) cout << "tmp_mark " << i << " " << tmp_mark[i] <<
    // endl;
    for (i = 1; i <= nal; i++) {
      if (alignment[i][apos_filtr[j]] == 0) continue;  // ignore gaps
      hwt_all[i] +=
          1.0 / position_num_letters / tmp_mark[alignment[i][apos_filtr[j]]];
    }
  }
  for (i = 1; i <= nal; i++) {
    sum_hwt_all += hwt_all[i];
  }
  for (i = 1; i <= nal; i++) {
    hwt_all[i] /= sum_hwt_all;
  }
  // cout << "sum_hwt_all: " << sum_hwt_all << " alilen_mat: " << alilen_mat <<
  // endl;
}

void subalign::printProfile() {
  int i, j;
  cout << "nal: " << nal << endl;
  cout << "alilen: " << alilen << endl;
  cout << "gap_threshold: " << gap_threshold << endl;
  cout << "gapt_threshold: " << gapt_threshold << endl;
  cout << "alilen_mat: " << alilen_mat << endl;
  cout << "n_eff: " << n_eff << endl;

  for (i = 1; i <= alilen; i++) {
    cout << "maskgapRegion apos_filtr sum_eff_let " << i << " "
         << maskgapRegion[i] << " " << apos_filtr[i] << " " << sum_eff_let[i]
         << endl;
  }
}

// Construct a sub alignment with sequences marked in array mark[]
subalign *subalign::sub2align(int *mark) {
  int i, j, k;
  int nal1 = 0, alilen1;
  char **aseq1, **aname1;
  int **alignment1;
  int mlen;
  int count;

  subalign *aln1 = new subalign();

  for (i = 1; i <= nal; i++) {
    if (mark[i]) nal1++;
  }
  aln1->nal = nal1;
  aln1->alilen = alilen;

  mlen = 0;
  for (i = 0; i < nal; i++) {
    if (mark[i + 1])
      if (mlen < strlen(aname[i])) {
        mlen = strlen(aname[i]);
      }
  }
  aln1->mnamelen = mlen;

  aln1->aseq = cmatrix(nal1, alilen);
  aln1->aname = cmatrix(nal1, mnamelen);
  aln1->alignment = imatrix(nal1, alilen);

  count = 0;
  for (i = 1; i <= nal; i++) {
    if (mark[i]) {
      strcpy(aln1->aseq[count], aseq[i - 1]);
      strcpy(aln1->aname[count], aname[i - 1]);
      for (j = 1; j <= alilen; j++) {
        aln1->alignment[count + 1][j] = alignment[i][j];
      }
      count++;
    }
  }

  /* testing
  for(i=1;i<=nal1;i++)
      cout << aln1.aname[i-1] << "  " << aln1.aseq[i-1] << endl;
  for(i=1;i<=nal1;i++) {
      for(j=1;j<=aln1.alilen;j++) {
          fprintf(stdout, "%d ", aln1.alignment[i][j]);
      }
      fprintf(stdout, "\n");
  } */

  return aln1;
}

// Construct a sub alignment with sequences marked in array mark[1..nal]
// mark1[1..alilen]
subalign *subalign::sub2align(int *mark, int *mark1) {
  int i, j, k;
  int nal1 = 0, alilen1 = 0;
  char **aseq1, **aname1;
  int **alignment1;
  int mlen;
  int count, countl;

  subalign *aln1 = new subalign();

  // sequence number
  for (i = 1; i <= nal; i++) {
    if (mark[i]) nal1++;
  }
  aln1->nal = nal1;
  // sequence length
  for (i = 1; i <= alilen; i++) {
    if (mark1[i]) alilen1++;
  }
  aln1->alilen = alilen1;
  // maximum name length
  mlen = 0;
  for (i = 0; i < nal; i++) {
    if (mark[i + 1])
      if (mlen < strlen(aname[i])) {
        mlen = strlen(aname[i]);
      }
  }
  aln1->mnamelen = mlen;

  aln1->aseq = cmatrix(nal1, alilen1 + 1);
  aln1->aname = cmatrix(nal1, mlen);
  aln1->alignment = imatrix(nal1, alilen1);

  count = 0;
  for (i = 1; i <= nal; i++) {
    if (mark[i]) {
      countl = 0;
      for (j = 1; j <= alilen; j++) {
        if (mark1[j]) {
          aln1->aseq[count][countl] = aseq[i - 1][j - 1];
          // cout  <<"|"<< aseq1[count][countl] <<"|";
          aln1->alignment[count + 1][countl + 1] = alignment[i][j];
          countl++;
        }
      }
      // cout << endl;
      // cout << alilen1 << "\t" << countl << endl;
      aln1->aseq[count][aln1->alilen] = '\0';
      // if(mark1[j-1]) aseq1[count][countl]='\0';
      // else aseq1[count][countl+1]='\0';
      // cout << aname[i-1] << "\t" << aseq1[count] << "|" << endl;
      strcpy(aln1->aname[count], aname[i - 1]);
      count++;
    }
  }

  /* testing  //
  for(i=1;i<=nal1;i++)
      cout << aln1->aname[i-1] << "  " << aln1->aseq[i-1] << endl;
  //
  for(i=1;i<=nal1;i++) {
      for(j=1;j<=aln1->alilen;j++) {
          fprintf(stdout, "%d ", aln1->alignment[i][j]);
      }
      fprintf(stdout, "\n");
  }*/

  return aln1;
}

subalign *subalign::purge_align_one_seq_name(char *seqname) {
  int i, j;

  int target_index = -1;
  for (i = 0; i < nal; i++) {
    if (strcmp(seqname, aname[i]) == 0) {
      target_index = i;
      break;
    }
  }
  if (target_index == -1) {
    cout << "cannot find the target name in subalign " << seqname << endl;
    return NULL;
  }

  int *mark_ = ivector(nal);
  int *mark1_ = ivector(alilen);

  for (i = 1; i <= nal; i++) mark_[i] = 1;
  for (i = 1; i <= alilen; i++) {
    if (aseq[target_index][i - 1] == '-')
      mark1_[i] = 0;
    else
      mark1_[i] = 1;
  }

  return sub2align(mark_, mark1_);
}

// convert the alignment in letters to alignment in numbers
void subalign::convertAseq2Alignment() {
  int i, j;

  if (!alignment) alignment = imatrix(nal, alilen);

  for (i = 1; i <= nal; i++) {
    for (j = 1; j <= alilen; j++) {
      alignment[i][j] = am2num(aseq[i - 1][j - 1]);
    }
  }
}

void subalign::add_sequence(char *name, char *seq) {
  int i, j;

  if (nal == 0) {
    nal++;
    alilen = strlen(seq);
    mnamelen = strlen(name) + 1;

    aname = cmatrix(nal, mnamelen);
    aseq = cmatrix(nal, alilen);
    strcpy(aname[nal - 1], name);
    strcpy(aseq[nal - 1], seq);
    alignment = imatrix(nal, alilen);
    for (i = 1; i <= nal; i++) {
      for (j = 1; j <= alilen; j++) {
        alignment[i][j] = am2num(aseq[i - 1][j - 1]);
      }
    }
  }

  else {
    nal++;
    if (mnamelen < strlen(name) + 1) mnamelen = strlen(name) + 1;
    char **new_aname = cmatrix(nal, mnamelen);
    char **new_aseq = cmatrix(nal, alilen);
    for (i = 0; i < nal - 1; i++) {
      strcpy(new_aname[i], aname[i]);
      strncpy(new_aseq[i], aseq[i], alilen);
      new_aseq[i][alilen] = '\0';
      // cout << new_aname[i] << "\t"<<new_aseq[i]<<endl;
    }
    strcpy(new_aname[nal - 1], name);
    strncpy(new_aseq[nal - 1], seq, alilen);
    new_aseq[nal - 1][alilen] = '\0';
    // cout << new_aname[nal-1] << "\t"<<new_aseq[nal-1]<<endl;

    int **new_alignment = imatrix(nal, alilen);

    for (i = 1; i <= nal; i++) {
      for (j = 1; j <= alilen; j++) {
        new_alignment[i][j] = am2num(new_aseq[i - 1][j - 1]);
      }
    }

    for (i = 0; i < nal - 1; i++) {
      delete[] aseq[i];
      delete[] aname[i];
    }
    for (i = 0; i <= nal - 1; i++) delete[] alignment[i];
    delete[] aseq;
    delete[] aname;
    delete[] alignment;

    aseq = new_aseq;
    aname = new_aname;
    alignment = new_alignment;
    // cout << "===========" << endl;
  }
}

subalign *oneSeq2subalign(char *seq, char *name) {
  int i, j;

  subalign *a = new subalign();

  a->aseq = cmatrix(1, strlen(seq) + 1);
  a->aname = cmatrix(1, strlen(name) + 1);
  a->alignment = imatrix(1, strlen(seq));
  for (i = 0; i <= strlen(seq); i++) a->aseq[0][i] = seq[i];
  for (i = 0; i <= strlen(name); i++) a->aname[0][i] = name[i];
  for (i = 1; i <= strlen(seq); i++) a->alignment[1][i] = am2num(seq[i - 1]);
  a->nal = 1;
  a->alilen = strlen(seq);
  a->mnamelen = strlen(name);

  // cout << a->alilen << endl;

  return a;
}

// purge an alignment by
// 1. removing highly similar sequences (id > high_id_thr)
// 2. removing highly divergent sequences to the first sequence (id <
// low_id_thr)
// 3. removing sequence fragment relative to the first sequence (gap_fraction >
// gap_fraction_thr)
// 4. selecting only a subset if remaining num > max_num_kept
// *** all positions are kept ***
subalign *subalign::purge_align(double low_id_thr, double high_id_thr,
                                int max_num_kept, double gap_fraction_thr) {
  int i, j, k;
  int pairs, id_pairs;
  double gap_fraction;

  int max_selected = 10000;

  int *mark = ivector(nal);
  for (i = 1; i <= nal; i++) mark[i] = 1;

  assert(low_id_thr <= 1);
  assert(high_id_thr <= 1);

  for (i = 1; i < nal; i++) {
    if (i > max_selected) {
      mark[i + 1] = 0;
      continue;
    }
    if (strcmp(aname[i], aname[i - 1]) == 0) {
      mark[i + 1] = 0;
      continue;
    }
  }
  if (low_id_thr > 0)
    for (i = 1; i < nal; i++) {
      pairs = 0;
      id_pairs = 0;
      for (k = 0; k < alilen; k++) {
        if ((aseq[i][k] != '-') && (aseq[0][k] != '-')) {
          pairs += 1;
          if (aseq[i][k] == aseq[0][k]) id_pairs += 1;
        }
      }
      if (1.0 * id_pairs / pairs < low_id_thr) mark[i + 1] = 0;
      if (mark[i + 1] == 0) continue;

      pairs = 0;
      id_pairs = 0;
      for (k = 0; k < alilen; k++) {
        if (aseq[0][k] != '-') {
          pairs += 1;
          if (aseq[i][k] == '-') id_pairs += 1;
        }
      }
      gap_fraction = 1.0 * id_pairs / pairs;
      if (gap_fraction > gap_fraction_thr) mark[i + 1] = 0;
    }

  for (i = 0; i < nal; i++) {
    if (i > max_selected) {
      mark[i + 1] = 0;
      continue;
    }
    if (mark[i + 1] == 0) continue;
    for (j = i + 1; j < nal; j++) {
      if (mark[j + 1] == 0) {
        continue;
      }
      pairs = 0;
      id_pairs = 0;
      for (k = 0; k < alilen; k++) {
        if ((aseq[i][k] != '-') && (aseq[j][k] != '-')) {
          pairs += 1;
          if (aseq[i][k] == aseq[j][k]) id_pairs += 1;
        }
      }
      if (1.0 * id_pairs / pairs > high_id_thr) mark[j + 1] = 0;
    }
  }

  // for(i=1;i<=nal;i++) { if(mark[i]==1) tmpcount+=1; if(tmpcount>50) mark[i] =
  // 0; }

  int tmpcount = 0;
  for (i = 1; i <= nal; i++) {
    if (mark[i] == 1) tmpcount++;
    if (tmpcount > max_num_kept) mark[i] = 0;
  }

  if (mydebug > 100)
    for (i = 1; i <= nal; i++) cout << i << "\t" << mark[i] << endl;

  subalign *xaln = sub2align(mark);
  cout << "************" << endl;

  delete[] mark;
  return xaln;
}

void subalign::prof_h_weight_all(double raw_gap_fraction_cutoff) {
  int i, j, k;
  int position_num_letters;
  int tmp_mark[21];
  prof_hwt_all = dvector(nal);
  for (i = 1; i <= nal; i++) prof_hwt_all[i] = 1e-06;
  double sum_hwt_all = 0;
  int raw_gap_counts;

  for (j = 1; j <= alilen; j++) {
    raw_gap_counts = 0;
    for (i = 1; i <= nal; i++) {
      if (!alignment[i][j]) raw_gap_counts++;
    }
    if (raw_gap_counts * 1.0 / nal > raw_gap_fraction_cutoff) continue;
    // cout << "j: " << j << endl;
    position_num_letters = 0;
    for (i = 0; i <= 20; i++) tmp_mark[i] = 0;
    // *** IMPORTANT ***
    // gap character is not considered as the 21 letter
    // this is to avoid giving fragments unreasonably high weights
    for (i = 1; i <= nal; i++) {
      // cout << j << " " << i << " " << alignment[i][apos_filtr[j]] << endl;
      if (alignment[i][j] == 0) continue;  // ignore gaps
      if (tmp_mark[alignment[i][j]]) {
        tmp_mark[alignment[i][j]]++;
      } else {
        position_num_letters += 1;
        tmp_mark[alignment[i][j]]++;
      }
    }
    // cout << "position_num_letters: " << position_num_letters << endl;
    // for(i=1;i<=20;i++) cout << "tmp_mark " << i << " " << tmp_mark[i] <<
    // endl;
    for (i = 1; i <= nal; i++) {
      if (alignment[i][j] == 0) continue;  // ignore gaps
      prof_hwt_all[i] += 1.0 / position_num_letters / tmp_mark[alignment[i][j]];
    }
  }
  for (i = 1; i <= nal; i++) {
    sum_hwt_all += prof_hwt_all[i];
  }
  for (i = 1; i <= nal; i++) {
    prof_hwt_all[i] /= sum_hwt_all;
  }
  // cout << "sum_hwt_all: " << sum_hwt_all << " alilen_mat: " << alilen_mat <<
  // endl;

  if (mydebug > 1) {
    cout << "henikoff weights: " << endl;
    for (i = 1; i <= nal; i++) {
      cout << i << "\t" << prof_hwt_all[i] << endl;
    }
    cout << endl;
  }
}

void subalign::prof_positions(double prof_gap_threshold_here) {
  int i, j;

  prof_gap_content = dvector(alilen);
  prof_pos = ivector(alilen);

  prof_len = 0;
  for (i = 1; i <= alilen; i++) {
    prof_gap_content[i] = 0;
    for (j = 1; j <= nal; j++) {
      if (alignment[j][i] == 0) prof_gap_content[i] += prof_hwt_all[j];
    }
    // if gap frequency is prof_gap_threshold_here, still consider it to be core
    // position this is for a common value of 0.5 and when there are only two
    // sequences, then every position is considered to be core positions
    if (prof_gap_content[i] <= prof_gap_threshold_here) {
      prof_len++;
      prof_pos[prof_len] = i;
    }
  }

  if (prof_len == 0) {
    cout << "Warning: number of positions with weighted gap fraction less than "
         << prof_gap_threshold_here << " is zero" << endl;
    cout << "    - Use all the positions for profile" << endl;
    cout << "    - Each sequence has the same weight" << endl;
    cout << endl;
    for (i = 1; i <= alilen; i++) {
      prof_pos[i] = i;
    }
    prof_len = alilen;
  }

  if (mydebug > 1) {
    cout << "prof positions: " << endl;
    cout << "prof length: " << prof_len << endl;
    for (i = 1; i <= prof_len; i++) {
      cout << i << "\t" << prof_pos[i] << endl;
    }
    cout << endl;
  }
}

void subalign::set_prof_gap_threshold(double gapthr) {
  prof_gap_threshold = gapthr;
}

void subalign::set_prof_raw_gap_threshold(double gapthr) {
  prof_raw_gap_threshold = gapthr;
}

void subalign::prof() {
  int i, j, k;

  set_prof_raw_gap_threshold(0.5);
  set_prof_gap_threshold(0.5);

  // determine the henikoff weight; prof_position
  prof_h_weight_all(prof_raw_gap_threshold);
  prof_positions(prof_gap_threshold);

  prof_effn = dmatrix(prof_len, 20);
  prof_freq = dmatrix(prof_len, 20);
  prof_sum_eff = dvector(prof_len);

  prof_get_effn(alignment, nal, alilen, prof_len, prof_effn, prof_sum_eff,
                prof_pos, &prof_nef);

  // pseudoCounts(prof_effn, prof_nef, prof_len, prof_freq);

  done_prof = 1;
}

void subalign::prof(double tmp_gap_threshold) {
  int i, j, k;

  set_prof_raw_gap_threshold(0.5);
  set_prof_gap_threshold(tmp_gap_threshold);

  // determine the henikoff weight; prof_position
  prof_h_weight_all(prof_raw_gap_threshold);
  prof_positions(prof_gap_threshold);

  prof_effn = dmatrix(prof_len, 20);
  // prof_freq = dmatrix(prof_len, 20);
  prof_sum_eff = dvector(prof_len);

  prof_get_effn(alignment, nal, alilen, prof_len, prof_effn, prof_sum_eff,
                prof_pos, &prof_nef);

  // pseudoCounts(prof_effn, prof_nef, prof_len, prof_freq);

  done_prof = 1;
}

void subalign::get_prof_freq(int get_freq, int take_log) {
  int i, j, k;

  prof_freq = dmatrix(prof_len, 20);
  switch (get_freq) {
    case 0:
      pseudoCounts(prof_effn, prof_nef, prof_len, prof_freq);
      break;
      /*
      case 1: pseudoCounts(prof_effn, prof_nef, prof_len, prof_freq, rqmatrix0,
      rbfreq0); break; case 2: pseudoCounts(prof_effn, prof_nef, prof_len,
      prof_freq, rqmatrix, rbfreq, prof_alphabet1);
      */
      /*case 1: pseudoCounts(prof_effn, prof_nef, prof_len, prof_freq,
      aa_pair[0], aa_bg[0]); break; case 2: pseudoCounts(prof_effn, prof_nef,
      prof_len, prof_freq, aa_pair, aa_bg, ss->alphabet1);
      */
  }

  if (take_log) {
    for (i = 1; i <= prof_len; i++) {
      for (j = 1; j <= 20; j++) {
        prof_freq[i][j] = log(prof_freq[i][j]);
      }
    }
  }

  done_prof_freq = 1;
}

/* from a given alignment with aa as numbers, computes effective aa counts
(PSIC->our formula) and marks the columns with EFFECTIVE content of gaps >
threshold (effgapmax) */
void subalign::prof_get_effn(int **ali, int n, int len, int prof_len,
                             double **n_effAa, double *sum_eff_let,
                             int *prof_pos, double *nef) {
  int i, j, k, l;
  int nsymbols_col, nsymbols;
  double nef_loc;
  int ele;
  double effnu[21];
  double sum_let;
  int *mark;
  int flagmark;

  int p;

  gap_content = new double[alilen + 1];

  if (mydebug > 1) cout << "prof_len: " << prof_len << endl;

  mark = new int[n + 11];
  nsymbols = 0;
  for (j = 1; j <= prof_len; j++) {
    p = prof_pos[j];

    nsymbols_col = 0;
    sum_let = 0;
    if (mydebug > 100) cout << "p: " << p << endl;

    for (k = 0; k <= 20; ++k) {
      /* Mark sequences that have amino acid  k (or gap, k=0) in this jth
       * position */
      flagmark = 0;
      for (i = 1; i <= n; ++i) {
        mark[i] = 0;

        ele = ali[i][p];
        if (ele == k) {
          mark[i] = 1;
          flagmark = 1;
        }
        ele = ali[i][p] - 25;
        if (ele == k) {
          mark[i] = 1;
          flagmark = 1;
        }
      }

      /* If aa k (or gap) is present in this position call compute k-th
       * effective count */
      if (flagmark == 1) {
        effnu[k] = effective_number_nogaps(ali, mark, n, 1, len);
        nsymbols_col++;

      } else {
        effnu[k] = 0.0;
      }

      if (k > 0) sum_let += effnu[k];
    }

    // calculate raw gap fraction
    if (mydebug > 100) {
      int gap_count = 0;
      for (i = 1; i <= n; i++) {
        if (ali[i][p] == 0) gap_count++;
      }
      // double gapfract = 1.0 * gap_count / n;
      // cout << j << "  " << 1.0*effnu[0]/(sum_let + effnu[0]) << "  " <<
      // gapfract << endl;
    }

    if (mydebug > 100)
      cout << j << "  " << aseq[0][j - 1] << " " << sum_let << "  " << effnu[0]
           << endl;
    for (k = 0; k <= 20; k++) {
      n_effAa[j][k] = effnu[k];
    }
    if (mydebug > 100)
      cout << j << " " << nsymbols_col << "  " << nsymbols << endl;
    nsymbols += nsymbols_col;
    if (mydebug > 100) cout << "nsymbols: " << nsymbols << endl;

    sum_eff_let[j] = sum_let;

    /*
    if ( sum_let > 0 && 1.0*effnu[0]/(sum_let + effnu[0]) < effgapmax ) {
            alilen_mat++;
            gap_content[alilen_mat] = 1.0*effnu[0]/(sum_let + effnu[0]);
            for (k=0; k<=20; k++) {
                    n_effAa[alilen_mat][k] = effnu[k];
            }
            sum_eff_let[alilen_mat] = sum_let;
            apos_filtr[alilen_mat] = j;
            nsymbols += nsymbols_col;

            if(1.0*effnu[0]/(sum_let + effnu[0]) < effgapRegionMin) {
                     maskgapRegion[alilen_mat] = 0;
            } else {
                    maskgapRegion[alilen_mat] = 1;
            }

    }*/
  }

  // cout << nsymbols << "\t" << alilen_mat << endl;

  nef_loc = 1.0 * nsymbols / prof_len;
  *nef = nef_loc;

  prof_average_sum_eff_let = 0;
  for (i = 1; i <= prof_len; i++) {
    prof_average_sum_eff_let += sum_eff_let[i];
  }
  prof_average_sum_eff_let /= prof_len;

  delete[] mark;

  /* for(i=1;i<=alilen_mat;i++) {
      for(j=1;j<=20;j++) {
          cout << n_effAa[i][j] << " ";
      }
      cout << endl;
  } */
}

//
void subalign::get_prof_map_ss(char *rep_name) {
  int i, j, k;

  int index = -1;
  for (i = 0; i < nal; i++) {
    if (strcmp(rep_name, aname[i]) == 0) {
      index = i;
      break;
    }
  }
  if (index == -1) {
    cout << "get_prof_map_ss error" << endl;
    cout << "cannot find a sequence with the same name as " << rep_name << endl;
    exit(0);
  }
  if (mydebug > 1) cout << "index: " << index << endl;

  int *tmp_iarr = ivector(alilen);

  int count = 0;
  for (i = 0; i < alilen; i++) {
    if (aseq[index][i] != '-') {
      count++;
      tmp_iarr[i + 1] = count;
    } else
      tmp_iarr[i + 1] = 0;
    if (mydebug > 1) cout << i << " " << tmp_iarr[i + 1] << endl;
  }

  prof_map_ss = ivector(prof_len);

  for (i = 1; i <= prof_len; i++) {
    prof_map_ss[i] = tmp_iarr[prof_pos[i]];
    if (mydebug > 1)
      cout << i << " " << prof_pos[i] << " " << prof_map_ss[i] << endl;
  }
  // cout << "=============" << endl;

  delete[] tmp_iarr;
}

void subalign::get_prof_alphabet1() {
  int i, j, k;

  prof_alphabet1 = ivector(prof_len);

  for (i = 1; i <= prof_len; i++) {
    if (prof_map_ss[i] == 0)
      prof_alphabet1[i] = 0;
    else
      prof_alphabet1[i] = ss->alphabet1[prof_map_ss[i]];
  }
  // cout << "done_get_prof_alphabet1" << endl;
}

void subalign::select_representative() {
  int i, j, k;

  // find a representative sequence for the group, make it the "aln"
  // right now, the representative is the longest sequence (excluding gaps)
  int tmp_index = 0;
  int tmp_count_aa = 0;
  int max_count_aa = 0;
  for (i = 0; i < nal; i++) {
    // if(strlen(aname[i])<=6) { tmp_index = i; break; }
    tmp_count_aa = 0;
    for (j = 0; j < alilen; j++) {
      if (aseq[i][j] != '-') tmp_count_aa++;
    }
    if (tmp_count_aa > max_count_aa) {
      max_count_aa = tmp_count_aa;
      tmp_index = i;
    }
  }
  /*
  cout << "nal: " << nal << endl;
  for(i=0;i<nal;i++) {
          cout << "|" << aname[i] << "|: " << strlen(aname[i]) << endl;
  }
  */
  char *tmp_name = new char[strlen(aname[tmp_index]) + 1];
  char *tmp_seq = new char[max_count_aa + 1];
  strcpy(tmp_name, aname[tmp_index]);
  int tmp_array_index = 0;
  for (i = 0; i < alilen; i++) {
    if (aseq[tmp_index][i] != '-') {
      tmp_seq[tmp_array_index] = aseq[tmp_index][i];
      tmp_array_index++;
    }
  }
  tmp_seq[tmp_array_index] = '\0';
  oneSeqAln = oneSeq2subalign(tmp_seq, tmp_name);
  // cout << "representative: " << tmp_name << endl;

  repres_name = new char[strlen(oneSeqAln->aname[0]) + 2];
  strcpy(repres_name, oneSeqAln->aname[0]);

  // return oneSeqAln;
}

// select a representative sequence with one of the best
// sum-of-henikoff-frequencies and minimum number of residues in the gappy
// region
int subalign::select_representative_henikoff(
    float core_position_nongap_freq_cutoff) {
  int i, j, k;

  // debug
  if (mydebug > 1100) printali(80);

  // 1. get the henikoff weighting
  // determine the henikoff weight; prof_position
  set_prof_raw_gap_threshold(0.5);
  set_prof_gap_threshold(core_position_nongap_freq_cutoff);
  prof_h_weight_all(prof_raw_gap_threshold);
  // for(i=1;i<=nal;i++) { cout << "h_weight: " << i << " " << prof_hwt_all[i]
  // << endl; }
  prof_positions(prof_gap_threshold);

  // 1.1 get the prof_pos1
  // prof_pos stores the indexes of ungappy regions in a consecutive way
  // prof_pos1 marks the ungappy positions with 1, and gappy positions with 0
  int *prof_pos1 = ivector(alilen);
  for (i = 1; i <= alilen; i++) {
    prof_pos1[i] = 0;
  }
  for (i = 1; i <= alilen; i++) {
    if (prof_pos[i] != 0) {
      prof_pos1[prof_pos[i]] = 1;
    }
  }

  // 2. get henikoff frequencies in each position, gap is treated as the 21st
  // letter
  float **hfreq = gmatrix<float>(alilen, 20);
  for (i = 1; i <= alilen; i++)
    for (j = 0; j <= 20; j++) hfreq[i][j] = 0;
  // float sum_freq = 0;
  for (i = 1; i <= alilen; i++) {
    // ignore the gappy regions, all the frequencies are 0 in these positions
    if (prof_pos1[i] == 0) continue;
    // sum_freq = 0;
    for (j = 1; j <= nal; j++) {
      hfreq[i][alignment[j][i]] += prof_hwt_all[j];
    }
  }

  // 3. select the sequence with one of the largest sum of hfreq across the
  // positions
  //      - determine the sum of freq for each sequence, find the maximum value
  //      (max_sum_freq)
  //      - consider those sequences that has a sum_of_freq > 0.95 *
  //      max_sum_freq
  //      - among these sequences, select the one with the lowest number of
  //      amino acids in gappy regions
  //      - if numbers of amino acids in gappy regions the same, select the one
  //      with the highest sumfreqs
  int tmp_index = -1;
  float sum_freq = 0;
  float max_sum_freq = 0;
  float *sumfreqs = gvector<float>(nal);
  int *nlet_gappy = ivector(nal);
  for (i = 1; i <= nal; i++) {
    // 3.1 determine sumfreqs
    sum_freq = 0;
    for (j = 1; j <= alilen; j++) {
      // only for amino acid frequences, ignore gaps
      if (alignment[i][j] == 0) continue;
      sum_freq += hfreq[j][alignment[i][j]];
    }
    sumfreqs[i] = sum_freq;
    if (sum_freq > max_sum_freq) {
      tmp_index = i;
      max_sum_freq = sum_freq;
    }
    // 3.2 determine number of letters in gappy regions for each sequence
    nlet_gappy[i] = 0;
    for (j = 1; j <= alilen; j++) {
      if (!prof_pos1[j])
        if (alignment[i][j]) nlet_gappy[i]++;
    }
    // debug *******************
    if (mydebug > 1100)
      cout << i << " " << aname[i - 1] << " " << sum_freq << " "
           << nlet_gappy[i] << endl;
  }
  assert(tmp_index > 0);
  assert(tmp_index <= nal);
  float max_sum_freq1 = 0;  // the current best among those valid positions
  int min_nlet_gappy = 1000000;
  for (i = 1; i <= nal; i++) {
    if (sumfreqs[i] < max_sum_freq * 0.9) continue;
    if (nlet_gappy[i] < min_nlet_gappy) {
      min_nlet_gappy = nlet_gappy[i];
      tmp_index = i;
      max_sum_freq1 = sumfreqs[i];
    } else if (nlet_gappy[i] == min_nlet_gappy) {
      if (sumfreqs[i] > max_sum_freq1) {
        tmp_index = i;
        max_sum_freq1 = sumfreqs[i];
      }
    }
  }
  assert(tmp_index <= nal);
  delete[] nlet_gappy;
  delete[] sumfreqs;
  tmp_index = tmp_index - 1;  // conform to the index of aname and aseq

  // 4. modify the representative sequence, so that core positions with gaps in
  // the
  //    representative is replaced by letters of the largest frequency
  float max_hfreq = 0;
  int max_index = 0;
  string tmpstr = "";
  for (i = 0; i < alilen; i++) {
    if (prof_pos1[i + 1] == 0) continue;
    if (aseq[tmp_index][i] == '-') {
      max_hfreq = 0;
      max_index = 0;
      for (j = 1; j <= 20; j++) {
        if (hfreq[i + 1][j] > max_hfreq) {
          max_hfreq = hfreq[i + 1][j];
          max_index = j;
        }
      }
      aseq[tmp_index][i] = am[max_index] + 32;
      tmpstr += aseq[tmp_index][i];
    }
  }

  // debug *******************
  if (mydebug > 1100) {
    cout << endl;
    cout << "henikoff weights" << endl;
    for (i = 1; i <= nal; i++) {
      cout << i << " " << prof_hwt_all[i] << endl;
    }
    cout << endl;
  }
  cout << "representative: " << aname[tmp_index] << endl;
  if (tmpstr.size()) cout << "filled segment: " << tmpstr << endl;
  for (i = 1; i <= alilen; i++) {
    cout << prof_pos1[i];
  }
  cout << endl;
  cout << aseq[tmp_index] << endl;
  if (mydebug > 1100) {
    printali(80);
    cout << "===============" << endl;
  }

  // 4. clear up the memory
  delete[] prof_hwt_all;
  delete[] prof_pos;
  delete[] prof_pos1;
  free_gmatrix<float>(hfreq, alilen, 20);

  return tmp_index;
}

void subalign::get_oneSeqAln(int tmp_index) {
  int i, j, k;

  // 5. store the sequence and get its subalign in oneSeqAln
  char *tmp_name = new char[strlen(aname[tmp_index]) + 1];
  int count_aa = 0;
  for (i = 1; i <= alilen; i++) {
    if (aseq[tmp_index][i] != '-') count_aa++;
  }
  char *tmp_seq = new char[count_aa + 1];
  strcpy(tmp_name, aname[tmp_index]);
  int tmp_array_index = 0;
  for (i = 0; i < alilen; i++) {
    if (aseq[tmp_index][i] != '-') {
      tmp_seq[tmp_array_index] = aseq[tmp_index][i];
      tmp_array_index++;
    }
  }
  tmp_seq[tmp_array_index] = '\0';
  oneSeqAln = oneSeq2subalign(tmp_seq, tmp_name);
  cout << "representative: " << tmp_name << endl;

  repres_name = new char[strlen(oneSeqAln->aname[0]) + 2];
  strcpy(repres_name, oneSeqAln->aname[0]);

  // return oneSeqAln;
}

// for representative sequence with repres_name
// if base_name.ss2 or base_name.horiz not in directory dir_name:
// make a directory in dir_name
// go to the new directory
// run runpsipred
// read psipred output files
// copy psipred output files to dir_name
void subalign::get_ss_prof(char *dir_name, char *runpsipred_command) {
  int i;
  char base_name[300];
  char cwd[300];

  if (strlen(dir_name) == 0) {
    strcpy(base_name, repres_name);
  } else {
    strcpy(base_name, dir_name);
    if (dir_name[strlen(dir_name) - 1] != '/') {
      strcat(base_name, "/");
    }
    strcat(base_name, repres_name);
  }

  // cout << base_name << " " << repres_name << endl;

  // read the psipred files if they exist
  ss = new ss_prof(base_name);
  if (ss->done_prof) {
    // check if the sequence in psipred output files is the same as the target
    // sequence
    if (ss->check_aa_seq(oneSeqAln->aseq[0]) == 1) return;
  }

  // set a temporary directory to run psipred:
  // dir_name/represname_randnum_psipred
  int rand_number = rand();
  char tmp_dir[500];
  sprintf(tmp_dir, "%s/%s_%d_psipred", dir_name, repres_name, rand_number);
  char command[500];
  sprintf(command, "mkdir %s", tmp_dir);
  system(command);

  // generate a fasta file in the temporary directory
  char fasta_name[500];
  sprintf(fasta_name, "%s/%s.fa", tmp_dir, repres_name);
  ofstream ofp(fasta_name, ios::out);
  ofp << ">";
  ofp << repres_name << endl;
  // here, change '[OJoj]' to 'A'
  for (i = 0; i < oneSeqAln->alilen; i++) {
    if ((oneSeqAln->aseq[0][i] == 'O') || (oneSeqAln->aseq[0][i] == 'o') ||
        (oneSeqAln->aseq[0][i] == 'U') || (oneSeqAln->aseq[0][i] == 'u')) {
      ofp << 'A';
    } else
      ofp << oneSeqAln->aseq[0][i];
  }
  ofp << endl;
  // old version: ofp << oneSeqAln->aseq[0] << endl;
  ofp.close();

  // remember the current directory
  getcwd(cwd, 300);
  if (mydebug > 1) cout << "cwd: " << cwd << endl;

  // go to temporary directory; run psipred
  chdir(tmp_dir);
  // sprintf(command, "runpsipred %s.fa", repres_name);
  sprintf(command, "%s %s.fa", runpsipred_command, repres_name);
  cout << command << endl;
  system(command);

  // get the psipred profile object
  ss->read_ss_files(repres_name);

  // copy the psipred files to psipred_dir
  sprintf(command, "cp %s.ss2 ../", repres_name);
  system(command);
  sprintf(command, "cp %s.horiz ../", repres_name);
  system(command);

  // clear the temporary directory
  chdir("../");
  // sprintf(command, "rm -r -f %s", tmp_dir);
  // system(command);

  // change directory back
  chdir(cwd);
}

// for representative sequence with repres_name
// if base_name.ss2 or base_name.horiz not in directory dir_name:
// make a directory in dir_name
// go to the new directory
// run runpsipred1 (assume psipmd already exist)
// read psipred output files
// copy psipred output files to dir_name
void subalign::get_ss_prof1(char *dir_name, char *query_name,
                            char *runpsipred1_command) {
  int i;
  char base_name[300];
  char cwd[300];

  if (strlen(dir_name) == 0) {
    strcpy(base_name, repres_name);
  } else {
    strcpy(base_name, dir_name);
    if (dir_name[strlen(dir_name) - 1] != '/') {
      strcat(base_name, "/");
    }
    strcat(base_name, repres_name);
  }

  // cout << base_name << " " << repres_name << endl;

  // read the psipred files if they exist
  ss = new ss_prof(base_name);
  if (ss->done_prof) {
    // check if the sequence in psipred output files is the same as the target
    // sequence
    if (ss->check_aa_seq(oneSeqAln->aseq[0]) == 1) return;
  }

  // set a temporary directory to run psipred:
  // dir_name/represname_randnum_psipred
  int rand_number = rand();
  char tmp_dir[500];
  if (dir_name[strlen(dir_name) - 1] != '/') {
    sprintf(tmp_dir, "%s/%s_%d_psipred", dir_name, repres_name, rand_number);
  } else {
    sprintf(tmp_dir, "%s%s_%d_psipred", dir_name, repres_name, rand_number);
  }
  char command[500];
  sprintf(command, "mkdir %s", tmp_dir);
  int a1 = system(command);
  cout << "a1: " << a1 << endl;

  // generate a fasta file in the temporary directory
  char fasta_name[500];
  sprintf(fasta_name, "%s/%s.fa", tmp_dir, repres_name);
  ofstream ofp(fasta_name, ios::out);
  ofp << ">";
  ofp << repres_name << endl;
  // here, change '[OJoj]' to 'A'
  for (i = 0; i < oneSeqAln->alilen; i++) {
    if ((oneSeqAln->aseq[0][i] == 'O') || (oneSeqAln->aseq[0][i] == 'o') ||
        (oneSeqAln->aseq[0][i] == 'U') || (oneSeqAln->aseq[0][i] == 'u')) {
      ofp << 'A';
    } else
      ofp << oneSeqAln->aseq[0][i];
  }
  ofp << endl;
  // old version: ofp << oneSeqAln->aseq[0] << endl;
  ofp.close();

  // remember the current directory
  char *status1 = getcwd(cwd, 300);
  if (mydebug > 1) cout << "cwd: " << cwd << endl;
  cout << "1: " << status1 << endl;

  // go to temporary directory; run psipred
  int status = chdir(tmp_dir);
  cout << "2: " << status << endl;
  sprintf(command, "cp ../%s.chk psitmp.chk", query_name);
  status = system(command);
  cout << "3: " << status << endl;
  status = system(command);
  cout << "4: " << status << endl;

  // here check if the sequence contains small letters
  run_psipred_check_lowercase_letters(repres_name, runpsipred1_command);
  // sprintf(command, "runpsipred %s.fa", repres_name);
  // sprintf(command, "%s %s.fa 2>/dev/null 1>/dev/null", runpsipred1_command,
  // repres_name); cout << command << endl; system(command);

  // get the psipred profile object
  ss->read_ss_files(repres_name);

  // copy the psipred files to psipred_dir
  sprintf(command, "cp %s.ss2 ../", repres_name);
  int ss2_status = system(command);
  if (ss2_status) {
    cout << "Error: .ss2 file not generated for " << query_name << endl;
    exit(1);
  }
  sprintf(command, "cp %s.horiz ../", repres_name);
  system(command);

  // clear the temporary directory
  // chdir("../../");
  chdir(cwd);
  // sprintf(command, "rm -r -f %s", tmp_dir);
  if (mydebug > 1) system("pwd");
  if (mydebug > 1) cout << "command: " << command << endl;
  // system(command);

  // change directory back
  chdir(cwd);
}

void subalign::checkSubalign() {
  int i, j, k;

  cout << "alilen: " << alilen << endl;
  cout << "nal: " << nal << endl;

  if (done_prof == 0) return;

  cout << "prof_len: " << prof_len << endl;
  cout << "prof_nef: " << prof_nef << endl;
  cout << "prof_average_sum_eff_let: " << prof_average_sum_eff_let << endl;
  cout << "prof_gap_threshold: " << prof_gap_threshold << endl;
  cout << "prof_raw_gap_threshold: " << prof_raw_gap_threshold << endl;
  cout << "num     prof_pos     prof_sum_eff   prof_gap_content" << endl;
  for (i = 1; i <= prof_len; i++) {
    cout << i << " " << prof_pos[i] << " " << aseq[0][prof_pos[i] - 1] << " "
         << prof_sum_eff[i] << " " << prof_gap_content[prof_pos[i]] << endl;
  }

  cout << "nal " << nal << endl;
  cout << "prof_hwt_all " << endl;
  for (i = 1; i <= nal; i++) {
    cout << i << " " << prof_hwt_all[i] << endl;
  }

  cout << endl;
  cout << "num prof_effn prof_freq" << endl;
  for (i = 1; i <= prof_len; i++) {
    cout << i << endl;
    for (j = 1; j <= 20; j++) {
      cout << am[j] << " " << setw(20) << prof_effn[i][j] << setw(20)
           << prof_freq[i][j] << endl;
    }
    cout << endl;
  }

  cout << endl;
}

void subalign::checkSSprof() {
  int i, j, k;

  cout << "secondary structure profile" << endl;
  cout << "repres_name:" << repres_name << endl;
  cout << "oneSeqAln:" << endl;
  oneSeqAln->printali(60);
  cout << "profpos prof_alphabet1   prof_map_ss:" << endl;
  for (i = 1; i <= prof_len; i++) {
    cout << i << " " << prof_alphabet1[i] << " " << prof_map_ss[i] << endl;
  }
}

void subalign::get_score_bg(int bg_type) {
  int i, j, k;

  bg_type_here = bg_type;

  score_bg_aa = gvector<float>(prof_len);
  score_bg_ss = gvector<float>(prof_len);
  for (i = 1; i <= prof_len; i++) {
    score_bg_aa[i] = 0;
    for (j = 1; j <= 20; j++) {
      if (bg_type == 0) {
        score_bg_aa[i] += prof_effn[i][j] * log_robinson_freq[j];
      } else if (bg_type == 1) {
        score_bg_aa[i] += prof_effn[i][j] * log_rbfreq0[j];
      } else if (bg_type == 2) {
        score_bg_aa[i] += prof_effn[i][j] * log_rbfreq[prof_alphabet1[i]][j];
      }
      if (mydebug > 1)
        cout << j << " " << prof_effn[i][j] << " " << log_robinson_freq[j]
             << " " << log_rbfreq0[j] << " " << log_rbfreq[prof_alphabet1[i]][j]
             << endl;
    }
    if (!ss) continue;
    score_bg_ss[i] = 0;
    for (j = 1; j <= 3; j++) {
      // ss->ssfreq[prof_alphabet1[i]][j] could be wrong!
      score_bg_ss[i] +=
          prof_sum_eff[i] * ss->ssfreq[prof_alphabet1[i]][j] * log_rssfreq[j];
      if (mydebug > 1)
        cout << j << " " << prof_sum_eff[i] << " "
             << ss->ssfreq[prof_alphabet1[i]][j] << " " << log_rssfreq[j]
             << endl;
    }
  }
}

// my own emission probability of a single residue
// e(i, a) = sum(q(i, a) * f(i, a))
// basically the inner product of profile frequency and background frequency
// in other words, the expected emission probability of a residue
// this is a bad scoring function; background amino acid emission frequencies
// should not depend on environment
void subalign::get_score_bg_mine(float **aa_loop, int bg_type) {
  int i, j, k;

  bg_type_here = bg_type;

  int *ss_alphabet;

  if (bg_type == 3)
    ss_alphabet = ss->sstype;
  else if (bg_type == 9)
    ss_alphabet = ss->alphabet1;

  score_bg_aa = gvector<float>(prof_len);
  score_bg_ss = gvector<float>(prof_len);
  if (mydebug > 1) cout << "test here: " << endl;
  for (i = 1; i <= prof_len; i++) {
    score_bg_aa[i] = 0;
    for (j = 1; j <= 20; j++) {
      score_bg_aa[i] += prof_freq[i][j] * aa_loop[ss_alphabet[i]][j];
      if (mydebug > 1)
        cout << i << " " << j << " " << prof_freq[i][j] << " "
             << aa_loop[ss_alphabet[i]][j] << endl;
    }
    score_bg_aa[i] = log(score_bg_aa[i]);
    score_bg_ss[i] = 0;
    // for(j=1;j<=3;j++) { score_bg_ss[i] += ss->ssfreq[i][j] *
    // ss_loop[ss_alphabet[i]][j]; }
    score_bg_ss[i] = ss->ssfreq[i][ss->sstype[i]];
    score_bg_ss[i] = log(score_bg_ss[i]);
  }
  done_score_bg = 1;
}

// my own emission probability of a single residue
// e(i, a) = sum(q(i, a) * f(i, a))
// basically the inner product of profile frequency and background frequency
// in other words, the expected emission probability of a residue
void subalign::get_score_bg_mine(float *aa_loop0, float *ss_loop0,
                                 int bg_type) {
  int i, j, k;

  bg_type_here = bg_type;

  int *ss_alphabet;

  if (bg_type == 3)
    ss_alphabet = ss->sstype;
  else if (bg_type == 9)
    ss_alphabet = ss->alphabet1;

  score_bg_aa = gvector<float>(prof_len);
  score_bg_ss = gvector<float>(prof_len);
  if (mydebug > 1) cout << "test here: " << endl;
  for (i = 1; i <= prof_len; i++) {
    score_bg_aa[i] = 0;
    // IMPORTANT: use raw frequencies derived from effn, without mixture with
    // pseudocount freqs expected emission probability given a effn vector
    for (j = 1; j <= 20; j++) {
      score_bg_aa[i] +=
          prof_effn[i][j] *
          aa_loop0[j];  // inner product of effn and loop frequencies
    }
    score_bg_aa[i] =
        log(score_bg_aa[i] / prof_sum_eff[i]);  // divide by sum effn
    score_bg_ss[i] = log(ss_loop0[ss_alphabet[i]]);
  }
  done_score_bg = 1;
}

void subalign::get_score_bg_mine2(float *aa_loop0, float *ss_loop0,
                                  int bg_type) {
  int i, j, k;

  bg_type_here = bg_type;

  int *ss_alphabet;
  if (bg_type == 3)
    ss_alphabet = ss->sstype;
  else if (bg_type == 9)
    ss_alphabet = ss->alphabet1;

  score_bg_aa = gvector<float>(prof_len);
  score_bg_ss = gvector<float>(prof_len);
  for (i = 1; i <= prof_len; i++) {
    score_bg_aa[i] = 0;
    for (j = 1; j <= 20; j++) {
      // IMPORTANT: aa_loop[i][j] is the logarithm of loop frequencies
      score_bg_aa[i] += aa_loop0[j] * prof_effn[i][j];
    }
    score_bg_aa[i] /= prof_sum_eff[i];
    score_bg_ss[i] = log(ss_loop0[ss_alphabet[i]]);
  }

  done_score_bg2 = 1;
}

void subalign::add_NC_terminal(char *query_name, char *query_seq) {
  int i, j, k;

  char *tmpstring = new char[alilen + 1];
  j = 0;
  for (i = 0; i < alilen; i++) {
    if (aseq[0][i] != '-') {
      tmpstring[j] = aseq[0][i];
      j++;
    }
  }
  tmpstring[j] = '\0';

  int right_numbers, left_numbers;

  right_numbers = string(query_seq).find(tmpstring, 0);
  if (right_numbers == string::npos) {
    cout << "Warning subalign: cannot find the sequence in " << query_name
         << endl;
    return;
  }
  left_numbers = strlen(query_seq) - right_numbers - strlen(tmpstring);

  assert(right_numbers >= 0);
  assert(left_numbers >= 0);

  // the first sequence of blast matches the query sequence exactly
  if ((right_numbers + left_numbers) == 0) return;

  string Nterm, Cterm, Nterm_gap, Cterm_gap;
  Nterm = string(query_seq).substr(0, right_numbers);
  Cterm =
      string(query_seq).substr(right_numbers + strlen(tmpstring), left_numbers);
  Nterm_gap = string("");
  for (i = 0; i < Nterm.size(); i++) Nterm_gap += "-";
  Cterm_gap = string("");
  for (i = 0; i < Cterm.size(); i++) Cterm_gap += "-";

  char **tmpseq = cmatrix(nal, alilen + right_numbers + left_numbers + 1);
  for (i = 0; i < nal; i++) {
    if (i == 0) {
      strcpy(tmpseq[i], Nterm.c_str());
      strcat(tmpseq[i], aseq[i]);
    } else {
      strcpy(tmpseq[i], Nterm_gap.c_str());
      strcat(tmpseq[i], aseq[i]);
    }
  }
  // cout << "nstr: " << Nterm_gap.c_str() << endl;
  // cout << "nstr: " << Nterm.c_str() << endl;
  if (left_numbers != 0) {
    // cout << "cstr: " <<  Cterm_gap.c_str() << endl;
    // cout << "cstr: " <<  Cterm.c_str() << endl;
    for (i = 0; i < nal; i++) {
      if (i == 0) {
        strcat(tmpseq[i], Cterm.c_str());
        tmpseq[i][alilen + right_numbers + left_numbers] = '\0';
      } else {
        strcat(tmpseq[i], Cterm_gap.c_str());
        tmpseq[i][alilen + right_numbers + left_numbers] = '\0';
      }
    }
  }
  for (i = 0; i < nal; i++) {
    delete[] aseq[i];
    aseq[i] = tmpseq[i];
    // cout << aseq[i] << endl;
  }

  // cout << "alilen: " << alilen << endl;
  alilen = alilen + right_numbers + left_numbers;
  // cout << "alilen: " << alilen << endl;
  for (i = 1; i <= nal; i++) {
    delete[] alignment[i];
    alignment[i] = new int[alilen + 1];
    for (j = 1; j <= alilen; j++) {
      alignment[i][j] = am2num(aseq[i - 1][j - 1]);
    }
  }
}

// assume a directory contains .fa file and psitmp.chk file
// the record in .fa file matches that in the psitmp.chk file
// however, .fa file sequence contains lowercase letters, which we would like to
// ignore for ss prediction create a new .fa file and a new psitmp.chk file, in
// which the small letters are removed run psipred command to get the results,
// and insert small letters into the result file, the inserted secondary
// structure types are 'C' and confidence level the lowest, equal probability
void run_psipred_check_lowercase_letters(char *repres_name,
                                         char *runpsipred1_command) {
  int i, j, k, l;
  char origseq[10000], newseq[10000];
  char command[1000];
  int psipred_status;

  // 1. read the fasta file, assume it contains only two lines: defline and
  // sequence
  char fastafile[500];
  sprintf(fastafile, "%s.fa", repres_name);
  ifstream fp(fastafile, ios::in);
  fp.getline(origseq, 10000);
  fp.getline(origseq, 10000);
  fp.close();

  // 2. check for lowercase letters, get newseq and lowercase_mark, a 0 and 1
  // array
  int origseqlen = strlen(origseq);
  int lowercase_mark[origseqlen];
  char *ss = newseq;
  for (i = 0; i < strlen(origseq); i++) {
    lowercase_mark[i] = 1;
    if (origseq[i] >= 'A')
      if (origseq[i] <= 'Z') {
        *ss = origseq[i];
        ss++;
        lowercase_mark[i] = 0;
      }
  }
  *ss = '\0';

  // if no small letters found, just run the command and return
  if (strlen(newseq) == origseqlen) {
    sprintf(command, "%s %s.fa 2>/dev/null 1>/dev/null", runpsipred1_command,
            repres_name);
    if (debug > -1) cout << command << endl;
    psipred_status = system(command);
    cout << "psipred_status: " << psipred_status << endl;
    if (psipred_status != 0) {
      cout << "Error running psipred for " << repres_name << endl;
      exit(1);
    }
    return;
  }

  // 3. read the psitmp.chk file, write it to psitmp.chk1
  FILE *pfp, *opfp;
  if ((pfp = fopen("psitmp.chk", "rb")) == NULL) {
    cout << "cannot read the psitmp.chk file " << repres_name << endl;
    exit(1);
  }
  if ((opfp = fopen("psitmp.chk1", "wb")) == NULL) {
    cout << "cannot read the psitmp.chk1 file for writing" << repres_name
         << endl;
    exit(1);
  }
  int chklen;
  fread(&chklen, sizeof(int), 1, pfp);
  assert(chklen == origseqlen);  // check lengthes match
  char chkseq[chklen + 1];
  fread(chkseq, sizeof(char), chklen, pfp);
  chkseq[chklen] = '\0';
  if (debug > 1) cout << "chkseq: " << chkseq << endl;
  int newseqlen = strlen(newseq);
  fwrite(&newseqlen, sizeof(int), 1, opfp);
  fwrite(&newseq, sizeof(char), newseqlen, opfp);
  double aa_prof[20];
  for (i = 1; i <= chklen; i++) {
    fread(aa_prof, sizeof(double), 20, pfp);
    if (debug > 1) cout << chkseq[i - 1] << " ";
    if (debug > 1)
      for (j = 0; j < 20; j++) {
        fprintf(stdout, "%4.3f ", aa_prof[j]);
      }
    if (debug > 1) cout << endl;
    if (lowercase_mark[i - 1] == 1) continue;
    fwrite(aa_prof, sizeof(double), 20, opfp);
  }
  fclose(pfp);
  fclose(opfp);

  // 4. update .fa and .chk files
  ofstream fafp(fastafile, ios::out);
  fafp << ">" << repres_name << "\n";
  fafp << newseq << "\n";
  fafp.close();
  system("cp psitmp.chk1 psitmp.chk");
  system("cp psitmp.chk1 tmp.chk1");

  // 5. run psipred
  sprintf(command, "%s %s.fa 2>/dev/null 1>/dev/null", runpsipred1_command,
          repres_name);
  if (debug > 1) cout << command << endl;
  psipred_status = system(command);
  cout << "psipred_status: " << psipred_status << endl;
  if (psipred_status != 0) {
    cout << "Error running psipred for " << repres_name << endl;
    exit(1);
  }

  // 6. update the .ss2 file and .horiz file
  sprintf(command, "mv %s.ss2 %s.ss2_1", repres_name, repres_name);
  system(command);
  sprintf(command, "mv %s.horiz %s.horiz_1", repres_name, repres_name);
  system(command);
  FILE *ssfp_1, *ssfp, *horizfp_1, *horizfp;
  if ((ssfp = fopen((string(repres_name) + ".ss2").c_str(), "w")) == NULL) {
    cout << "canoot open the .ss2 file for write " << repres_name << endl;
    exit(0);
  }
  if ((ssfp_1 = fopen((string(repres_name) + ".ss2_1").c_str(), "r")) == NULL) {
    cout << "canoot open the .ss2_1 file for read " << repres_name << endl;
    exit(0);
  }
  char ss2line[200];
  while (!feof(ssfp_1)) {
    fgets(ss2line, 200, ssfp_1);
    if (strlen(ss2line) == 1) {
      fputs(ss2line, ssfp);
      break;
    }
    if (ss2line[0] == '#') {
      fputs(ss2line, ssfp);
      continue;
    }
  }
  for (i = 0; i < origseqlen; i++) {
    if (lowercase_mark[i] == 0) {
      fgets(ss2line, 200, ssfp_1);
      for (ss = ss2line; *ss == ' '; ss++)
        ;
      for (ss = ss; isdigit(*ss); ss++)
        ;
      sprintf(ss2line, "%4d%s", i + 1, ss);
      fputs(ss2line, ssfp);
    } else {
      sprintf(ss2line, "%4d%2c C   0.333  0.333  0.334\n", i + 1,
              toupper(origseq[i]));
      fputs(ss2line, ssfp);
    }
  }
  fclose(ssfp);
  fclose(ssfp_1);

  if ((horizfp = fopen((string(repres_name) + ".horiz").c_str(), "w")) ==
      NULL) {
    cout << "canoot open the .horiz file for write " << repres_name << endl;
    exit(0);
  }
  if ((horizfp_1 = fopen((string(repres_name) + ".horiz_1").c_str(), "r")) ==
      NULL) {
    cout << "canoot open the .horiz_1 file for read " << repres_name << endl;
    exit(0);
  }
  char confline[10000], predline[10000], aaline[10000];
  strcpy(confline, "");
  strcpy(aaline, "");
  strcpy(predline, "");
  while (!feof(horizfp_1)) {
    fgets(ss2line, 200, horizfp_1);
    if (ss2line[strlen(ss2line) - 1] == '\n')
      ss2line[strlen(ss2line) - 1] = '\0';
    if (strncmp(ss2line, "Conf:", 5) == 0) {
      strcat(confline, ss2line + 6);
    }
    if (strncmp(ss2line, "Pred:", 5) == 0) {
      strcat(predline, ss2line + 6);
    }
    if (strncmp(ss2line, "  AA:", 5) == 0) {
      strcat(aaline, ss2line + 6);
    }
  }
  if (debug > 1) cout << confline << endl;
  if (debug > 1) cout << newseq << endl;
  assert(strlen(confline) == newseqlen);
  assert(strlen(predline) == newseqlen);
  assert(strlen(aaline) == newseqlen);
  char confline1[10000], predline1[10000], aaline1[10000];
  strcpy(confline1, "");
  strcpy(aaline1, "");
  strcpy(predline1, "");
  int count = 0;
  for (i = 0; i < origseqlen; i++) {
    if (lowercase_mark[i]) {
      confline1[i] = '0';
      predline1[i] = 'C';
      aaline1[i] = toupper(origseq[i]);
    } else {
      confline1[i] = confline[count];
      predline1[i] = predline[count];
      aaline1[i] = aaline[count];
      count++;
    }
  }
  confline1[i] = '\0';
  aaline1[i] = '\0';
  predline[i] = '\0';
  fclose(horizfp_1);

  fputs("# PSIPRED HFORMAT (PSIPRED V2.5 by David Jones)", horizfp);
  fputs("\n\n", horizfp);
  for (i = 0; i <= (strlen(confline1) - 1) / 60; i++) {
    fputs("Conf: ", horizfp);
    for (j = i * 60; (j < i * 60 + 60) && (j < origseqlen); j++) {
      putc(confline1[j], horizfp);
    }
    fputs("\n", horizfp);
    fputs("Pred: ", horizfp);
    for (j = i * 60; (j < i * 60 + 60) && (j < origseqlen); j++) {
      putc(predline1[j], horizfp);
    }
    fputs("\n", horizfp);
    fputs("  AA: ", horizfp);
    for (j = i * 60; (j < i * 60 + 60) && (j < origseqlen); j++) {
      putc(aaline1[j], horizfp);
    }
    fputs("\n", horizfp);
    fputs("      ", horizfp);
    for (j = i * 60; (j < i * 60 + 60) && (j < origseqlen); j++) {
      if ((j + 1) % 10 == 0) {
        sprintf(ss2line, "%10d", j + 1);
        fputs(ss2line, horizfp);
      }
    }
    fputs("\n", horizfp);
    fputs("\n", horizfp);
  }
  fclose(horizfp);
}
