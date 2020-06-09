#include "constraint.h"

static int debug = 11;

constraint::~constraint() {
  if (nseqs != 0) {
    delete[] rep_index;
    delete[] inside_index;
    delete[] startp;
  }
}

/*
void constraint::read_seqs(char *file_list_name) {

        int i, j, k;
        char fileName[200];

        ifstream fp(file_list_name, ios::in);

        while(fp>>fileName) {
                if(strlen(fileName)==0) { break;}
                seqs sa;
                sa.readFasta(fileName, 0);
                sa.isAlign = 1;
                sv.push_back(sa);
                sa.rep_index = new int [sa.nseqs+1];
        }
}
*/

void constraint::assign_prealigned(vector<tnode *> *prealn) {
  prealigned = prealn;
}

void constraint::assign_originalseq(sequences *os) { oseq = os; }

void constraint::allocate_index() {
  rep_index = new int[nseqs + 1];
  inside_index = new int[nseqs + 1];
  startp = new int[nseqs + 1];
  for (int i = 0; i <= nseqs; i++) {
    rep_index[i] = inside_index[i] = startp[i] = 0;
  }
}

// rep_index:
// >0: correct match, selected as representative sequence
// -1: name does not match
// -2: sequence does not match
// -3: redundant sequences (only one sequence is selected for each prealigned
// group) inside_index: 0: representative sequence itself >0: non-representative
// sequence -1: when rep_index is <0
void constraint::checkNamesSequences() {
  int i, j, k;
  tnode *tn;
  // cout << "nseqs: " << nseqs << endl;

  // find the name index for each sequence; record the two indexes
  for (i = 1; i <= nseqs; i++) {
    // cout << "name: " << name[i] << endl;
    rep_index[i] = -1;
    for (j = 0; j < prealigned->size(); j++) {
      if (rep_index[i] != -1) break;
      tn = (*prealigned)[j];
      if (strcmp(tn->aln->aname[0], name[i].c_str()) == 0) {
        rep_index[i] = j;
        inside_index[i] = 0;
        continue;
      }
      if (!tn->similarSet) continue;
      // cout << "nal: " << tn->similarSet->nal << endl;
      for (k = 1; k <= tn->similarSet->nal; k++) {
        // cout << tn->similarSet->aname[k-1] << endl;
        if (strcmp(tn->similarSet->aname[k - 1], name[i].c_str()) == 0) {
          // cout << "Here" << endl;
          rep_index[i] = j;
          // cout << "Here" << endl;
          inside_index[i] = k;
          // cout << "Here" << endl;
        }
      }
    }
    // cout << "rep_index: " << rep_index[i] << endl;
  }
  // cout << "here" << endl;
  // check if the sequence match the original sequence
  for (i = 1; i <= nseqs; i++) {
    if (rep_index[i] == -1) {
      cout << "The sequence with name " << name[i]
           << " does not match any names in original sequences" << endl;
      continue;
    }
    tn = (*prealigned)[rep_index[i]];
    string a1;
    if (inside_index[i] == 0) {
      a1 = tn->aln->aseq[0];
    } else
      a1 = removegap(tn->similarSet->aseq[inside_index[i] - 1], 1);
    string a2 = removegap(seq[i], 1);
    string::size_type loc = a1.find(a2, 0);
    if (loc != string::npos) {
      startp[i] = loc;
    } else {
      cout << "the sequence with name " << name[i]
           << " does not match the one in the original sequence" << endl;
      cout << a1 << endl;
      cout << a2 << endl;
      startp[i] = -1;
      rep_index[i] = -2;
    }
  }
  // remove redundancy: representatives (inside_index: 0) are in priority,
  // otherwise, remove redundancy: the first occurance is kept, the others
  // removed for sequences within the same prealigned group
  int mark[prealigned->size()];
  for (i = 0; i < prealigned->size(); i++) mark[i] = 0;
  // mark representatives
  for (i = 1; i <= nseqs; i++) {
    if (rep_index[i] < 0) continue;
    if (inside_index[i] == 0) mark[rep_index[i]] = 1;
  }
  // keep one non-representative sequence (the first one) if representative
  // sequence is not available
  for (i = 1; i <= nseqs; i++) {
    if (rep_index[i] < 0) continue;
    if (inside_index[i] == 0) continue;
    if (mark[rep_index[i]] != 0) {
      rep_index[i] = -3;
    } else {
      mark[rep_index[i]] = 1;
    }
  }
  if (debug > 1) {
    printIndexSeq();
  }
}

// print sequences with information
void constraint::printIndexSeq() {
  for (int i = 1; i <= nseqs; i++) {
    fprintf(stdout, "%20s %3d %3d %3d %s\n", name[i].c_str(), rep_index[i],
            inside_index[i], startp[i], seq[i].c_str());
  }
}

// remove gap letters in a string
string constraint::removegap(string a1, int touppercase) {
  int i, j;
  string a2 = "";
  for (i = 0; i < a1.size(); i++) {
    if ((a1[i] != '-') && (a1[i] != '.') && (a1[i] != ' ')) {
      if (touppercase)
        a2 += toupper(a1[i]);
      else
        a2 += a1[i];
    }
  }
  return a2;
}

// for sequences with index ci and index cj, find alignment matrix for
// their representative sequences
float **constraint::get_constraint_matrix(int ci, int cj) {
  int i, j, k;
  int *correspond1 = NULL, *correspond2 = NULL;
  int replen1, replen2;
  char *tmpseq1, *tmpseq2;
  int len1, len2;  // length of target sequence
  int ind1, ind2;  // index of representative sequence in similarSet
  int m1, m2;

  // the dimensions of the score matrix
  replen1 = (*prealigned)[rep_index[ci]]->aln->alilen;
  replen2 = (*prealigned)[rep_index[cj]]->aln->alilen;
  // cout << "replen1: " << replen1 << endl;
  // cout << "replen2: " << replen2 << endl;
  float **score_matrix = gmatrix<float>(replen1, replen2);
  for (i = 1; i <= replen1; i++)
    for (j = 1; j <= replen2; j++) score_matrix[i][j] = 0;

  // find the representative correspondance
  // case 1: the sequence itself is representative sequence
  // in this case, correspondance is simply the array: 1, 2, 3, ...
  if (inside_index[ci] == 0) {
    correspond1 = new int[replen1 + 1];
    for (i = 1; i <= replen1; i++) correspond1[i] = i;
    // cout << "this place 1" << endl;
  }
  if (inside_index[cj] == 0) {
    correspond2 = new int[replen2 + 1];
    for (i = 1; i <= replen2; i++) correspond2[i] = i;
    // cout << "this place 2" << endl;
  }
  // case 2: the sequence is not a representative sequence
  // in this case, the correspondance is derived from the pairwise alignment of
  // target sequence and representative sequence in the similarSet
  if (inside_index[ci] != 0) {
    // find the index of the representative in the similarSet
    for (i = 0; i < (*prealigned)[rep_index[ci]]->similarSet->nal; i++) {
      if (strcmp((*prealigned)[rep_index[ci]]->aln->aname[0],
                 (*prealigned)[rep_index[ci]]->similarSet->aname[i]) == 0)
        ind1 = i;
    }
    tmpseq1 =
        (*prealigned)[rep_index[ci]]->similarSet->aseq[inside_index[ci] - 1];
    tmpseq2 = (*prealigned)[rep_index[ci]]->similarSet->aseq[ind1];
    correspond1 = get_correspond(tmpseq1, tmpseq2, len1);
  }
  if (inside_index[cj] != 0) {
    // find the index of the representative in the similarSet
    for (i = 0; i < (*prealigned)[rep_index[cj]]->similarSet->nal; i++) {
      if (strcmp((*prealigned)[rep_index[cj]]->aln->aname[0],
                 (*prealigned)[rep_index[cj]]->similarSet->aname[i]) == 0)
        ind2 = i;
    }
    tmpseq1 =
        (*prealigned)[rep_index[cj]]->similarSet->aseq[inside_index[cj] - 1];
    tmpseq2 = (*prealigned)[rep_index[cj]]->similarSet->aseq[ind2];
    correspond2 = get_correspond(tmpseq1, tmpseq2, len1);
  }

  // now derive the matrix
  int min_alilen = seq[ci].size();
  if (min_alilen > seq[cj].size()) min_alilen = seq[cj].size();
  // cout << "min_length: " << min_alilen << endl;
  m1 = startp[ci];
  m2 = startp[cj];
  for (i = 0; i < min_alilen; i++) {
    if (seq[ci][i] != '-') {
      m1++;
    }
    if (seq[cj][i] != '-') {
      m2++;
    }
    if ((seq[ci][i] != '-') && (seq[cj][i] != '-')) {
      // cout << "MMMMM m1: " << seq[ci][i] << " m2: " << seq[cj][i] << endl;
      // ignore the cases for small letters
      if (seq[ci][i] != toupper(seq[ci][i])) {
        ;
      } else if (seq[cj][i] != toupper(seq[cj][i])) {
        ;
      } else if ((correspond1[m1] != -1) && (correspond2[m2] != -1)) {
        score_matrix[correspond1[m1]][correspond2[m2]] = 1.0;
        cout << "m1: " << correspond1[m1] << " m2: " << correspond2[m2] << endl;
      }
    }
  }

  // cout << "here" << endl;
  if (correspond1) delete[] correspond1;
  if (correspond2) delete[] correspond2;
  return score_matrix;
}

// for two aligned sequences, find residue correspondances
int *get_correspond(char *tmpseq1, char *tmpseq2, int &len) {
  int i, j;
  int *correspond;
  int m1, m2;

  len = 0;
  for (i = 0; i < strlen(tmpseq1); i++) {
    if (tmpseq1[i] != '-') len++;
  }
  correspond = new int[len + 1];
  // get the index
  m1 = 0;
  m2 = 0;
  for (i = 0; i < strlen(tmpseq1); i++) {
    if (tmpseq1[i] != '-') {
      m1++;
    }
    if (tmpseq2[i] != '-') {
      m2++;
    }
    if (tmpseq1[i] != '-') {
      if (tmpseq2[i] != '-') {
        correspond[m1] = m2;
      } else {
        correspond[m1] = -1;
      }
    }
  }
  // cout << tmpseq1 << endl;
  // cout << tmpseq2 << endl;
  for (i = 1; i <= len; i++) {
    ;
    // cout << "correspond: " << i << " " << correspond[i] << endl;
  }
  return correspond;
}

vector<constraint *> read_multiple_constraint(char *filename) {
  int i, j, k;

  ifstream fp(filename, ios::in);

  vector<constraint *> vcons;
  int count = 0;
  while (fp.good()) {
    constraint *cons = new constraint();
    cons->readFasta(fp, 0);
    vcons.push_back(cons);
    count++;
  }

  fp.close();

  return vcons;
}
