#include "sequences.h"

static int n6 = 46656;
static int n5 = 7776;
static int n4 = 1296;
static int n3 = 216;
static int n2 = 36;
static int n1 = 6;

#define MIN_HERE(x, y) (((x) < (y)) ? (x) : (y))

sequences::sequences() {
  nseqs = 0;
  distMat = 0;
}

void sequences::get_map() {
  int i;

  for (i = 1; i <= nseqs; i++) {
    name2seq.insert(pair<string, string>(name[i], seq[i]));
    name2index.insert(pair<string, int>(name[i], i));
  }

  // for(i=1;i<=nseqs;i++) { cout << i << " " << name2index.find(name[i])->first
  // << " " << name2index.find(name[i])->second << endl; }
}

// copy constructor
sequences::sequences(const sequences &s) {
  int i;
  string str;

  nseqs = s.nseqs;
  isAlign = s.isAlign;
  name.erase(name.begin(), name.end());
  seq.erase(seq.begin(), seq.end());
  // unlike subalign, index begins from 1 and not 0
  name.push_back("");
  seq.push_back("");
  for (i = 1; i <= nseqs; i++) {
    str = s.name[i];
    name.push_back(str);
    // name.push_back(s.name[i]);
    str = s.seq[i];
    seq.push_back(str);
    // seq.push_back(s.seq[i]);
  }

  distMat = 0;
}

// copy constructor, input is a subalign
sequences::sequences(const subalign &s) {
  int i;

  nseqs = s.nal;
  isAlign = 1;
  name.erase(name.begin(), name.end());
  seq.erase(seq.begin(), seq.end());
  // unlike subalign, index begins from 1 and not 0
  name.push_back("");
  seq.push_back("");
  for (i = 0; i < s.nal; i++) {
    name.push_back(s.aname[i]);
    seq.push_back(s.aseq[i]);
  }

  distMat = 0;
}

// read input from a fasta format file; remove non-alphabetic letters
// read input from a fasta format file; remove non-alphabetic letters
void sequences::readFasta(char *inputFastaName, int zap, int disallow_one_seq) {
  int i, j, k;
  string tmpstring;

  nseqs = 0;
  isAlign = 0;
  name.push_back("");
  seq.push_back("");

  ifstream fastaFile(inputFastaName);
  if (!fastaFile) {
    cout << "Failed to read the fasta file: " << inputFastaName << endl;
    exit(1);
  }

  string disallowed_letters = string(" \t?'\"`&|\\{}*()/[]$;");
  while (fastaFile) {
    fastaFile >> tmpstring;
    // cout << tmpstring << "===" << endl;
    if (tmpstring[0] == '>') {
      char tmp1[10000];
      fastaFile.getline(tmp1, 10000);
      if (tmp1[strlen(tmp1) - 1] == 13) tmp1[strlen(tmp1) - 1] = '\0';
      if ((strlen(tmp1) == 1) && (tmp1[0] == 13)) {
        ;
      }
      // cout << strlen(tmp1) << "|" << tmp1 << "|" << endl;
      else
        tmpstring = tmpstring + tmp1;
      tmpstring = tmpstring.substr(0, 100);
      // cout << "|||||||||||"<< tmp1 << "------"<< endl;
      nseqs++;
      // cout << nseqs << endl;
      string tmpstring1 = tmpstring.substr(1);
      tmpstring = "";
      for (i = 0; i < tmpstring1.length(); i++) {
        string::size_type pos =
            disallowed_letters.find(tmpstring1.substr(i, 1), 0);
        if (pos != string::npos) {
          tmpstring += '_';
        } else {
          tmpstring += tmpstring1[i];
        }
      }
      // name.push_back(tmpstring.substr(0, 15) );
      name.push_back(tmpstring.substr(0, 25));
      seq.push_back("");
      // cout << tmpstring << "===" << endl;
    } else {
      seq[nseqs] += tmpstring;
    }
    tmpstring = "";
  }

  if (disallow_one_seq)
    if (nseqs == 1) {
      cout << "Error: your input data contains only one sequence." << endl;
      cout << "       Please input at least two sequences for alignment "
              "construction."
           << endl;
      exit(0);
    }

  for (i = 1; i <= nseqs; i++) {
    if (zap) zapLetters(seq[i]);
    toUpper(seq[i]);
    // cout << name[i] << endl;
    // cout << seq[i] << endl;
  }

  // check for identical file names
  int same_name_check = 0;
  cout << endl;
  for (i = 1; i <= nseqs; i++) {
    if (name[i].size() == 0) {
      cout << endl;
      cout << "Error: sequence number " << i << " does not have a name."
           << endl;
      cout << "Please provide a name for every sequence." << endl;
      cout << endl;
      exit(0);
    }
    for (j = i + 1; j <= nseqs; j++) {
      if (name[i] == name[j]) {
        cout << "sequence number " << i << " and sequence number " << j
             << " have the same name: " << name[i] << endl;
        same_name_check = 1;
      }
    }
  }
  if (same_name_check) {
    cout << endl;
    cout << "Error message:" << endl;
    cout << "       - there should be no two sequences with the same name in "
            "your dataset"
         << endl;
    cout << "       - please change your sequence names and try again" << endl;
    exit(0);
  }

  // check for zero length sequences
  int zero_len_check = 0;
  for (i = 1; i <= nseqs; i++) {
    if (seq[i].length() == 0) {
      cout << "Error message: " << endl;
      cout << "      - The sequence with name '" << name[i]
           << "' has a zero length" << endl;
      zero_len_check = 1;
    }
  }
  if (zero_len_check) {
    exit(0);
  }
}

// read input from a fasta format file; remove non-alphabetic letters
void sequences::readFasta(ifstream &fastaFile, int zap) {
  int i, j, k;
  string tmpstring;

  nseqs = 0;
  isAlign = 0;
  name.push_back("");
  seq.push_back("");

  // ifstream fastaFile(inputFastaName);
  /*
  if(!fastaFile) {
          cout << "Failed to read the fasta file: " << inputFastaName << endl;
          exit(1);
  }*/

  string disallowed_letters = string(" \t?'\"`&|\\{}*()/[]$;");
  while (fastaFile) {
    fastaFile >> tmpstring;
    // cout << tmpstring << "===" << endl;
    // here add the condition to stop
    if (tmpstring[0] == '@') break;
    if (tmpstring[0] == '>') {
      char tmp1[10000];
      fastaFile.getline(tmp1, 10000);
      if (tmp1[strlen(tmp1) - 1] == 13) tmp1[strlen(tmp1) - 1] = '\0';
      if ((strlen(tmp1) == 1) && (tmp1[0] == 13)) {
        ;
      }
      // cout << strlen(tmp1) << "|" << tmp1 << "|" << endl;
      else
        tmpstring = tmpstring + tmp1;
      tmpstring = tmpstring.substr(0, 100);
      // cout << "|||||||||||"<< tmp1 << "------"<< endl;
      nseqs++;
      // cout << nseqs << endl;
      string tmpstring1 = tmpstring.substr(1);
      tmpstring = "";
      for (i = 0; i < tmpstring1.length(); i++) {
        string::size_type pos =
            disallowed_letters.find(tmpstring1.substr(i, 1), 0);
        if (pos != string::npos) {
          tmpstring += '_';
        } else {
          tmpstring += tmpstring1[i];
        }
      }
      // name.push_back(tmpstring.substr(0, 15) );
      name.push_back(tmpstring.substr(0, 25));
      seq.push_back("");
      // cout << tmpstring << "===" << endl;
    } else {
      seq[nseqs] += tmpstring;
    }
    tmpstring = "";
  }

  /*
  if(nseqs==1) {
          cout << "Error: your input data contains only one sequence." << endl;
          cout << "       Please input at least two sequences for alignment
  construction." << endl; exit(0);
  }
  */

  for (i = 1; i <= nseqs; i++) {
    if (zap) zapLetters(seq[i]);
    toUpper(seq[i]);
    // cout << name[i] << endl;
    // cout << seq[i] << endl;
  }

  // check for identical file names
  int same_name_check = 0;
  cout << endl;
  for (i = 1; i <= nseqs; i++) {
    if (name[i].size() == 0) {
      cout << endl;
      cout << "Error: sequence number " << i << " does not have a name."
           << endl;
      cout << "Please provide a name for every sequence." << endl;
      cout << endl;
      exit(0);
    }
    for (j = i + 1; j <= nseqs; j++) {
      if (name[i] == name[j]) {
        cout << "sequence number " << i << " and sequence number " << j
             << " have the same name: " << name[i] << endl;
        same_name_check = 1;
      }
    }
  }
  if (same_name_check) {
    cout << endl;
    cout << "Error message:" << endl;
    cout << "       - there should be no two sequences with the same name in "
            "your dataset"
         << endl;
    cout << "       - please change your sequence names and try again" << endl;
    exit(0);
  }

  // check for zero length sequences
  int zero_len_check = 0;
  for (i = 1; i <= nseqs; i++) {
    if (seq[i].length() == 0) {
      cout << "Error message: " << endl;
      cout << "      - The sequence with name '" << name[i]
           << "' has a zero length" << endl;
      zero_len_check = 1;
    }
  }
  if (zero_len_check) {
    exit(0);
  }
}

// zapLetters: remove letters that are not alphabets such as gaps
inline void sequences::zapLetters(string &str) {
  int i;
  int num = 0;
  string tmpstr;
  for (i = 0; i < str.length(); i++) {
    if ((str[i] >= 'A') && (str[i] <= 'Z')) {
      tmpstr += str[i];
      continue;
    }
    if ((str[i] >= 'a') && (str[i] <= 'z')) {
      tmpstr += str[i];
    }
  }

  str = tmpstr;
  // cout << "|" << str << "|" << endl;
  char c1;

  for (i = 0; i < str.length(); i++) {
    c1 = str[i];
    switch (c1) {
      // blastpgp does not work with [oOjJ], doesnot include them in checkpoint
      // file psipred thus misses these letters and creates problems
      case 'o':
        str[i] = 'A';
        continue;
      case 'O':
        str[i] = 'A';
        continue;
      case 'j':
        str[i] = 'A';
        continue;
      case 'J':
        str[i] = 'A';
        continue;
        // case 'u':str[i] = 'A'; continue;
        // case 'U':str[i] = 'A'; continue;
        // case 'x':str[i] = 'A'; continue;
        // case 'X':str[i] = 'A'; continue;
        // case 'b':str[i] = 'A'; continue;
        // case 'B':str[i] = 'A'; continue;
        // case 'z':str[i] = 'A'; continue;
        // case 'Z': str[i] = 'A'; continue;
    }
  }

  // cout << "|" << str << "|" << endl;

  // cout << tmpstr << endl;
}

// toUpper: change the letters to upper case
inline void sequences::toUpper(string &str) {
  int i;
  for (i = 0; i < str.length(); i++) {
    str[i] = toupper(str[i]);
  }
}

// print out the sequences
void sequences::printSeqs() {
  int i;
  cout << "number of sequences: " << nseqs << endl;
  for (i = 1; i <= nseqs; i++) {
    // zapLetters(seq[i]);
    cout << ">" << name[i] << endl;
    cout << seq[i] << endl;
  }
}

void sequences::output_fasta(char *fasta_file_name) {
  int i, j, k;

  ofstream ofp(fasta_file_name, ios::out);
  if (!ofp) {
    cout << "cannot open " << fasta_file_name << " for writing fasta records"
         << endl;
  }
  for (i = 1; i <= nseqs; i++) {
    ofp << ">" << name[i] << endl << seq[i] << endl;
  }
  ofp.close();
}

// convert to Dayhoff 6-alphabet
void sequences::toDayhoff6() {
  int i, j;
  for (i = 1; i <= nseqs; i++) {
    for (j = 0; j < seq[i].length(); j++) {
      seq[i][j] = am2d6(seq[i][j]);
    }
  }
}

// derive a K-mer map(table, hash) for a sequence(str)
dayhoff6Table sequences::seqToD6t(string str, int K) {
  int i, j;

  dayhoff6Table tmpd6t;

  for (i = 0; i < str.length() - K + 1; i++) {
    tmpd6t[str.substr(i, K)]++;
  }

  return tmpd6t;
}

// generate the vector of K-mer maps
void sequences::generateD6t(int K) {
  int i;
  dayhoff6Table d;

  d6t.erase(d6t.begin(), d6t.end());
  d6t.push_back(d);
  for (i = 1; i <= nseqs; i++) {
    d6t.push_back(seqToD6t(seq[i], K));
  }
}

// the difference between two dayhoff6tables
int sequences::diffCountD2t(dayhoff6Table t1, dayhoff6Table t2) {
  int i, j;
  int count = 0;

  dayhoff6Table::iterator iter1 = t1.begin();

  while (iter1 != t1.end()) {
    // cout << iter1->first << "\t" << iter1->second << endl;
    count += abs(iter1->second - t2[iter1->first]);
    iter1++;
  }

  dayhoff6Table::iterator iter2 = t2.begin();
  while (iter2 != t2.end()) {
    if (t1[iter2->first] == 0) {
      count += iter2->second;
    }
    iter2++;
  }
  // cout << "Count: " << count << endl;
  return count;
}

void sequences::get_kmer_array(int K) {
  int i, j, k;
  int tmpindex;

  kmer_array = gmatrix<short int>(nseqs, n6);
  for (i = 1; i <= nseqs; i++) {
    for (j = 1; j <= n6; j++) {
      kmer_array[i][j] = 0;
    }
  }

  for (i = 1; i <= nseqs; i++) {
    if (seq[i].length() < 6) continue;
    int tmpintseq[seq[i].length()];
    for (j = 0; j < seq[i].length(); j++) {
      tmpintseq[j] = am2d6(seq[i][j]) - 64 - 1;
      // seq[i][j] = am2d6(seq[i][j]);
    }

    for (j = 1; j <= seq[i].length() - K + 1; j++) {
      tmpindex = 0;
      tmpindex = n5 * tmpintseq[j - 1] + n4 * tmpintseq[j] +
                 n3 * tmpintseq[j + 1] + n2 * tmpintseq[j + 2] +
                 n1 * tmpintseq[j + 3] + tmpintseq[j + 4];
      kmer_array[i][tmpindex + 1]++;
    }
  }
}

void sequences::get_kmer_distance(int K) {
  int i, j, k;
  short int commonCount, minLength;
  float log10 = 0.434294482;
  float tmpnum;

  if (!distMat) {
    distMat = gmatrix<float>(nseqs, nseqs);
  }
  // vector<int> a1; // vector implementation of non-zero indexes, slower
  int *a2;  // array implementation of non-zero indexes, faster
  int count;
  for (i = 1; i <= nseqs; i++) {
    // cout << "i: " << i << " " << seq[i].length() << endl;
    // array implementation of non-zero indexes, faster
    a2 = new int[seq[i].length()];
    count = 0;
    for (k = 1; k <= n6; k++) {
      if (kmer_array[i][k]) {
        a2[count] = k;
        count++;
      }
    }

    /*
    //vector implementation of non-zero indexes, slower
    a1.clear();
    for(k=1;k<=n6;k++) if(kmer_array[i][k]) a1.push_back(k);
    if(i==1) {
            for(k=0;k<a1.size();k++) {
                    cout << a1[k] << " " << kmer_array[i][a1[k]] << endl;
            }
    }
    */

    // cout << "HERE" << endl;
    for (j = i + 1; j <= nseqs; j++) {
      commonCount = 0;
      // if(seq[i].length() > seq[j].length() ) { minLength = seq[j].length(); }
      // else minLength = seq[i].length();
      minLength = MIN_HERE(seq[i].length(), seq[j].length());
      if (minLength < 6) minLength = 6;
      // cout << j << " minLength: " << minLength << endl;
      // cout << "len1: " << seq[i].length() << " len2: " << seq[j].length() <<
      // endl; cout << "count: " << count << endl;

      // array implementation of non-zero indexes, faster
      for (k = 0; k < count; k++) {
        // cout << "k: " << k << endl;
        // cout << kmer_array[i][a2[k]] << endl;
        // cout << kmer_array[j][a2[k]] << endl;
        commonCount += MIN_HERE(kmer_array[i][a2[k]], kmer_array[j][a2[k]]);
      }

      /*
      //vector implementation of non-zero indexes, slower
      for(k=0;k<a1.size();k++) {
              commonCount += MIN_HERE(kmer_array[i][a1[k]],
      kmer_array[j][a1[k]]);
      }
      */

      /*
      // original implementation, bad
      for(k=1;k<=n6;k++) {
              if(kmer_array[i][k] || kmer_array[j][k])
              commonCount += MIN_HERE(kmer_array[i][k], kmer_array[j][k]);
      }*/

      // tmpnum = 1.0/(minLength-K+1);
      // distMat[i][j] = log(0.1 + tmpnum * commonCount) * log10;
      // cout << "HERERER" << endl;
      // cout << log(0.1 + 1.0 * commonCount/(minLength-K+1) ) * log10 << endl;
      distMat[i][j] =
          log(0.1 + 1.0 * commonCount / (minLength - K + 1)) * log10;
      // cout << "TEERERER" << endl;
      if (distMat[i][j] >= 0)
        distMat[i][j] = 0;
      else
        distMat[i][j] = 0 - distMat[i][j];
      ////cout << "TEERERER" << endl;
      distMat[j][i] = distMat[i][j];
      // cout << "TEERERER" << endl;
      // cout << i << endl;
      // cout << j << endl;
      ////cout << distMat[i][j] << endl;
      // cout << "i: " << i << " j: " << j << "  " << distMat[i][j] << endl;
    }
    distMat[i][i] = 0;
    delete[] a2;  // array implementation of non-zero indexes, faster
  }

  // printDistMat();

  free_gmatrix<short int>(kmer_array, nseqs, n6);
}

// the common counts between two dayhoff6tables
int sequences::commonCountD2t(dayhoff6Table t1, dayhoff6Table t2) {
  int i, j;
  int count = 0;

  dayhoff6Table::iterator iter1 = t1.begin();

  while (iter1 != t1.end()) {
    // cout << iter1->first << "\t" << iter1->second << endl;
    // exit(0);
    if (iter1->second > t2[iter1->first]) {
      count += t2[iter1->first];
    } else
      count += iter1->second;
    iter1++;
  }

  // dayhoff6Table::iterator iter2 = t2.begin();
  // while(iter2 != t2.end() ) {
  //	if(t1[iter2->first] == 0) { count += iter2->second; }
  //	iter2++;
  //}
  // cout << "Count: " << count << endl;
  return count;
}

void sequences::d6t2DistMat(int K) {
  int i, j;
  int minLength;

  if (!distMat) {
    distMat = gmatrix<float>(nseqs, nseqs);
  }

  for (i = 1; i <= nseqs; i++) {
    for (j = 1; j <= nseqs; j++) {
      if (seq[i].length() > seq[j].length()) {
        minLength = seq[j].length();
      } else
        minLength = seq[i].length();
      // cout << i << "\t" << j << "\t" << minLength << endl;
      // distMat[i][j] = 1; // - 1.0 *
      // commonCountD2t(d6t[i],d6t[j])/(minLength-K+1);
      // distMat[i][j] = 1 - 1.0 *
      // commonCountD2t(d6t[i],d6t[j])/(minLength-K+1);
      distMat[i][j] = log(0.1 + 1.0 * commonCountD2t(d6t[i], d6t[j]) /
                                    (minLength - K + 1)) /
                      log(10.0);
      if (distMat[i][j] >= 0)
        distMat[i][j] = 0;
      else
        distMat[i][j] = 0 - distMat[i][j];
    }
  }
}

void sequences::printDistMat() {
  int i, j;

  cout << "  " << nseqs << endl;
  for (i = 1; i <= nseqs; i++) {
    cout << left << setw(10) << name[i].substr(0, 9) << " ";
    for (j = 1; j <= nseqs; j++) {
      cout << left << setw(10) << setprecision(6) << distMat[i][j];
    }
    cout << endl;
  }
}

void sequences::seqIdentity2DistMat() {  // dist = 1 - (percentage seq. ident.)

  int i, j;

  if (!isAlign) {
    cout << "warning:sequences may not be aligned" << endl;
  }
  if (!distMat) {
    distMat = gmatrix<float>(nseqs, nseqs);
  }
  for (i = 1; i <= nseqs; i++) {
    // cout << seq[i].length() << endl;
    distMat[i][i] = 0;
    for (j = i + 1; j <= nseqs; j++) {
      distMat[i][j] = distMat[j][i] = 1.0;
      distMat[i][j] = distMat[j][i] =
          1 - seqIdentity(seq[i], seq[j], seq[i].length());
    }
  }
}

double seqIdentity(const string &c1, const string &c2, int length) {
  int i, j;
  int ip = 0;
  int len1 = 0, len2 = 0;
  for (i = 0; i < length; i++) {
    if ((c1[i] == c2[i]) && isalpha(c1[i]) && isalpha(c2[i])) {
      ip++;
    }
    if (isalpha(c1[i])) len1++;
    if (isalpha(c2[i])) len2++;
  }
  return (1.0 * ip / ((len1 < len2) ? len1 : len2));
}
