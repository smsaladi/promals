#include "ss_prof.h"
#include "regularizer.h"

ss_prof::ss_prof() {
  sslen = 0;
  valid_file = 0;
}

ss_prof::ss_prof(char *base_name) {
  sslen = 0;
  valid_file = 0;
  read_ss_files(base_name);
}

void ss_prof::read_ss_files(char *base_name) {
  int i, j;

  done_prof = 0;

  char horiz_file[200];
  char ss_file[200];
  char line[200], tmpstr[200];
  char tmptype;

  strcpy(horiz_file, base_name);
  strcat(horiz_file, ".horiz");
  strcpy(ss_file, base_name);
  strcat(ss_file, ".ss2");

  // read the ss file
  sslen = 0;
  ifstream ssfp(ss_file, ios::in);
  if (!ssfp) {
    cout << "Secondary structure prediction file does not exist" << endl;
    cout << " --- " << ss_file << endl;
    return;
  }

  ssfp.getline(line, 200);
  ssfp.getline(line, 200);
  while (ssfp) {
    ssfp.getline(line, 200);
    if (ssfp) {
      sslen++;
    }
  }
  ssfp.close();

  sstype = ivector(sslen);
  seq = cvector(sslen + 1);
  ssfreq = gmatrix<float>(sslen, 3);
  ssrel = ivector(sslen);
  alphabet1 = ivector(sslen);

  // important, set ssfreq[0] to be the background frequency vector
  ssfreq[0][1] = rssfreq[1];
  ssfreq[0][2] = rssfreq[2];
  ssfreq[0][3] = rssfreq[3];

  ifstream ssfp1;
  ssfp1.open(ss_file, ios::in);
  ssfp1.getline(line, 200);
  ssfp1.getline(line, 200);
  i = 0;
  float freq_all = 0;
  while (ssfp1.good()) {
    i++;
    freq_all = 0;
    ssfp1 >> tmpstr;
    if (!ssfp1) break;
    // cout << tmpstr << " ";
    ssfp1 >> seq[i - 1];
    ssfp1 >> tmptype;
    // cout << tmptype << " ";
    if (tmptype == 'H')
      sstype[i] = 1;
    else if (tmptype == 'E')
      sstype[i] = 2;
    else
      sstype[i] = 3;
    ssfp1 >> ssfreq[i][3];  // frequencies in the order of "C", "H", "E"
    ssfreq[i][3] += 0.01;
    freq_all += ssfreq[i][3];  // my order is "H", "E", "C"
    ssfp1 >> ssfreq[i][1];
    ssfreq[i][1] += 0.01;
    freq_all += ssfreq[i][1];
    ssfp1 >> ssfreq[i][2];
    ssfreq[i][2] += 0.01;
    freq_all += ssfreq[i][2];

    ssfreq[i][1] /= freq_all;
    ssfreq[i][2] /= freq_all;
    ssfreq[i][3] /= freq_all;

    if (ssfp1.eof()) break;
  }
  seq[sslen] = '\0';
  ssfp1.close();

  // read the horiz file
  i = 0;
  ifstream horizfp(horiz_file, ios::in);
  if (!horizfp) {
    cout << "Warning: secondary structure prediction file does not exist"
         << endl;
    cout << " --- " << horiz_file << endl;
    return;
  }
  while (horizfp.good()) {
    horizfp.getline(line, 200);
    if (strncmp(line, "Conf:", 5) == 0) {
      for (j = 6; j < strlen(line); j++) {
        i++;
        // cout << line[j];
        ssrel[i] = (int)(line[j]) - 48;
        // cout << i << " " << ssrel[i] << endl;
      }
    }
  }
  horizfp.close();

  get_alphabet1();

  valid_file = 1;
  done_prof = 1;
}

ss_prof::~ss_prof() {
  delete[] sstype;
  delete[] seq;
  delete[] ssrel;
  free_gmatrix<float>(ssfreq, sslen, 3);
  delete[] alphabet1;
}

void ss_prof::print_ss_info() {
  int i;

  cout << "secondary structure profile" << endl;
  cout << "num aa ss alphabet1 rel   freq1    freq2    freq3" << endl;

  // cout << sslen << endl;
  for (i = 1; i <= sslen; i++) {
    cout << setw(4) << i << setw(3) << seq[i - 1] << setw(3) << sstype[i]
         << right << setw(5) << alphabet1[i] << "     " << setw(4) << left
         << ssrel[i] << fixed << setprecision(6) << " " << ssfreq[i][1] << " "
         << ssfreq[i][2] << " " << ssfreq[i][3] << endl;
  }
}

void ss_prof::get_alphabet1() {
  int i, j;

  int hcount = 0, ccount = 0, ecount = 0;
  int previous = -1;
  for (i = 1; i <= sslen; i++) {
    if (sstype[i] == 1) {
      hcount++;
      if (sstype[i] != previous) {
        ccount = 0;
        ecount = 0;
        alphabet1[i] = 1;
        if (previous == 2)
          alphabet1[i - 1] = 6;
        else if (previous == 3)
          alphabet1[i - 1] = 9;
      } else {
        alphabet1[i] = 2;
      }
    }
    if (sstype[i] == 2) {
      ecount++;
      if (sstype[i] != previous) {
        hcount = 0;
        ccount = 0;
        alphabet1[i] = 4;
        if (previous == 1)
          alphabet1[i - 1] = 3;
        else if (previous == 3)
          alphabet1[i - 1] = 9;
      } else {
        alphabet1[i] = 5;
      }
    }
    if (sstype[i] == 3) {
      ccount++;
      if (sstype[i] != previous) {
        hcount = 0;
        ecount = 0;
        alphabet1[i] = 7;
        if (previous == 1)
          alphabet1[i - 1] = 3;
        else if (previous == 2)
          alphabet1[i - 1] = 6;
      } else {
        alphabet1[i] = 8;
      }
    }
    previous = sstype[i];
  }
}

int ss_prof::check_aa_seq(char *target_seq) {
  int i, j;

  // cout << "checking sequence" << endl;
  // cout << target_seq << endl;
  // cout << seq << endl;
  //
  char *uppercase_seq = cvector(strlen(target_seq) + 1);
  strcpy(uppercase_seq, target_seq);
  for (i = 0; i < strlen(uppercase_seq); i++)
    uppercase_seq[i] = toupper(uppercase_seq[i]);

  if (strcmp(uppercase_seq, seq) == 0) {
    delete[] uppercase_seq;
    return 1;
  } else {
    delete[] uppercase_seq;
    return 0;
  }
}

ss_prof *ss_prof::sub_ss_prof(int start, int end) {
  int i, j;

  ss_prof *tmpss = new ss_prof();

  // cout << "start: " << start << " end: " << end << endl;
  tmpss->sslen = end - start + 1;
  tmpss->sstype = ivector(tmpss->sslen);
  tmpss->seq = cvector(tmpss->sslen + 1);
  tmpss->ssfreq = gmatrix<float>(tmpss->sslen, 3);
  tmpss->ssrel = ivector(tmpss->sslen);
  tmpss->alphabet1 = ivector(tmpss->sslen);

  for (i = 1; i <= tmpss->sslen; i++) {
    tmpss->sstype[i] = sstype[i + start - 1];
    tmpss->seq[i - 1] = seq[i + start - 1 - 1];
    tmpss->ssrel[i] = ssrel[i + start - 1];
    tmpss->alphabet1[i] = alphabet1[i + start - 1];
    for (j = 1; j <= 3; j++) {
      tmpss->ssfreq[i][j] = ssfreq[i + start - 1][j];
    }
    // cout <<
    // tmpss->sstype[i]<<tmpss->seq[i-1]<<tmpss->ssrel[i]<<tmpss->alphabet1[i]<<"
    // " << tmpss->ssfreq[i][1] << " " << tmpss->ssfreq[i][2] <<
    // tmpss->ssfreq[i][3] << endl;
  }
  tmpss->seq[i - 1] = '\0';
  // important, set ssfreq[0] to be the background frequency vector
  tmpss->ssfreq[0][1] = rssfreq[1];
  tmpss->ssfreq[0][2] = rssfreq[2];
  tmpss->ssfreq[0][3] = rssfreq[3];
  tmpss->valid_file = 1;
  tmpss->done_prof = 1;

  return tmpss;
}

void ss_prof::print_result() {
  int i, j;

  cout << "result of ssprof" << endl;
  cout << "sslen: " << sslen << endl;
  cout << "# seq sstype ssrel alphabet1 freq1 freq2 freq3 " << endl;
  for (i = 1; i <= sslen; i++) {
    cout << i << ": " << seq[i - 1] << " " << sstype[i] << " " << ssrel[i]
         << " ";
    cout << alphabet1[i] << " " << ssfreq[i][1] << " " << ssfreq[i][2] << " "
         << ssfreq[i][3] << endl;
  }
  cout << endl;
}
