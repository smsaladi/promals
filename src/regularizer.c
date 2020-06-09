#include "regularizer.h"
#include "util.h"

// transition probabilities from or to begin and end states
float rsn;
float rsb;
float rnb;
float rnn;
float rec;
float ret;
float rcc;
float rct;

// transition probabilities from or to match states
float rbm[12];
float rme[12];

// transition probabilities among match, insert, delete states
float **rtrans0;
float ***rtrans;
// float rtrans[10][4][4];

// background residue pair emission probabilites
float **rqmatrix0;
float ***rqmatrix;
// float rqmatrix[10][21][21];

// background residue emission probabilites
float *rbfreq0;
float **rbfreq;
// float rbfreq[10][21];
float *log_rbfreq0, **log_rbfreq;

// background predicted secondary structure type frequencies
// from a database of psipred files
/*
473314 1
282112 2
543747 3
0.36431945553055672
0.21714736990377725
0.41853317456566602

real secondary structure type frequencies from a database of dssp files
273188 -
14099 B
275140 E
48047 G
423538 H
229 I
117718 S
147936 T

H = 471814
E = 289239
C = 538842

0.36296316240927151
0.22250951038353098
0.41452732720719748

*/
float rssfreq[4], log_rssfreq[4];

int done_read_regularizer = 0;

void read_regularizer(char *file_name) {
  int i, j, k, l;

  rssfreq[1] = 0.36431945553055672;
  rssfreq[2] = 0.21714736990377725;
  rssfreq[3] = 0.41853317456566602;
  for (i = 1; i <= 3; i++) log_rssfreq[i] = log(rssfreq[i]);

  // allocate memories
  rtrans = new float **[10];
  for (i = 0; i < 10; i++) rtrans[i] = gmatrix<float>(3, 3);
  rqmatrix = new float **[10];
  for (i = 0; i < 10; i++) rqmatrix[i] = gmatrix<float>(20, 20);
  rbfreq = gmatrix<float>(9, 20);
  log_rbfreq = gmatrix<float>(9, 20);

  char line[1000];

  ifstream fp(file_name, ios::in);

  if (!fp) {
    cout << "Cannot read the parameter file: " << file_name << endl;
    exit(0);
  }

  // sn, sb, nb, nn, ec, et, cc, ct
  while (fp.good()) {
    fp >> line;
    if (strncmp(line, "sn:", 3) == 0) {
      fp >> line;
      fp >> rsn;
    }
    if (strncmp(line, "sb:", 3) == 0) {
      fp >> line;
      fp >> rsb;
    }
    if (strncmp(line, "nb:", 3) == 0) {
      fp >> line;
      fp >> rnb;
    }
    if (strncmp(line, "nn:", 3) == 0) {
      fp >> line;
      fp >> rnn;
    }
    if (strncmp(line, "ec:", 3) == 0) {
      fp >> line;
      fp >> rec;
    }
    if (strncmp(line, "et:", 3) == 0) {
      fp >> line;
      fp >> ret;
    }
    if (strncmp(line, "cc:", 3) == 0) {
      fp >> line;
      fp >> rcc;
    }
    if (strncmp(line, "ct:", 3) == 0) {
      fp >> line;
      fp >> rct;
    }
    if (strncmp(line, "bm:", 3) == 0) {
      break;
    }
  }

  // bm and me
  j = 10;
  for (i = 1; i <= j; i++) {
    fp >> line;
    fp >> line;
    fp >> rbm[i];
  }
  if (j == 10) {
    fp >> line;
    fp >> rbm[11];
  }
  while (fp.good()) {
    fp >> line;
    if (strncmp(line, "me:", 3) == 0) {
      break;
    }
  }
  for (i = 1; i <= j; i++) {
    fp >> line;
    fp >> line;
    fp >> rme[i];
  }
  if (j == 10) {
    fp >> line;
    fp >> rme[11];
  }

  // transitions
  j = 0;
  while (fp.good()) {
    fp >> line;
    if (strcmp(line, "env") == 0) {
      fp >> line;
      fp >> line;
      fp >> line;
      fp >> line;
      fp >> rtrans[j][1][1];
      fp >> rtrans[j][1][2];
      fp >> rtrans[j][1][3];
      fp >> line;
      fp >> line;
      fp >> line;
      fp >> rtrans[j][2][1];
      fp >> rtrans[j][2][2];
      fp >> rtrans[j][2][3];
      fp >> line;
      fp >> line;
      fp >> line;
      fp >> rtrans[j][3][1];
      fp >> rtrans[j][3][2];
      fp >> rtrans[j][3][3];
      j++;
      if (j == 10) break;
    }
  }

  rtrans0 = rtrans[0];

  // match probability matrices
  j = 0;
  while (fp.good()) {
    fp >> line;
    if (strcmp(line, "matrix") == 0) {
      fp.getline(line, 1000);
      fp.getline(line, 1000);
      fp.getline(line, 1000);

      for (k = 1; k <= 20; k++) {
        fp >> line;
        for (l = 1; l <= 20; l++) {
          fp >> rqmatrix[j][k][l];
        }
      }
      fp.getline(line, 1000);
      fp.getline(line, 1000);
      fp.getline(line, 1000);
      for (k = 1; k <= 20; k++) {
        fp >> line;
        fp >> rbfreq[j][k];
        log_rbfreq[j][k] = log(rbfreq[j][k]);
      }

      j++;
      if (j >= 10) break;
    }
  }

  rbfreq0 = rbfreq[0];
  log_rbfreq0 = log_rbfreq[0];
  rqmatrix0 = rqmatrix[0];

  fp.close();
}

void print_regularizer() {
  int i, j, k;

  cout << "rsn: " << rsn << endl;
  cout << "rsb: " << rsb << endl;
  cout << "rnb: " << rnb << endl;
  cout << "rnn: " << rnn << endl;
  cout << "rec: " << rec << endl;
  cout << "ret: " << ret << endl;
  cout << "rcc: " << rcc << endl;
  cout << "rct: " << rct << endl;

  cout << "rbm:" << endl;
  for (i = 1; i <= 11; i++) {
    cout << i << "  " << rbm[i] << endl;
  }

  cout << "rme:" << endl;
  for (i = 1; i <= 11; i++) {
    cout << i << "  " << rme[i] << endl;
  }

  cout << "rtrans:" << endl;
  for (i = 0; i <= 9; i++) {
    cout << i << endl;
    for (j = 1; j <= 3; j++) {
      for (k = 1; k <= 3; k++) {
        cout << rtrans[i][j][k] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }

  cout << "rmatrix and rbfreq:" << endl;
  for (i = 0; i <= 9; i++) {
    cout << "matrix: " << i << endl;
    for (j = 1; j <= 20; j++) {
      for (k = 1; k <= 20; k++) {
        cout << rqmatrix[i][j][k] << " ";
      }
      cout << endl;
    }
    cout << endl;

    cout << "rbfreq: " << i << endl;
    for (j = 1; j <= 20; j++) {
      cout << rbfreq[i][j] << endl;
    }
    cout << endl;
  }

  done_read_regularizer = 1;
}
