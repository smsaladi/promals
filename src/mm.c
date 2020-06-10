#include "all.h"
#include "mm.h"

// static int debug = 1;

MM::MM() {
  DD = NULL;
  HH = NULL;
  RR = NULL;
  SS = NULL;
  gS = NULL;
  displ = NULL;

  debug = 0;
  last_print = 0;
  print_ptr = 1;
  endgappenalties = 1;
}

MM::MM(int m, int n, int _g, int _gh) {
  int i;

  M = m;
  N = n;
  g = _g;
  gh = _gh;

  if (g < 0) warning("gap opening penalty should be a negative integhr");
  if (gh < 0) warning("gap extension penalty should be a negative integhr");

  debug = 0;

  DD = new int[N + 1];
  HH = new int[N + 1];
  RR = new int[N + 1];
  SS = new int[N + 1];
  gS = new int[N + 1];
  displ = new int[M + N + 2];
  last_print = 0;
  print_ptr = 1;
  endgappenalties = 1;
  cout << " Using this constructer" << endl;
}

MM::~MM() {
  int i;

  if (DD) delete[] DD;
  if (HH) delete[] HH;
  if (SS) delete[] SS;
  if (gS) delete[] gS;
  if (RR) delete[] RR;
  if (displ) delete[] displ;
}

void MM::setM(int x) { M = x; }

void MM::setN(int x) { N = x; }

void MM::set_g(int x) {
  g = x;
  if (g < 0) warning("gap opening penalty should be a negative integhr");
}

void MM::set_h(int x) {
  gh = x;
  if (gh < 0) warning("gap extension penalty should be a negative integhr");
}

void MM::dp() {
  int i;

  if (endgappenalties)
    maxscore = diff(0, 0, M, N, g, g);
  else
    maxscore = diff1(0, 0, M, N, g, g);

  if (debug > 1) {
    cout << "displ: " << endl;
    for (i = 1; i <= print_ptr - 1; i++) {
      cout << i << " " << displ[i] << endl;
    }
  }
}

void MM::dp(int **smat) {
  int i, j;
  scoreMat = smat;
  // for(i=1;i<=M;i++) { for(j=1;j<=N;j++) { cout << scoreMat[i][j] << " "; }
  // cout << endl; } g = 10; gh =10; cout << "g: " << g << "\th: " << gh << endl;

  if (!DD) DD = new int[N + 1];
  if (!HH) HH = new int[N + 1];
  if (!RR) RR = new int[N + 1];
  if (!SS) SS = new int[N + 1];
  if (!gS) gS = new int[N + 1];
  if (!displ) displ = new int[M + N + 2];
  last_print = 0;
  print_ptr = 1;
  endgappenalties = 1;
  maxscore = diff1(0, 0, M, N, g, g);
}

void MM::del(int k) {
  if (last_print < 0)
    last_print = displ[print_ptr - 1] -= k;
  else
    last_print = displ[print_ptr++] = -(k);
}

void MM::add(int v) {
  if (last_print < 0) {
    displ[print_ptr - 1] = v;
    displ[print_ptr++] = last_print;
  } else
    last_print = displ[print_ptr++] = v;
}

void MM::palign() { displ[print_ptr++] = last_print = 0; }

void MM::warning(const char *s) { fprintf(stderr, "%s\n", s); }

int MM::diff1(int A, int B, int m, int n, int go1, int go2) {
  int midi, midj, type;
  int midh;

  static int t, tl, go, h;

  {
    static int i, j;
    static int hh, f, e, s;

    /* Boundary cases: m <= 1 or n == 0 */
    if (debug > 2)
      fprintf(stdout, "A %d B %d m %d n %d midi %d go1 %d go2 %d\n", (int)A,
              (int)B, (int)m, (int)n, (int)m / 2, (int)go1, (int)go2);

    /* if sequence B is empty....                                            */

    if (n <= 0) {
      /* if sequence A is not empty.... */

      if (m > 0) {
        /* delete residues A[1] to A[m] */

        del(m);
      }
      return (-gap_penalty1(A, B, m));
    }

    /* if sequence A is empty....                                            */

    if (m <= 1) {
      if (m <= 0) {
        /* insert residues B[1] to B[n] */

        add(n);
        return (-gap_penalty2(A, B, n));
      }

      /* if sequence A has just one residue.... */

      if (go1 == 0)
        midh = -gap_penalty1(A + 1, B + 1, n);
      else
        midh = -gap_penalty2(A + 1, B, 1) - gap_penalty1(A + 1, B + 1, n);
      midj = 0;
      for (j = 1; j <= n; j++) {
        hh = -gap_penalty1(A, B + 1, j - 1) +
             scoreMat[A + 1][B + j];  // w(A+1, B+j)
        -gap_penalty1(A + 1, B + j + 1, n - j);
        if (hh > midh) {
          midh = hh;
          midj = j;
        }
      }

      if (midj == 0) {
        add(n);
        del(1);
      } else {
        if (midj > 1) add(midj - 1);
        palign();
        if (midj < n) add(n - midj);
      }
      return midh;
    }

    /* Divide sequence A in half: midi */

    midi = m / 2;
    // cout << "Here" << endl;

    /* In a forward phase, calculate all HH[j] and HH[j] */

    HH[0] = 0;
    t = -open_penalty1(A, B + 1);
    tl = -ext_penalty1(A, B + 1);
    // cout << t <<"\t" <<  tl << endl;
    for (j = 1; j <= n; j++) {
      HH[j] = t = t + tl;
      DD[j] = t - open_penalty2(A + 1, B + j);
    }

    if (go1 == 0)
      t = 0;
    else
      t = -open_penalty2(A + 1, B);
    tl = -ext_penalty2(A + 1, B);
    for (i = 1; i <= midi; i++) {
      s = HH[0];
      HH[0] = hh = t = t + tl;
      f = t - open_penalty1(A + i, B + 1);

      for (j = 1; j <= n; j++) {
        go = open_penalty1(A + i, B + j);
        h = ext_penalty1(A + i, B + j);
        if ((hh = hh - go - h) > (f = f - h)) f = hh;
        go = open_penalty2(A + i, B + j);
        h = ext_penalty2(A + i, B + j);
        if ((hh = HH[j] - go - h) > (e = DD[j] - h)) e = hh;
        hh = s + scoreMat[A + i][B + j];  // w(A+i, B+j);// w(B+j, A+i);
        if (f > hh) hh = f;
        if (e > hh) hh = e;

        s = HH[j];
        HH[j] = hh;
        DD[j] = e;
      }
      // cout << endl;
    }

    DD[0] = HH[0];

    /* In a reverse phase, calculate all RR[j] and SS[j] */

    RR[n] = 0;
    tl = 0;
    for (j = n - 1; j >= 0; j--) {
      go = -open_penalty1(A + m, B + j + 1);
      tl -= ext_penalty1(A + m, B + j + 1);
      RR[j] = go + tl;
      SS[j] = RR[j] - open_penalty2(A + m, B + j);
      gS[j] = open_penalty2(A + m, B + j);
    }

    tl = 0;
    for (i = m - 1; i >= midi; i--) {
      s = RR[n];
      if (go2 == 0)
        go = 0;
      else
        go = -open_penalty2(A + i + 1, B + n);
      tl -= ext_penalty2(A + i + 1, B + n);
      RR[n] = hh = go + tl;
      t = open_penalty1(A + i, B + n);
      f = RR[n] - t;

      for (j = n - 1; j >= 0; j--) {
        go = open_penalty1(A + i, B + j + 1);
        h = ext_penalty1(A + i, B + j + 1);
        if ((hh = hh - go - h) > (f = f - h - go + t)) f = hh;
        t = go;
        go = open_penalty2(A + i + 1, B + j);
        h = ext_penalty2(A + i + 1, B + j);
        hh = RR[j] - go - h;
        if (i == (m - 1)) {
          e = SS[j] - h;
        } else {
          e = SS[j] - h - go + open_penalty2(A + i + 2, B + j);
          gS[j] = go;
        }
        if (hh > e) e = hh;
        hh = s + scoreMat[A + i + 1]
                         [B + j + 1];  // w(A+i+1,B+j+1); //w(B+j+1, A+i+1);
        if (f > hh) hh = f;
        if (e > hh) hh = e;

        s = RR[j];
        RR[j] = hh;
        SS[j] = e;
      }
    }
    SS[n] = RR[n];
    gS[n] = open_penalty2(A + midi + 1, B + n);

    /* find midj, such that HH[j]+RR[j] or DD[j]+SS[j]+gap is the maximum */

    midh = HH[0] + RR[0];
    midj = 0;
    type = 1;
    for (j = 0; j <= n; j++) {
      hh = HH[j] + RR[j];
      if (hh >= midh)
        if (hh > midh || (HH[j] != DD[j] && RR[j] == SS[j])) {
          midh = hh;
          midj = j;
        }
    }

    for (j = n; j >= 0; j--) {
      hh = DD[j] + SS[j] + gS[j];
      if (hh > midh) {
        midh = hh;
        midj = j;
        type = 2;
      }
    }
  }

  /* Conquer recursively around midpoint                                   */

  if (type == 1) { /* Type 1 gaps  */
    if (debug > 2) fprintf(stdout, "Type 1,1: midj %d\n", (int)midj);
    diff1(A, B, midi, midj, go1, 1);
    if (debug > 2) fprintf(stdout, "Type 1,2: midj %d\n", (int)midj);
    diff1(A + midi, B + midj, m - midi, n - midj, 1, go2);
  } else {
    if (debug > 2) fprintf(stdout, "Type 2,1: midj %d\n", (int)midj);
    diff1(A, B, midi - 1, midj, go1, 0);
    del(2);
    if (debug > 2) fprintf(stdout, "Type 2,2: midj %d\n", (int)midj);
    diff1(A + midi + 1, B + midj, m - midi - 1, n - midj, 0, go2);
  }

  return midh; /* Return the score of the best alignment */
}

/* calculate the score for opening a gap at residues A[i] and B[j]       */

inline int MM::open_penalty1(int i, int j) {
  if (!endgappenalties && (i == 0 || i == M))
    return (0);
  else
    return (g);
}

/* calculate the score for extending an existing gap at A[i] and B[j]    */

inline int MM::ext_penalty1(int i, int j) {
  if (!endgappenalties && (i == 0 || i == N))
    return (0);
  else
    return (gh);
}

/* calculate the score for a gap of length k, at residues A[i] and B[j]  */

inline int MM::gap_penalty1(int i, int j, int k) {
  int gp;
  if (k <= 0) return (0);
  if (!endgappenalties && (i == 0 || i == M)) return (0);
  gp = g + gh * k;
  return (gp);
}
/* calculate the score for opening a gap at residues A[i] and B[j]       */

inline int MM::open_penalty2(int i, int j) {
  if (!endgappenalties && (j == 0 || j == N))
    return (0);
  else
    return (g);
}

/* calculate the score for extending an existing gap at A[i] and B[j]    */

inline int MM::ext_penalty2(int i, int j) {
  if (!endgappenalties && (j == 0 || j == N))
    return (0);
  else
    return (gh);
}

/* calculate the score for a gap of length k, at residues A[i] and B[j]  */

inline int MM::gap_penalty2(int i, int j, int k) {
  int gp;
  if (k <= 0) return (0);
  if (!endgappenalties && (j == 0 || j == N)) return (0);
  gp = g + gh * k;
  return (gp);
}

int MM::diff(int A, int B, int m, int n, int go1, int go2) {
  int midi, midj, type;
  int midh;

  static int t, tl;
  int go = g, ge = gh;

  {
    static int i, j;
    static int hh, f, e, s;

    /* Boundary cases: m <= 1 or n == 0 */
    if (debug > 2)
      fprintf(stdout, "A %d B %d m %d n %d midi %d go1 %d go2 %d\n", (int)A,
              (int)B, (int)m, (int)n, (int)m / 2, (int)go1, (int)go2);

    /* if sequence B is empty....                                            */

    if (n <= 0) {
      /* if sequence A is not empty.... */

      if (m > 0) {
        /* delete residues A[1] to A[m] */

        del(m);
      }
      return (-(go + (m - 1) * ge));
      /*return(-gap_penalty1(A,B,m)); */
    }

    /* if sequence A is empty....                                            */

    if (m <= 1) {
      if (m <= 0) {
        /* insert residues B[1] to B[n] */

        add(n);
        return (-(go + (n - 1) * ge));
        /* return(-gap_penalty2(A,B,n)); */
      }

      /* if sequence A has just one residue.... */

      if (go1 == 0) /*midh =  -gap_penalty1(A+1,B+1,n); */
        midh = -(ge + go + n * ge);
      else
        midh = -(go2 + ge + go + n * ge);
      /*midh =  -gap_penalty2(A+1,B,1)-gap_penalty1(A+1,B+1,n); */

      midj = 0;
      midh = -(go1 + ge) - (go + n * ge);
      if (midh < -(go2 + ge) - (go + n * ge)) {
        midh = -(go2 + ge) - (go + n * ge);
      }
      for (j = 1; j <= n; j++) {
        /*hh = -gap_penalty1(A,B+1,j-1) + prfscore(A+1,B+j)
         * -gap_penalty1(A+1,B+j+1,n-j); */
        if (j > 1)
          hh = -(go + (j - 1) * ge) +
               scoreMat[A + 1][B + j];  // w(A+1, B+j); //w(B+j,A+1);
                                        // //prfscore(A+1, B+j);
        else
          hh = scoreMat[A + 1][B + j];  // w(A+1, B+j); //w(B+j, A+1);
                                        // //prfscore(A+1, B+j);
        if (j < n) hh += -(go + (n - j) * ge);
        /* hh = -(go+(j-1)*ge) + absprfscore(setnum, A+1, B+j) -(go+(n-j)*ge);
         */
        if (hh > midh) {
          midh = hh;
          midj = j;
        }
      }

      if (midj == 0) {
        add(n);
        del(1);
      } else {
        if (midj > 1) add(midj - 1);
        palign();
        if (midj < n) add(n - midj);
      }
      return midh;
    }

    /* Divide sequence A in half: midi */

    midi = m / 2;

    /* In a forward phase, calculate all HH[j] and HH[j] */

    HH[0] = 0;
    t = -go;
    tl = -ge;
    for (j = 1; j <= n; j++) {
      HH[j] = t = t - ge;
      DD[j] = t - go;
    }

    t = -go1;
    for (i = 1; i <= midi; i++) {
      s = HH[0];
      HH[0] = hh = t = t - ge;
      f = t - go;

      for (j = 1; j <= n; j++) {
        f = f - ge;
        if ((hh - go - ge) > f) f = hh - go - ge;
        DD[j] = DD[j] - ge;
        if (DD[j] < HH[j] - go - ge) DD[j] = HH[j] - go - ge;
        hh = s + scoreMat[A + i][B + j];  // w(A+i, B+j); //w(B+j, A+i);
                                          // //prfscore(A+i, B+j);
        if (f > hh) hh = f;
        if (DD[j] > hh) hh = DD[j];

        s = HH[j];
        HH[j] = hh;
      }
    }

    DD[0] = HH[0];

    /* In a reverse phase, calculate all RR[j] and SS[j] */

    RR[n] = 0;
    t = -go;
    for (j = n - 1; j >= 0; j--) {
      RR[j] = t = t - ge;
      SS[j] = t - go;
    }

    t = -go2;
    for (i = m - 1; i >= midi; i--) {
      s = RR[n];
      RR[n] = hh = t = t - ge;
      f = t - go;

      for (j = n - 1; j >= 0; j--) {
        f = f - ge;
        if (hh - go - ge > f) f = hh - go - ge;
        SS[j] = SS[j] - ge;
        if (RR[j] - go - ge > SS[j]) SS[j] = RR[j] - go - ge;

        hh =
            s +
            scoreMat[A + i + 1][B + j + 1];  // w(A+i+1, B+j+1); //w(B+j+1,
                                             // A+i+1); //prfscore(A+i+1, B+j+1);
        if (f > hh) hh = f;
        if (SS[j] > hh) hh = SS[j];

        s = RR[j];
        RR[j] = hh;
      }
    }
    SS[n] = RR[n];

    /* find midj, such that HH[j]+RR[j] or DD[j]+SS[j]+gap is the maximum */

    midh = HH[0] + RR[0];
    midj = 0;
    type = 1;
    for (j = 0; j <= n; j++) {
      hh = HH[j] + RR[j];
      if (hh >= midh)
        if (hh > midh || (HH[j] != DD[j] && RR[j] == SS[j])) {
          midh = hh;
          midj = j;
        }
    }

    for (j = n; j >= 0; j--) {
      hh = DD[j] + SS[j] - go;
      if (hh > midh) {
        midh = hh;
        midj = j;
        type = 2;
      }
    }
  }

  /* Conquer recursively around midpoint                                   */

  /*fprintf(stdout, "%d %d %d %d %d %d %d \n", A, midi, m, B, midj, n, type);*/

  if (type == 1) { /* Type 1 gaps  */
    if (debug > 2) fprintf(stdout, "Type 1,1: midj %d\n", (int)midj);
    diff(A, B, midi, midj, go1, go);
    if (debug > 2) fprintf(stdout, "Type 1,2: midj %d\n", (int)midj);
    diff(A + midi, B + midj, m - midi, n - midj, go, go2);
  } else {
    if (debug > 2) fprintf(stdout, "Type 2,1: midj %d\n", (int)midj);
    diff(A, B, midi - 1, midj, go1, 0);
    del(2);
    if (debug > 2) fprintf(stdout, "Type 2,2: midj %d\n", (int)midj);
    diff(A + midi + 1, B + midj, m - midi - 1, n - midj, 0, go2);
  }

  return midh; /* Return the score of the best alignment */
}
