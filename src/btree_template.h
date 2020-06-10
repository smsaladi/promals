#ifndef _tree_template__
#define _tree_template__

#include <vector>

#include "amino.h"
#include "hmm_multim.h"
#include "hmm_psipred.h"
#include "mm.h"
#include "sequences.h"
#include "sparsematrix.h"
#include "util.h"

static int debug = 1;

template <typename TNODE>
class btree {
 public:
  TNODE *root;

  // constructors
  btree();
  btree(const btree &a);
  ~btree();

  sequences myseq;

  int size;
  bool unrooted;
  char *treename;

  // array of pointers to the TNODEs
  TNODE **v;

  // read tree from a file
  void readTree(char *treefile);
  void writeTree(char *outfile);
  void writeTopology(char *outfile);

  // reroot the tree
  void reroot(TNODE *);
  // root the tree between r and its parent
  void rootTree(TNODE *, double);
  void unrootTree();
  void leastSquareRoot(TNODE **rootNode, double *dist, double *var, int &Num);

  // UPGMA, given a distance matrix
  void UPGMA(float **distMat, vector<string> seq, vector<string> name,
             int nseqs);

  // progressive alignment by prof-prof hidden markov model
  void progressiveAlignHMM(TNODE *n);

  // progressive alignment by prof-prof hidden markov model: fast stage
  void progressiveAlignHMM_FastStage(TNODE *n, float distCutoff);
  vector<int> progressiveAlignHMM_FastStage_mafft(TNODE *n, float distCutoff);
  // generate a vector of pointers to TNODES that correpond to pre-aligned
  // groups
  vector<TNODE *> preAligned;
  void obtainPreAligned(
      TNODE *);  // definition of a pre-aligned: the node itself is
                 // aligned, but its parent is not aligned
  void
  getProfileForPreAligned();  // calculate the profile for the vector preAligned
  sparseMatrix ***smat,
      ***smat1;  // sparse matrix for installing posterior probabilities
  void profileConsistency();
  void profileConsistency_local();
  void profileConsistency_glocal(float wg);
  void profileConsistency_multim(hmm_parameters *params);
  void profileConsistency_multim(hmm_parameters *params, double **dist_matrix,
                                 double max_dist_cutoff);
  void profileConsistency_multim2(hmm_parameters *params1,
                                  hmm_parameters *params2, sequences *tmpSeq2,
                                  double id_cutoff);
  // void profileConsistency_multim(char *params);
  void printAlignmentFromAbs(TNODE *n, char *outfilename,
                             vector<subalign *> similaraln,
                             vector<char *> repnames);
  void printAlignmentFromAbs(TNODE *n);
  void profileConsistency_profilehmm(hmm_parameters *params, char *ss_dir_name,
                                     int use_ss, float ss_weight,
                                     double **dist_matrix,
                                     double max_dist_cutoff);

  // consistency derived from a set of alignments
  void alignment_consistency(vector<sequences> seq_vector);

  // progressive alignment using consistency function: slow stage
  void computeConsistencyAlignment(TNODE *a);
  int **computeConsistMatrix(TNODE *a, TNODE *b);
  void relaxConsistMatrix();
  void relaxConsistMatrix(double **dist_matrix, double max_dist_cutoff);
  void relaxConsistMatrix(double **dist_matrix, double max_dist_cutoff,
                          double min_cutoff);
  int *computePairwiseAlignment(int **scoreMat, int m, int n,
                                int &len);  // by dynamic programming
  void computeConsistencyAlignment(TNODE *a, float divergent_cutoff);
  void profileConsistency_psipred(hmm_psipred_parameters *params,
                                  double **dist_matrix, double max_dist_cutoff,
                                  int use_homologs);
  void profileConsistency_psipred_sum_of_pairs(hmm_psipred_parameters *params,
                                               double **dist_matrix,
                                               double max_dist_cutoff,
                                               int use_homologs);

  static int btreecount;

  void get_descendants(TNODE *r);
  void assign_weights(TNODE *r);
  void refine_align_new(subalign *x);

  void iterativeRefinement(int maxround);
  void recomputeAlignment(TNODE *a);

 private:
  void readTreefile(ifstream &tf, TNODE *rt);
  void writeTreefile(ofstream &otf, TNODE *rt);
  void writeTopologyfile(ofstream &otf, TNODE *rt);
  void map_v();
  void map_v_recurse(TNODE *r, int &index);

  void checkRoot();
};

// below are functions of btree
template <typename TNODE>
int btree<TNODE>::btreecount = 0;

template <typename TNODE>
btree<TNODE>::btree() {
  root = 0;
  size = 0;
  v = 0;
  treename = 0;
}

template <typename TNODE>
btree<TNODE>::~btree() {
  if (v) delete[] v;
  if (treename) delete[] treename;
}

template <typename TNODE>
void btree<TNODE>::readTree(char *treefile) {
  int i, j, k;

  ifstream tf(treefile, ios::in);
  if (!tf) {
    cout << "Tree file " << treefile << " cannot be read" << endl;
    exit(0);
  }

  i = strlen(treefile);
  treename = cvector(i + 2);
  strcpy(treename, treefile);

  root = new TNODE;
  size++;
  // cout << size << endl;
  root->rootFlag = 1;
  readTreefile(tf, root);

  map_v();
}

template <typename TNODE>
void btree<TNODE>::readTreefile(ifstream &tf, TNODE *rt) {
  int i, j, k;
  char cfirst;
  char name[50];
  double branchlen;

  // ignore the spaces
  while ((tf.peek() == ' ') || (tf.peek() == '\n') || (tf.peek() == '\t')) {
    tf.get();
  }
  cfirst = tf.peek();
  //
  // cout <<"|"<< cfirst << "|"<<endl;
  //
  if (cfirst == '(') {  // the childL
    cfirst = tf.get();
    rt->childL = new TNODE;
    size++;
    rt->childL->parent = rt;
    readTreefile(tf, rt->childL);
  } else if (cfirst == ',') {  // the childR or the childT
    tf.get();
    if (!rt->childR) {
      rt->childR = new TNODE;
      size++;
      rt->childR->parent = rt;
      readTreefile(tf, rt->childR);
      // unrooted = false;
    } else {
      if (rt->childT) {
        cout << "One node contains more than three children in the tree file"
             << endl;
        exit(0);
      }
      rt->childT = new TNODE;
      size++;
      rt->childT->parent = rt;
      readTreefile(tf, rt->childT);
      // unrooted = true;
    }
  }
  /* else if(cfirst == ')') {
      tf.get();
      readTreefile(tf, rt->parent();
  } */
  else if (cfirst == ';') {
    tf.close();
    return;
  } else {
    // 1. get the node name
    while ((tf.peek() == ' ') || (tf.peek() == '\n') || (tf.peek() == '\t')) {
      tf.get();
    }

    i = 0;
    cfirst = tf.peek();
    if (cfirst == ')') {
      tf.get();
    }
    while ((tf.peek() == ' ') || (tf.peek() == '\n') || (tf.peek() == '\t')) {
      tf.get();
    }
    if (tf.peek() == ';') return;

    cfirst = tf.peek();
    while ((cfirst != '(') &&
           (cfirst != ')' && (cfirst != ',') && cfirst != ':')) {
      name[i] = tf.get();
      i++;
      if (i >= 49) {
        name[i] = '\0';
        cout << "name in the tree too long: |" << name << "|" << endl;
        exit(0);
      }
      cfirst = tf.peek();
    }
    name[i] = '\0';
    // strcpy(rt->name, name);
    rt->name = name;
    // cout << "==== " << rt->name << " ==="<<endl;

    // 2. get the node len
    if (cfirst == ':') {
      tf.get();
      tf >> rt->branchlen;
    }
    // cout << "==== " << rt->name << " === " << rt->branchlen << " === " <<
    // rt->nodecount <<endl;
    readTreefile(tf, rt->parent);
  }
}

template <typename TNODE>
void btree<TNODE>::writeTree(char *outfile) {
  int i, j, k;
  ofstream otf(outfile, ios::out);
  if (!otf) {
    cout << outfile << " cannot be written into" << endl;
    exit(0);
  }
  writeTreefile(otf, root);
  otf.close();
}

template <typename TNODE>
void btree<TNODE>::writeTreefile(ofstream &otf, TNODE *rt) {
  int i, j, k;

  // if rt is a leaf node: print leaf node name + branchlen
  if ((!rt->childL) && (!rt->childR) && (!rt->childT)) {
    otf << rt->name << ":" << rt->branchlen;
  }

  // if rt is an internal node: print "(", childrens, ")", and branchlen
  else {
    otf << "(";
    if (rt->childL) {
      writeTreefile(otf, rt->childL);
    }
    if (rt->childR) {
      otf << ",";
      writeTreefile(otf, rt->childR);
    }
    if (rt->childT) {
      otf << ",";
      writeTreefile(otf, rt->childT);
    }
    otf << ")";
    // exception is the root, where the branchlen is not printed
    if (!rt->rootFlag) otf << ":" << rt->branchlen;
    // instead, the end signal ";" is printed
    else {
      otf << ";" << endl;
    }
  }
}

template <typename TNODE>
void btree<TNODE>::writeTopology(char *outfile) {
  int i, j, k;
  ofstream otf(outfile, ios::out);
  if (!otf) {
    cout << outfile << " cannot be written into" << endl;
    exit(0);
  }
  writeTopologyfile(otf, root);
  otf.close();
}

template <typename TNODE>
void btree<TNODE>::writeTopologyfile(ofstream &otf, TNODE *rt) {
  int i, j, k;

  // if rt is a leaf node: print leaf node name + branchlen
  if ((!rt->childL) && (!rt->childR) && (!rt->childT)) {
    otf << rt->name;
    // otf << rt->name << ":" << rt->branchlen ;
  }

  // if rt is an internal node: print "(", childrens, ")", and branchlen
  else {
    otf << "(";
    if (rt->childL) {
      writeTopologyfile(otf, rt->childL);
    }
    if (rt->childR) {
      otf << ",";
      writeTopologyfile(otf, rt->childR);
    }
    if (rt->childT) {
      otf << ",";
      writeTopologyfile(otf, rt->childT);
    }
    otf << ")";
    // exception is the root, where the branchlen is not printed
    // do not write the distance
    // if(!rt->rootFlag) otf << ":" << rt->branchlen;
    // instead, the end signal ";" is printed
    if (rt->rootFlag) {
      otf << ";" << endl;
    }
  }
}

// map the array of pointers to the nodes
// this routine is called in readTree
template <typename TNODE>
void btree<TNODE>::map_v() {
  int i, j, k;
  int index;

  // before map v, check rooted tree or not
  if (root->childT)
    unrooted = true;
  else
    unrooted = false;

  // allocation of the array of pointers to the TNODEs
  if (size == 0)
    return;
  else {
    if (v) delete[] v;
    v = new TNODE *[size + 1];
  }

  for (i = 0; i <= size; i++) v[i] = 0;

  index = 1;
  map_v_recurse(root, index);
}

template <typename TNODE>
void btree<TNODE>::map_v_recurse(TNODE *r, int &index) {
  int i, j, k;

  // cout << index << endl;

  v[index] = r;
  r->n = index;
  index++;
  if (r->childL) {
    map_v_recurse(r->childL, index);
  }
  if (r->childR) {
    map_v_recurse(r->childR, index);
  }
  if (r->childT) {
    map_v_recurse(r->childT, index);
  }
}

// change the root to another TNODE
// seems easy; but actually require some design
// recursion might be easier to write
template <typename TNODE>
void btree<TNODE>::reroot(TNODE *nr) {
  int i, j, k;
  TNODE *tmpn = nr;
  TNODE *tmpp = nr->parent;
  TNODE *tmppo;
  double tmpbranchlen, tmpbranchlen1;

  checkRoot();

  // for a rooted tree, make the root disppear
  TNODE *tmproot = root;
  if (root->childT == 0) {
    if (root->childL->childL) {  // left child is not a leaf node
      root->childL->childT = root->childR;
      root->childR->parent = root->childL;
      root->childR->branchlen += root->childL->branchlen;
      root->childL->rootFlag = 1;
      root->childL->parent = 0;
      root = root->childL;
      root->branchlen = 0;
      delete tmproot;
      unrooted = true;
      size--;
    } else {
      root->childR->childT = root->childL;
      root->childL->parent = root->childR;
      root->childL->branchlen += root->childR->branchlen;
      root->childR->rootFlag = 1;
      root->childR->parent = 0;
      root = root->childR;
      root->branchlen = 0;
      delete tmproot;
      unrooted = true;
      size--;
    }
  }

  if (nr == root) {
    // cout << "reroot: already at the root\n";
    return;
  }

  // nr should not be a leaf node
  if (nr->childL == 0) {
    cout << "reroot error msg: nr is a leaf node and cannot be the root"
         << endl;
    return;
  }

  while (tmpn != root) {
    if (tmpn == nr) {  // if tmpn is the new root; make the third child
      // 1. make the n's childT point to n's parent: this is my convention
      //    keep the childR and childL unchanged
      tmpn->childT = tmpn->parent;

      // judge if tmpp is already the root
      if (tmpp->parent) {  // not the root yet
        tmppo = tmpp->parent;
        // 2. make p's parent point to n
        tmpp->parent = tmpn;

        // 3. make p's childL or childR to be NULL
        if (tmpp->childL == tmpn)
          tmpp->childL = 0;
        else
          tmpp->childR = 0;

        // 4. update the branch length
        tmpbranchlen = tmpp->branchlen;
        tmpp->branchlen = tmpn->branchlen;
        tmpn->branchlen = 0;

        // 5. to an upper level: tmpn is changed to tmpp and tmpp is changed to
        // tmpp's orignial parent
        tmpn = tmpp;
        tmpp = tmppo;
      } else {  // tmpp is already the root; make the third child be NULL
        tmpp->parent = tmpn;
        if (tmpp->childT == tmpn) {
          tmpp->childT = 0;
        } else if (tmpp->childR == tmpn) {
          tmpp->childR = tmpp->childT;
          tmpp->childT = 0;
        } else {
          tmpp->childL = tmpp->childT;
          tmpp->childT = 0;
        }
        // tmpbranchlen = tmpp->branchlen;
        tmpp->branchlen = tmpn->branchlen;
        tmpn->branchlen = 0;
        tmpn = tmpp;  // the loop ends here since tmpp is the root
      }
    } else {
      // 1. make n' childR or childL point to p
      if (tmpn->childL) {
        tmpn->childR = tmpp;
      } else {
        tmpn->childL = tmpp;
      }

      // judge if tmpp is already the root
      if (tmpp->parent) {  // not the root yet
        tmppo = tmpp->parent;
        // 2. make p's parent point to n
        tmpp->parent = tmpn;

        // 3. make p's childL or childR to be NULL
        if (tmpp->childL == tmpn)
          tmpp->childL = 0;
        else
          tmpp->childR = 0;

        // 4. update the branch length
        tmpbranchlen1 = tmpp->branchlen;
        tmpp->branchlen = tmpbranchlen;
        tmpbranchlen = tmpbranchlen1;

        // 5. to an upper level: tmpn is changed to tmpp and tmpp is changed to
        // tmpp's orignial parent
        tmpn = tmpp;
        tmpp = tmppo;
      } else {  // tmpp is already the root; make the third child be NULL
        tmpp->parent = tmpn;
        if (tmpp->childT == tmpn) {
          tmpp->childT = 0;
        } else if (tmpp->childR == tmpn) {
          tmpp->childR = tmpp->childT;
          tmpp->childT = 0;
        } else {
          tmpp->childL = tmpp->childT;
          tmpp->childT = 0;
        }
        // exchange tmpp->branchlen and tmpbranchlen
        tmpbranchlen1 = tmpp->branchlen;
        tmpp->branchlen = tmpbranchlen;
        tmpbranchlen = tmpbranchlen1;
        tmpn = tmpp;  // the loop ends here since tmpp is the root
      }
    }
  }  // end of while

  // update the root information
  root->rootFlag = 0;  // old root
  root = nr;           // new root
  root->rootFlag = 1;  // new root
  root->parent = 0;    // the parent of the root is NULL
  root->branchlen = 0;

  // update the v array
  map_v();

  checkRoot();
}

// root the tree by adding a new node between r and r->parent
// the distance from the new root to r is len
template <typename TNODE>
void btree<TNODE>::rootTree(TNODE *r, double len) {
  int i, j, k;
  TNODE *rp;

  // cout << "----\n";
  // cout << root->branchlen << endl;
  checkRoot();
  // cout << "++++\n";

  // r cannot be a leaf node or the root
  /* if(r->childL==0) {
          cout << "rootTree error msg: r cannot be a leaf node" << endl;
          cout << "rootTree not executed" << endl;
          return;
  }
  else */
  if (r->childT) {
    cout << "rootTree error msg: r cannot be the root" << endl;
    cout << "rootTree not executed" << endl;
    return;
  }

  // len should not be longer than r->branchlen
  if (len > r->branchlen) {
    cout << "rootTree error msg: len larger than r->branchlen" << len << "\t"
         << r->branchlen << endl;
    cout << "rootTree not execuated" << endl;
    return;
  }

  // reroot the tree to rp
  rp = r->parent;
  reroot(rp);

  // add the new root
  TNODE *newRoot = new TNODE;
  if (rp->childR == r) {
    rp->childR = rp->childT;
  } else if (rp->childL == r) {
    rp->childL = rp->childT;
  }
  rp->childT = 0;
  rp->parent = newRoot;
  r->parent = newRoot;
  newRoot->childR = r;
  newRoot->childL = rp;
  newRoot->childT = 0;
  newRoot->branchlen = 0;
  rp->branchlen = r->branchlen - len;
  r->branchlen = len;
  size++;

  /*
          // reroot the tree to r
          rp = r->parent;
          reroot(r); // after rerooting, rp becomes r->childT (my convention)

          // add the new root
          TNODE *newRoot = new TNODE;
          newRoot->childL = r;
          newRoot->childR = rp;
          newRoot->childT = 0;
          r->parent = newRoot;
          r->childT = 0;
          rp->parent = newRoot;
          r->branchlen = len;
          rp->branchlen = rp->branchlen - len;
          size++;
  */

  // update the root information
  root->rootFlag = 0;
  root = newRoot;
  root->rootFlag = 1;
  root->parent = 0;

  // update the v array
  map_v();
  // TNODE **v1 = new TNODE * [size+1];
  // for(i=1;i<size;i++) v1[i] = v[i];
  // v1[size] = root;
  // root->n = size;
  // delete [] v;
  // v = new tonde[size+1];
  // for(i=1;i<=size;i++) v[i] = v1[i];
  // delete [] v1;

  checkRoot();

  return;
}

// make a rooted tree unrooted
template <typename TNODE>
void btree<TNODE>::unrootTree() {
  int i, j, k;
  TNODE *r = root->childL;  // the new root

  checkRoot();

  if (root->childT) {  // already an unrooted tree
    return;
  }

  if (root->childL->childL) {
    root->childR->parent = root->childL;
    root->childL->childT = root->childR;
    root->childR->branchlen += root->childL->branchlen;
    root->childL->branchlen = 0;
    root->childR = root->childL = 0;
    root->rootFlag = 0;
  } else {
    root->childL->parent = root->childR;
    root->childR->childT = root->childL;
    root->childL->branchlen += root->childR->branchlen;
    root->childR->branchlen = 0;
    root->childL = root->childR = 0;
    root->rootFlag = 0;
  }

  delete root;
  root = r;
  root->rootFlag = 1;
  root->parent = 0;
  size--;
  unrooted = true;

  map_v();

  checkRoot();

  return;
}

template <typename TNODE>
void btree<TNODE>::checkRoot() {
  if (root->parent) {
    cout << "Root has a non-null parent" << endl;
    exit(1);
  }
  if (root->branchlen != 0) {
    cout << "Root has a non-zero branchlen" << endl;
    exit(1);
  }
  if (unrooted && (!root->childT)) {
    cout << "Unrooted tree should have childT" << endl;
    exit(1);
  }
  if ((!unrooted) && (root->childT)) {
    cout << "A rooted tree should not have childT" << endl;
    exit(1);
  }
}

// find in each branch if there exists a point other than the end points that
// minimizes the variance of the distances of all leaf nodes to that point Na:
// number of leaf node in subtree a; Nb: number of leaf node in subtree b a[],
// b[]: distance from the leaf node to the root a[i] + x and b[j] - x

//    1 \ 			      / 3
//	 \  _root		     /
//	  \/________._______________/
//	  /			    \                                     .
//	 /			     \____ 4
//    2 /  |--- x --|		      \                                   .
// 				       5
//      |a |----------b------------------|

//	Na = 2;  Nb = 3

// average distance from the leaf nodes to the root is:
//     ave_dist = t + s * x
//   t = (sum{1}^{Na}(a[i]) + sum{1}^{Nb}(b[j]) )/ (Na+Nb)
//   s = (Na - Nb)/(Na + Nb)

// variance * (Na + Nb) = sum{1}^{Na}(a[i]+x-t-s*x) + sum{1}^{Nb}(b[j]-x-t-s*x)

// derivative(variance) = 0 to get the x(min) = argmin(variance)

template <typename TNODE>
void btree<TNODE>::leastSquareRoot(TNODE **rootNode, double *dist, double *var,
                                   int &Num) {
  int i, j, k;
  TNODE **origV;
  TNODE *origR;
  int origSize;

  double *a, *b;
  int Na, Nb;
  double t, s;
  double x;

  // if it is a rooted tree, unroot
  if (!unrooted) {
    unrootTree();
  }
  origSize = size;

  // record the orignial vector of TNODEs in the tree and the root position
  origV = new TNODE *[size + 1];
  for (i = 1; i <= size; i++) origV[i] = v[i];
  origR = root;
  a = dvector(size + 1);
  b = dvector(size + 1);

  k = 1;
  for (i = 1; i <= origSize; i++) {
    var[k] = 0;

    for (j = 1; j <= origSize + 1; j++) {
      a[j] = 0;
      b[j] = 0;
    }
    Na = Nb = 0;
    t = s = 0;
    x = 0;

    if (origV[i] == origR) {
      continue;
    }  // ignore the root

    // root the tree to origV[i]
    // cout << "i: " << i << endl;
    rootTree(origV[i], origV[i]->branchlen / 2);
    // char treename[20];
    // sprintf(treename, "t%d.tre", i);
    // writeTree(treename);

    // calculate the distances from the leaf nodes to the root
    for (j = 1; j <= size; j++) {
      // cout << j << "\t" << a[j] << "\t" << b[j] << endl;
    }

    cal_len2root(root->childL, a);
    cal_len2root(root->childR, b);
    for (j = 1; j <= size; j++) {
      if (v[j]->childL) {  // internal node ignore
        a[j] = 0;
        b[j] = 0;
        continue;
      }
      // judge if j belong to a or b
      TNODE *tmpNode = v[j];
      while (tmpNode != root) {
        if (tmpNode == root->childL) {
          Na++;
          t += a[j];
          break;
        } else if (tmpNode == root->childR) {
          Nb++;
          b[j] += (root->childL->branchlen + root->childR->branchlen);
          t += b[j];
          break;
        } else {
          tmpNode = tmpNode->parent;
        }
      }
      // cout << j << "\t" << a[j] << "\t" << b[j] << endl;
      /* if(a[j]!=0) {
              Na++;
              t += a[j];
      }
      if(b[j]!=0) {
              Nb++;
              b[j] += (root->childL->branchlen + root->childR->branchlen);
              t += b[j];
      } */
    }
    t /= (Na + Nb);
    s = 1.0 * (Na - Nb) / (Na + Nb);

    // cout << "Na: " << Na << " Nb: " << Nb << " t: " << t << " s: " << s <<
    // endl;

    for (j = 1; j <= size; j++) {
      if (v[j]->childL) continue;
      // judge if j belong to a or b
      TNODE *tmpNode = v[j];
      while (tmpNode != root) {
        if (tmpNode == root->childL) {
          x += (a[j] - t) * (1 - s);
          // cout << "j: a " << j << " " << (a[j]-t)*(1-s) << endl;
          break;
        } else if (tmpNode == root->childR) {
          x -= (b[j] - t) * (1 + s);
          // cout << "j: b " << j << " " << (b[j]-t)*(1+s) << endl;
          break;
        } else {
          tmpNode = tmpNode->parent;
        }
      }
    }
    double y = (1 + s) * (1 + s) * Nb;
    y += (1 - s) * (1 - s) * Na;
    y = 0 - y;
    x /= y;
    // cout << "x: " << x << " " << root->childL->branchlen +
    // root->childR->branchlen << endl; x / = (0-(1+s)*(1+s)*Nb-(1-s)*(1-s)*Na );

    if ((x >= 0) && x < (root->childL->branchlen + root->childR->branchlen)) {
      rootNode[k] = root->childL;
      dist[k] = x;
      for (j = 1; j <= size; j++) {
        if (a[j] != 0) {
          var[k] += ((a[j] - t) + (1 - s) * dist[k]) *
                    ((a[j] - t) + (1 - s) * dist[k]);
        }
        if (b[j] != 0) {
          var[k] += ((b[j] - t) - (1 + s) * dist[k]) *
                    ((b[j] - t) - (1 + s) * dist[k]);
        }
      }
      cout << k << "\t" << dist[k] << "\t" << var[k] << endl;
      root->childR->branchlen =
          root->childL->branchlen + root->childL->branchlen - x;
      root->childL->branchlen = x;
      char outtreename[20];
      sprintf(outtreename, "%s%d.tre", treename, k);
      writeTree(outtreename);
      k++;
    }

    // cout << "i again: " << i << endl;
    reroot(origR);
    // cout << root->branchlen << endl;
    // cout << "=========" << endl;

  }  // end of i
  Num = k - 1;
}

template <typename TNODE>
void btree<TNODE>::UPGMA(float **distMat, vector<string> seq,
                         vector<string> name, int nseqs) {
  int i, j, k;
  double minElement;
  int ci = 0, cj = 0;
  // int mark[nseqs+1];
  // double d[nseqs+1];

  int *mark;
  double *d;
  mark = ivector(nseqs);
  d = dvector(nseqs);

  for (i = 1; i <= nseqs; i++) {
    mark[i] = 1;
    d[i] = 0;
  }

  // allocate an array of tnode and assign the names and seqs
  TNODE **tnodeArray = new TNODE *[nseqs * 2];  // new
  v = tnodeArray;                               // new
  for (i = 1; i <= 2 * nseqs - 1; i++) {
    tnodeArray[i] = new TNODE;
  }
  size = 2 * nseqs - 1;
  for (i = 1; i <= nseqs; i++) {
    // copy the names
    // if(name[i].length()>=49) name[i].copy(tnodeArray[i]->name,49);
    // else name[i].copy(tnodeArray[i]->name, name[i].length());
    tnodeArray[i]->name = name[i];
    tnodeArray[i]->n = i;  // new
    // cout << tnodeArray[i]->name << endl;

    // copy the arrays: skipping the gaps
    tnodeArray[i]->aseq = cvector(seq[i].length() + 1);
    int p = -1;
    for (j = 0; j < seq[i].length(); j++) {
      if (isalpha(seq[i][j])) {
        p++;  // skipping the gaps
        tnodeArray[i]->aseq[p] = seq[i][j];
      }
    }
    tnodeArray[i]->aseq[p + 1] = '\0';
    if (debug > 1)
      cout << "TNODE sequence: " << tnodeArray[i]->name << "\t"
           << tnodeArray[i]->aseq << endl;
    if (debug > 1) cout << seq[i] << endl;

    // get the subalign from the seq and name
    tnodeArray[i]->getSubalign();
    // cout << i << " " << tnodeArray[i]->aligned << endl;
  }

  tnodeArray = new TNODE *[nseqs + 1];
  for (i = 1; i <= nseqs; i++) tnodeArray[i] = v[i];

  double **newMat = dmatrix(nseqs, nseqs);
  for (i = 1; i <= nseqs; i++)
    for (j = 1; j <= nseqs; j++) newMat[i][j] = distMat[i][j];

  if (nseqs == 1) {
    root = tnodeArray[1];
    tnodeArray[1]->rootFlag = true;
  }

  int current_tnode_index = nseqs + 1;  // new
  for (k = 1; k < nseqs; k++) {
    // find the minimum distance
    minElement = 100000;
    for (i = 1; i <= nseqs; i++) {
      if (mark[i] == 0) continue;
      for (j = 1; j < i; j++) {
        if (mark[j] == 0) continue;
        if (newMat[i][j] < minElement) {
          minElement = newMat[i][j];
          ci = i;
          cj = j;
        }
      }
    }
    // cout << "minElement: " << minElement << endl;
    // TNODE *a = new TNODE;
    TNODE *a = v[current_tnode_index];  // new
    a->n = current_tnode_index;         // new
    // cout << "current_tnode_index " << current_tnode_index << endl;
    current_tnode_index += 1;  // new
    a->childL = tnodeArray[ci];
    a->childR = tnodeArray[cj];
    tnodeArray[ci]->parent = a;
    tnodeArray[cj]->parent = a;
    tnodeArray[ci]->branchlen =
        minElement / 2 - d[ci];  //(minElement + d[cj] - d[ci])/2;
    tnodeArray[cj]->branchlen =
        minElement / 2 - d[cj];  //(minElement + d[ci] - d[cj])/2;
    // cout << tnodeArray[ci]->branchlen << "\t" << tnodeArray[cj]->branchlen <<
    // endl;
    for (i = 1; i <= nseqs; i++) {
      if (mark[i] == 0) continue;
      if (i == cj) continue;
      if (i == ci) continue;
      newMat[ci][i] = newMat[i][ci] =
          (newMat[ci][i] * mark[ci] + newMat[cj][i] * mark[cj]) /
          (mark[ci] + mark[cj]);
    }
    mark[ci] = mark[ci] + mark[cj];
    mark[cj] = 0;
    tnodeArray[ci] = a;
    d[ci] = minElement / 2;  //(minElement + d[cj] + d[ci])/2;
    if (k == (nseqs - 1)) {
      a->rootFlag = true;
      root = a;
    }
  }

  // cout << "finished here" << endl;
  // writeTree("tmptree.tre");
  //
  // for(i=1;i<nseqs*2;i++) { cout << i << " " << v[i]->aligned << endl; }

  delete[] mark;
  delete[] d;
  delete[] tnodeArray;
  free_dmatrix(newMat, nseqs, nseqs);

  // writeTree("tmp.tree");
}

#endif
