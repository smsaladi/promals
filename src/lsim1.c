#include <stdio.h>
#include <math.h>
#include "pcma.h" 
/*#include "new.h"*/

extern SN ****glib;
extern streeptr *groupptr;
extern char *am;
extern sint debug;

/* extern char name0[], name1[]; */
/* extern int match, mismh; */
/* extern char *sq, sqnam[], *seqc0, *seqc1; */
/* extern char ttitle[], ltitle[]; */
int min0,min1,max0,max1;
int smin0, smin1;
int markx = 0;
int gscore;

#define min(x,y) ((x)<=(y) ? (x) : (y))


static int **vv;			/* substitution scores */
static int q, r;			/* gap penalties */
static int qr;				/* qr = q + r */

#ifdef FAR_PTR
typedef struct ONE
	{ int COL ;  struct ONE  far * NEXT ;}
	pair, far * pairptr;
pairptr far *row, z;		/* for saving used aligned pairs */
#define PAIRNULL (pairptr)NULL
#else
typedef struct ONE { int COL ;  struct ONE  *NEXT ;} pair, *pairptr;
pairptr *row, z; 			/* for saving used aligned pairs */
#define PAIRNULL (pairptr)NULL
#endif
static int tt;

typedef struct LNODE
	{ int  SCORE;
	  int  STARI;
	  int  STARJ;
	  int  ENDI;
	  int  ENDJ;
	  int  TOP;
	  int  BOT;
	  int  LLEFT;
	  int  LRIGHT; }  vertex,
#ifdef FAR_PTR
 far *vertexptr;
#else
     *vertexptr;
#endif
		
vertexptr  *LIST;			/* an array for saving k best scores */
vertexptr  low = 0;			/* lowest score node in LIST */
vertexptr  most = 0;			/* latestly accessed node in LIST */
static int numnode;			/* the number of nodes in LIST */

static int *CC, *DD;			/* saving matrix scores */
static int *RR, *SS, *EE, *FF; 		/* saving start-points */
static int *HH, *WW;		 	/* saving matrix scores */
static int *II, *JJ, *XX, *YY; 		/* saving start-points */
static int  m1, mm, n1, nn;		/* boundaries of recomputed area */
static int  rl, cl;			/* left and top boundaries */
static int  lmin;			/* minimum score in LIST */
static int flag;			/* indicate if recomputation necessary*/

/* DIAG() assigns value to x if (ii,jj) is never used before */
#define DIAG(ii, jj, x, value)				\
{ for ( tt = 1, z = row[(ii)]; z != PAIRNULL; z = z->NEXT )	\
    if ( z->COL == (jj) )				\
      { tt = 0; break; }				\
  if ( tt )						\
    x = ( value );					\
}

/* replace (ss1, xx1, yy1) by (ss2, xx2, yy2) if the latter is large */
#define ORDER(ss1, xx1, yy1, ss2, xx2, yy2)		\
{ if ( ss1 < ss2 )					\
    { ss1 = ss2; xx1 = xx2; yy1 = yy2; }		\
  else							\
    if ( ss1 == ss2 )					\
      { if ( xx1 < xx2 )				\
	  { xx1 = xx2; yy1 = yy2; }			\
	else						\
	  if ( xx1 == xx2 && yy1 < yy2 )		\
	    yy1 = yy2;					\
      }							\
}

/* The following definitions are for function diff() */

int  diff(), display();
static int  zero = 0;				/* int type zero        */

#define gap(k)  ((k) <= 0 ? 0 : q+r*(k))	/* k-symbol indel score */

static int *sapp;				/* Current script append ptr */
static int  last;				/* Last script op appended */

static int I, J;				/* current positions of A ,B */
static int no_mat; 				/* number of matches */ 
static int no_mis; 				/* number of mismatches */ 
static int al_len; 				/* length of alignment */
						/* Append "Delete k" op */
#define DEL(k)				\
{ I += k;				\
  al_len += k;				\
  if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
}
						/* Append "Insert k" op */
#define INS(k)				\
{ J += k;				\
  al_len += k;				\
  if (last < 0)				\
    { sapp[-1] = (k); *sapp++ = last; }	\
  else					\
    last = *sapp++ = (k);		\
}

						/* Append "Replace" op */
#define REP 				\
{ last = *sapp++ = 0; 			\
  al_len += 1;				\
}

#ifndef FAR_PTR
#define FCKALLOC lckalloc
#else
#define FCKALLOC flckalloc
#endif

/* SIM(A,B,M,N,K,V,Q,R) reports K best non-intersecting alignments of
   the segments of A and B in order of similarity scores, where
   V[a][b] is the score of aligning a and b, and -(Q+R*i) is the score
   of an i-symbol indel.  						*/

void SIM(int M,int N,int K,int **V,int Q,int R,int nseq,int gi, int gj)
{
  int endi, endj, stari, starj;	/* endpoint and startpoint */ 
  int  score;   			/* the max score in LIST */
  int count;				/* maximum size of list */	
  register  int  i, j;			/* row and column indices */
  char *lckalloc();			/* space-allocating function */
#ifdef FAR_PTR
  char far *flckalloc();
#endif
  int  *S;				/* saving operations for diff */
  int nc, nd, ns, nident;		/* for display */
  int tmp;				/* for switching min0,min1 */
  vertexptr cur; 			/* temporary pointer */
  vertexptr findmax();	 		/* return the largest score node */
  double percent;

  /* JP */
  int first, second;
  SN **snd;
  int *aln_path1, *aln_path2;
  int tmpfirst, tmpsecond;
  int n, nblocks;
  int ii, jj;
  int tpairs, idpairs;
  int k;
  int ngaps, nonegaps;
	
  /* allocate space for all vectors */
  j = (N + 1) * sizeof(int);
  CC = ( int * ) lckalloc(j);
  DD = ( int * ) lckalloc(j);
  RR = ( int * ) lckalloc(j);
  SS = ( int * ) lckalloc(j);
  EE = ( int * ) lckalloc(j);
  FF = ( int * ) lckalloc(j);
  i = (M + 1) * sizeof(int);
  HH = ( int * ) lckalloc(i);
  WW = ( int * ) lckalloc(i);
  II = ( int * ) lckalloc(i);
  JJ = ( int * ) lckalloc(i);
  XX = ( int * ) lckalloc(i);
  YY = ( int * ) lckalloc(i);
  S = ( int * ) lckalloc(min(i,j)*5/4);
#ifdef FAR_PTR
  row = ( pairptr far * ) FCKALLOC( (M + 1) * sizeof(pairptr));
#else
  row = ( pairptr * ) lckalloc( (M + 1) * sizeof(pairptr));
#endif

  /* set up list for each row */
  if (nseq == 2) for ( i = 1; i <= M; i++ ) row[i]= PAIRNULL;
  else {
	  z = ( pairptr )FCKALLOC((int)sizeof(pair)*M);
	  for ( i = 1; i <= M; i++,z++) {
		  row[i] = z;
		  z->COL = i;			
		  z->NEXT = PAIRNULL;
	  }
  }

  vv = V;
  /*for(i=1;i<=M;i++) {
     for(j=1;j<=N;j++) {
  	fprintf(stdout, "%d ", vv[i][j]);
     }
     fprintf(stdout, "\n");
  }*/
  q = Q;
  r = R;
  qr = q + r;

  /*LIST = ( vertexptr * ) lckalloc( K * sizeof(vertexptr));*/
  LIST = ckalloc(K*sizeof(vertexptr));
  for ( i = 0; i < K ; i++ )
    /*LIST[i] = ( vertexptr ) FCKALLOC( (int) sizeof(vertex));*/
    LIST[i] = ckalloc(sizeof(vertex));
  
  numnode = lmin = 0;
  big_pass(M,N,K,nseq);
  
  /* Report the K best alignments one by one. After each alignment is
     output, recompute part of the matrix. First determine the size
     of the area to be recomputed, then do the recomputation         */
  
  for ( count = K - 1; count >= 0 ; count-- )
    { if ( numnode == 0 )
	lfatal("The number of alignments computed is too large");
      cur = findmax();	/* Return a pointer to a node with max score*/
      score = cur->SCORE;
      stari = ++cur->STARI;
      starj = ++cur->STARJ;
      endi = cur->ENDI;
      endj = cur->ENDJ;
      m1 = cur->TOP;
      mm = cur->BOT;
      n1 = cur->LLEFT;
      nn = cur->LRIGHT;
      rl = endi - stari + 1;
      cl = endj - starj + 1;
	if(debug>1){fprintf(stdout, "score: %d\n", score); 
	   fprintf(stdout, "rl: %d; stari: %d; endi: %d\n", rl, stari, endi);
	   fprintf(stdout, "cl: %d; starj: %d; endj: %d\n", cl, starj, endj);
	}
	
      /*I = stari - 1;
      J = starj - 1; */
      I = stari-1;
      J = starj-1;
      sapp = S;
      last = 0;
      al_len = 0;
      no_mat = 0;
      no_mis = 0;
      /*diff(&A[stari]-1, &B[starj]-1,rl,cl,q,q);*/
	/*fprintf(stdout, "%d %d %d %d \n", stari, starj, rl, cl);*/
	fflush(stdout); /*rl--; cl--;*/
	diff(stari-1, starj-1, rl,cl,q,q);
	/*fprintf(stdout, "============\n");fflush(stdout);*/

      min0 = stari;
      min1 = starj;
      max0 = stari+rl-1;
      max1 = starj+cl-1;
	/*for(i=0;i<rl;i++) {
	   fprintf(stdout, "%d ", *(S+i));} */


      fflush(stdout);

      /* JP: determine the alignment paths */
      aln_path1 = ckalloc( (al_len+1) * sizeof(int));
      aln_path2 = ckalloc( (al_len+1) *sizeof(int) );
      first = stari; second = starj;
      ii=0;jj=0;
      for(i=0;i<rl,(ii<al_len&&jj<al_len);i++) {
        if(*(S+i)==0) {
            ii++;jj++;
            aln_path1[ii] = 0;
            aln_path2[jj] = 0;
        }
        if(*(S+i)>0) {
            for(j=1;j<=*(S+i);j++) {
                jj++;ii++;
                aln_path2[jj]=0;
                aln_path1[ii] = 1;
            }
        }
        if(*(S+i)<0) {
            for(j=1;j<=(-*(S+i));j++) {
                ii++;jj++;
                aln_path2[jj]=1;
                aln_path1[ii]=0;
            }
        }
      }

      /* JP: calculate the average sequence identity between the two blocks */
      tpairs = 0; idpairs = 0;
      first = stari; second = starj;
      for(i=1;i<=al_len;i++) {

	if(aln_path1[i]==1) { second++; continue; }
	if(aln_path2[i]==1) { first++; continue; }

	
	for(j=1;j<=groupptr[gi]->seqnum;j++) {
	  for(k=1;k<=groupptr[gj]->seqnum;k++) {

		if(groupptr[gi]->seq[j][first] && groupptr[gj]->seq[k][second]) {
		    tpairs++;
		    if(am[groupptr[gi]->seq[j][first]]==am[groupptr[gj]->seq[k][second]]) 
			idpairs++;
		}
	  }
	}
	first++; second++;
      }


      /* JP calculate the gap fraction */
      ngaps = 0; nonegaps = 0;
      for(i=1;i<=al_len;i++) {
	if(aln_path1[i]==1) {ngaps++; continue; }
	if(aln_path2[i]==1) {ngaps++; continue; }
	nonegaps++;
      }

      if(debug>1) fprintf(stdout, "gaps %d nonegaps %d", ngaps, nonegaps);
		


      /* JP NOTICE: normalized score or not */
      /*score = score/al_len;  */ /* using average score */
      /*score = (int) (100.0*idpairs/tpairs); */
      /*score = (int) (score/sqrt(al_len) );*/
      score = endi - stari + 1;
      if(score > endj -starj +1) score = endj - starj + 1;
      score =  (int) ( (log(score) * 100.0*idpairs/tpairs)*(1.0*nonegaps/(nonegaps+ngaps)) ) ;
	
            
      /* JP: lib generation  */
      first = stari; second = starj;
      for(i=0;i<rl, first<=endi;i++) {
	if(*(S+i)==0) {
	   /* fprintf(stdout, "first: %d; second: %d\n", first, second); */

	   if(glib[gi][gj][first] == NULL) {
  	     glib[gi][gj][first] = SNavail();
	     glib[gi][gj][first]->ind = second;
	     glib[gi][gj][first]->sbe = score;
	   }
	   else {
	     snd = &glib[gi][gj][first];
	     /*fprintf(stdout, "snd: %d; glib:  %d\n", snd, glib[gi][gj][first]);*/
	     while(*snd) {
		 if((*snd)->ind == second) {
		    (*snd)->sbe += score;
		    break;
		  }
		  else {
		     snd = &((*snd)->next);
		     /*fprintf(stdout, "snd: %d\n", snd);*/
		  }
	     }
	     if(!(*snd)) {
		(*snd) = SNavail();
		(*snd)->ind = second;
		(*snd)->sbe = score;
		/*fprintf(stdout, "====%d %d %d \n", first, (*snd)->ind, (*snd)->sbe);*/
	     }
	   }
	   /* AddSbe(glib[gi][gj][first], second, score/rl);*/
	   if(glib[gj][gi][second]==NULL) {
	     glib[gj][gi][second] = SNavail();
	     glib[gj][gi][second]->ind = first;
	     glib[gj][gi][second]->sbe = score;
	   }
	   else {
             snd = &glib[gj][gi][second];
             while(*snd) {
                 if((*snd)->ind == first) {
                    (*snd)->sbe += score;
                    break;
                  }
                  else snd = &((*snd)->next);
             }
             if(!*snd) {
                *snd = SNavail();
                (*snd)->ind = first;
                (*snd)->sbe = score;
             }
           }
	   /* AddSbe(glib[gj][gi][second], first, score/rl); */
	   first++; second++;
	}
	if(*(S+i) > 0) {
	   second+=*(S+i);
	   /*first++; second++; */
	}
	if(*(S+i) < 0) {
	   first-=*(S+i);
	   /*first++; second++; */
	}
      }

      /* JP: print out the subalignments */
   if(debug > 1) {

      first = stari;second = starj;
      nblocks = (al_len-1)/80+1;
      fprintf(stdout, "subalignment %d score %d average score %d average id %d\n", K-count, score, score/al_len, (int) (100.0*idpairs/tpairs) );
      fprintf(stdout, "tpairs %d idpairs %d \n", tpairs, idpairs); 
      for(n=0;n<nblocks;n++) {
	fflush(stdout);
	fprintf(stdout, "first: %d\n", first);
	for(j=1;j<=groupptr[gi]->seqnum;j++) {
	    tmpfirst = first;
	    fprintf(stdout, "%s\t", groupptr[gi]->name[j]);
	    for(i=n*80+1;(i<=(n+1)*80&&i<=al_len);i++){
		if(aln_path1[i]==0) {
		   fprintf(stdout, "%c", am[groupptr[gi]->seq[j][tmpfirst]]);
		   tmpfirst++;
		}
		else {fprintf(stdout, "-");}
	    }
	    fprintf(stdout, "\n");
	}
	first = tmpfirst;
	
	fprintf(stdout, "second: %d\n", second);
	for(j=1;j<=groupptr[gj]->seqnum;j++){
	    tmpsecond = second;
	    fprintf(stdout, "%s\t", groupptr[gj]->name[j]);
	    for(i=n*80+1;(i<=(n+1)*80&&i<=al_len);i++) {
		if(aln_path2[i]==0) {
		   fprintf(stdout, "%c", am[groupptr[gj]->seq[j][tmpsecond]]);
		   tmpsecond++;
		}
		else{fprintf(stdout,"-");}
	    }
	    fprintf(stdout, "\n");
	}
	second = tmpsecond;
	fprintf(stdout, "\n\n");
      }
   }

      ckfree(aln_path1);
      ckfree(aln_path2);

      /*first = i;
      for(i=0;i<first;i++) {
           fprintf(stdout, "%d ", *(S+i));}
	fprintf(stdout, "\n");*/
	
      
      if ( count )
	{ flag = 0;
	  locate(nseq);
	  if ( flag )
	    small_pass(count,nseq);
	}
    }

  /* JP: reinitiate and free memories */
  low = 0; most = 0;
  for ( i = 0; i < K ; i++ )
    LIST[i] = ckfree(LIST[i]);
  LIST = ckfree(LIST);

  ckfree(CC);
  ckfree(DD);
  ckfree(RR);
  ckfree(SS);
  ckfree(EE);
  ckfree(FF);
  ckfree(HH);
  ckfree(WW);
  ckfree(II);
  ckfree(JJ);
  ckfree(XX);
  ckfree(YY);
  ckfree(S);
  ckfree(row);

  
}

/* A big pass to compute K best classes */

big_pass(M,N,K,nseq) int M,N,K,nseq;
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */ 
  register  int  di, dj;		/* end-point associated with d */
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  int  *va;				/* pointer to vv(A[i], B[j]) */
  int   addnode();			/* function for inserting a node */

	
	/* Compute the matrix and save the top K best scores in LIST
	   CC : the scores of the current row
	   RR and EE : the starting point that leads to score CC
	   DD : the scores of the current row, ending with deletion
	   SS and FF : the starting point that leads to score DD        */
 	/* Initialize the 0 th row */
	for ( j = 1; j <= N ; j++ )
	  {  CC[j] = 0;
	     RR[j] = 0;
	     EE[j] = j;
	     DD[j] = - (q);
	     SS[j] = 0;
	     FF[j] = j;
	  }
	for ( i = 1; i <= M; i++) 
	  {  c = 0;				/* Initialize column 0 */
	     f = - (q);
	     ci = fi = i;
 	     /*va = vv[A[i]];*/
	     if ( nseq == 2 )
	       { p = 0;
	         pi = i - 1;
	         cj = fj = pj = 0;
	       }
	     else
	       { p = CC[i];
		 pi = RR[i];
		 pj = EE[i];
	         cj = fj = i;
	       }
	     for ( j = (nseq == 2 ? 1 : (i+1)) ; j <= N ; j++ )  
	       {  f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
		  DIAG(i, j, c, p+vv[i][j])
		  /*fprintf(stdout, "I %d J %d c %d\n", i,j,c);*/
		  /*DIAG(i, j, c, p+va[B[j]])*/		/* diagonal */
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > lmin ){	/* add the score into list */
		    /*fprintf(stdout, "c: %d; lmin before: %d;  ", c, lmin);*/
		    /*fprintf(stdout, "K: %d", K);*/
		    lmin = addnode(c, ci, cj, i, j, K, lmin);
		    /*fprintf(stdout, "lmin after: %d\n", lmin); */
		  }
	        }
	  }
}

/* Determine the left and top boundaries of the recomputed area */

locate(nseq) int nseq;
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */ 
  register  int  di, dj;		/* end-point associated with d */
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  int  cflag, rflag;			/* for recomputation */
  int  *va;				/* pointer to vv(A[i], B[j]) */
  int   addnode();			/* function for inserting a node */
  int  limit;				/* the bound on j */

	/* Reverse pass
	   rows
	   CC : the scores on the current row
	   RR and EE : the endpoints that lead to CC
	   DD : the deletion scores 
	   SS and FF : the endpoints that lead to DD

	   columns
	   HH : the scores on the current columns
	   II and JJ : the endpoints that lead to HH
	   WW : the deletion scores
	   XX and YY : the endpoints that lead to WW
	*/
	for ( j = nn; j >= n1 ; j-- )
          {  CC[j] = 0;
	     EE[j] = j;
	     DD[j] = - (q);
	     FF[j] = j;
	     if ( nseq == 2 || j > mm )
                RR[j] = SS[j] = mm + 1;
	     else
                RR[j] = SS[j] = j;
	  }

        for ( i = mm; i >= m1; i-- )
	  {  c = p = 0;
	     f = - (q);
	     ci = fi = i;
	     pi = i + 1;
	     cj = fj = pj = nn + 1;
	     /*va = vv[A[i]];*/
	     if ( nseq == 2 || n1 > i )
		limit = n1;
	     else
		limit = i + 1;
	     for ( j = nn; j >= limit ; j-- )  
	       {  f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
		  /* DIAG(i, j, c, p+va[B[j]]) */
		  DIAG(i, j, c, p+vv[i][j])		/* diagonal */
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > lmin )
		    flag = 1;
	        }
	     if ( nseq == 2 || i < n1 )
	       { HH[i] = CC[n1];
	         II[i] = RR[n1];
	         JJ[i] = EE[n1];
	         WW[i] = DD[n1];
	         XX[i] = SS[n1];
	         YY[i] = FF[n1];
	       }
	  }
      
  for ( rl = m1, cl = n1; ; )
    { for ( rflag = cflag = 1; ( rflag && m1 > 1 ) || ( cflag && n1 > 1 ) ;  )
        { if ( rflag && m1 > 1 )	/* Compute one row */
            { rflag = 0;
	      m1--;
      	      c = p = 0;
	      f = - (q);
	      ci = fi = m1;
	      pi = m1 + 1;
	      cj = fj = pj = nn + 1;
	      /*va = vv[A[m1]]; */
	      for ( j = nn; j >= n1 ; j-- )  
	        { f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
		  /*DIAG(m1, j, c, p+va[B[j]])*/		/* diagonal */
		  DIAG(m1, j, c, p+vv[m1][j]);
		  if ( c <= 0 )
		    { c = 0; ci = m1; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > lmin )
		     flag = 1;
		  if ( ! rflag && ( ci > rl && cj > cl || di > rl && dj > cl
	 		                            || fi > rl && fj > cl ) )
		      rflag = 1;
	        }
	      HH[m1] = CC[n1];
	      II[m1] = RR[n1];
	      JJ[m1] = EE[n1];
	      WW[m1] = DD[n1];
	      XX[m1] = SS[n1];
	      YY[m1] = FF[n1];
	      if ( ! cflag && ( ci > rl && cj > cl || di > rl && dj > cl
			     || fi > rl && fj > cl ) )
	         cflag = 1;
	    }

	  if ( nseq == 1 && n1 == (m1 + 1) && ! rflag )
	     cflag = 0;
	  if ( cflag && n1 > 1 )	/* Compute one column */
	    { cflag = 0;
	      n1--;
	      c = 0;
	      f = - (q);
	      cj = fj = n1;
	      /*va = vv[B[n1]]; */
	      if ( nseq == 2 || mm < n1 )
		{ p = 0;
	          ci = fi = pi = mm + 1;
	          pj = n1 + 1;
		  limit = mm;
		}
	      else
		{ p = HH[n1];
		  pi = II[n1];
		  pj = JJ[n1];
	          ci = fi = n1;
		  limit = n1 - 1;
		}
	      for ( i = limit; i >= m1 ; i-- )  
	        { f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = HH[i] - qr; 
		  ci = II[i];
		  cj = JJ[i];
		  d = WW[i] - r;
		  di = XX[i];
		  dj = YY[i];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
	          /*DIAG(i, n1, c, p+va[A[i]]) */
		  DIAG(i, n1, c, p+vv[i][n1])
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = n1; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = HH[i];
		  HH[i] = c;
		  pi = II[i];
		  pj = JJ[i];
		  II[i] = ci;
		  JJ[i] = cj;
		  WW[i] = d;
		  XX[i] = di;
		  YY[i] = dj;
		  if ( c > lmin )
		     flag = 1;
	          if ( ! cflag && ( ci > rl && cj > cl || di > rl && dj > cl
		               || fi > rl && fj > cl ) )
		     cflag = 1;
	        }
	      CC[n1] = HH[m1];
	      RR[n1] = II[m1];
	      EE[n1] = JJ[m1];
	      DD[n1] = WW[m1];
	      SS[n1] = XX[m1];
	      FF[n1] = YY[m1];
	      if ( ! rflag && ( ci > rl && cj > cl || di > rl && dj > cl
		                                 || fi > rl && fj > cl ) )
	         rflag = 1;
	    }
	}
      if ( m1 == 1 && n1 == 1 || no_cross() )
	 break;
   }
  m1--;
  n1--;
}

/* recompute the area on forward pass */
small_pass(count,nseq) int count, nseq;
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */ 
  register  int  di, dj;		/* end-point associated with d */
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  int  *va;				/* pointer to vv(A[i], B[j]) */
  int   addnode();			/* function for inserting a node */
  int  limit;				/* lower bound on j */

	for ( j = n1 + 1; j <= nn ; j++ )
	  {  CC[j] = 0;
	     RR[j] = m1;
	     EE[j] = j;
	     DD[j] = - (q);
	     SS[j] = m1;
	     FF[j] = j;
	  }
	for ( i = m1 + 1; i <= mm; i++) 
	  {  c = 0;				/* Initialize column 0 */
	     f = - (q);
	     ci = fi = i;
	     /*va = vv[A[i]];*/
	     if ( nseq == 2 || i <= n1 )
	       { p = 0;
	         pi = i - 1;
	         cj = fj = pj = n1;
		 limit = n1 + 1;
	       }
	     else
	       { p = CC[i];
		 pi = RR[i];
		 pj = EE[i];
	         cj = fj = i;
		 limit = i + 1;
	       }
	     for ( j = limit ; j <= nn ; j++ )  
	       {  f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
		  /*DIAG(i, j, c, p+va[B[j]])*/		/* diagonal */
		  DIAG(i, j, c, p+vv[i][j])
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > lmin )	/* add the score into list */
		    lmin = addnode(c, ci, cj, i, j, count, lmin);
	        }
	  }
}

/* Add a new node into list.  */

int addnode(c, ci, cj, i, j, K, cost)  int c, ci, cj, i, j, K, cost;
{ int found;				/* 1 if the node is in LIST */
  register int d;

  found = 0;
  if ( most != 0 && most->STARI == ci && most->STARJ == cj )
    found = 1;
  else
     for ( d = 0; d < numnode ; d++ )
	{ most = LIST[d];
	  if ( most->STARI == ci && most->STARJ == cj )
	    { found = 1;
	      break;
	    }
        }
  if ( found )
    { if ( most->SCORE < c )
        { most->SCORE = c;
          most->ENDI = i;
          most->ENDJ = j;
        }
      if ( most->TOP > i ) most->TOP = i;
      if ( most->BOT < i ) most->BOT = i;
      if ( most->LLEFT > j ) most->LLEFT = j;
      if ( most->LRIGHT < j ) most->LRIGHT = j;
    }
  else
    { if ( numnode == K )	/* list full */
	 most = low;
      else
         most = LIST[numnode++];
      most->SCORE = c;
      most->STARI = ci;
      most->STARJ = cj;
      most->ENDI = i;
      most->ENDJ = j;
      most->TOP = most->BOT = i;
      most->LLEFT = most->LRIGHT = j;
    }
  if ( numnode == K )
    { if ( low == most || ! low ) 
        { for ( low = LIST[0], d = 1; d < numnode ; d++ )
            if ( LIST[d]->SCORE < low->SCORE )
              low = LIST[d];
	}
      return ( low->SCORE ) ;
    }
  else
    return cost;
}

/* Find and remove the largest score in list */

vertexptr findmax()
{ vertexptr  cur;
  register int i, j;

  for ( j = 0, i = 1; i < numnode ; i++ )
    if ( LIST[i]->SCORE > LIST[j]->SCORE )
       j = i;
  cur = LIST[j];
  if ( j != --numnode )
    { LIST[j] = LIST[numnode];
      LIST[numnode] =  cur;
    }
  most = LIST[0];
  if ( low == cur ) low = LIST[0];
  return ( cur );
}

/* return 1 if no node in LIST share vertices with the area */

no_cross()
{ vertexptr  cur;
  register int i;

      for ( i = 0; i < numnode; i++ )
	{ cur = LIST[i];
	  if ( cur->STARI <= mm && cur->STARJ <= nn && cur->BOT >= m1-1 && 
	       cur->LRIGHT >= n1-1 && ( cur->STARI < rl || cur->STARJ < cl ))
	     { if ( cur->STARI < rl ) rl = cur->STARI;
	       if ( cur->STARJ < cl ) cl = cur->STARJ;
	       flag = 1;
	       break;
	     }
	}
      if ( i == numnode )
	return 1;
      else
	return 0;
}

/* diff(A,B,M,N,tb,te) returns the score of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */

int diff(A,B,M,N,tb,te) int A, B, M, N; int tb, te;

{ int   midi, midj, type;	/* Midpoint, type, and cost */
  int midc;

{ register int   i, j;
  register int c, e, d, s;
           int t, *va;
#ifdef FAR_PTR
	char far * flckalloc();
#else
	char  *lckalloc();
#endif

/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0)
    { if (M > 0) DEL(M)
      return - gap(M);
    }
  if (M <= 1)
    { if (M <= 0)
        { INS(N);
          return - gap(N);
        }
      if (tb > te) tb = te;
      midc = - (tb + r + gap(N) );
      midj = 0;
      /*va = vv[A[1]];*/
      for (j = 1; j <= N; j++)
        {  for ( tt = 1, z = row[I+1]; z != PAIRNULL; z = z->NEXT )	
              if ( z->COL == j+J )			
	         { tt = 0; break; }		
           if ( tt )			
            { /*c = va[B[j]] - ( gap(j-1) + gap(N-j) );*/
		c= vv[A+1][B+j] - ( gap(j-1) + gap(N-j) );
		/*fprintf(stdout, "1 j: 1 %d; score: %d\n", j, vv[1][j]); */
              if (c > midc)
               { midc = c;
                 midj = j;
               }
	    }
	}
      if (midj == 0)
        { INS(N) DEL(1) }
      else
        { if (midj > 1) INS(midj-1)
          REP
	  /*if ( A[1] == B[midj] )
	     no_mat += 1;
	  else
	     no_mis += 1;
	  */
	  /* mark (A[I],B[J]) as used: put J into list row[I] */	
          I++; J++;
	  z = ( pairptr ) FCKALLOC( (int) sizeof(pair));
          z->COL = J;			
          z->NEXT = row[I];				
	  row[I] = z;
          if (midj < N) INS(N-midj)
        }
      return midc;
    }

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;			/* Forward phase:                          */
  CC[0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
  t = -q;
  for (j = 1; j <= N; j++)
    { CC[j] = t = t-r;
      DD[j] = t-q;
    }
  t = -tb;
  for (i = 1; i <= midi; i++)
    { s = CC[0];
      CC[0] = c = t = t-r;
      e = t-q;
      /*va = vv[A[i]];*/
      for (j = 1; j <= N; j++)
        { if ((c = c - qr) > (e = e - r)) e = c;
          if ((c = CC[j] - qr) > (d = DD[j] - r)) d = c;
	  /*DIAG(i+I, j+J, c, s+va[B[j]]) */
	  DIAG(i+I, j+J, c, s+vv[A+i][B+j])
	  /*fprintf(stdout, "i j: %d %d; score: %d\n", i, j, vv[A+i][B+j]); */
          if (c < d) c = d;
          if (c < e) c = e;
          s = CC[j];
          CC[j] = c;
          DD[j] = d;
        }
    }
  DD[0] = CC[0];

  RR[N] = 0;			/* Reverse phase:                          */
  t = -q;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
  for (j = N-1; j >= 0; j--)
    { RR[j] = t = t-r;
      SS[j] = t-q;
    }
  t = -te;
  for (i = M-1; i >= midi; i--)
    { s = RR[N];
      RR[N] = c = t = t-r;
      e = t-q;
      /*va = vv[A[i+1]];*/
      for (j = N-1; j >= 0; j--)
        { if ((c = c - qr) > (e = e - r)) e = c;
          if ((c = RR[j] - qr) > (d = SS[j] - r)) d = c;
	  /*DIAG(i+1+I, j+1+J, c, s+va[B[j+1]]) */
	  DIAG(i+1+I, j+1+J, c, s+vv[A+i+1][B+j+1])
	  /*fprintf(stdout, "i+1, j+1: %d %d; score: %d\n", i+1, j+1, vv[i+1][j+1]);*/
          if (c < d) c = d;
          if (c < e) c = e;
          s = RR[j];
          RR[j] = c;
          SS[j] = d;
        }
    }
  SS[N] = RR[N];

  midc = CC[0]+RR[0];		/* Find optimal midpoint */
  midj = 0;
  type = 1;
  for (j = 0; j <= N; j++)
    if ((c = CC[j] + RR[j]) >= midc)
      if (c > midc || CC[j] != DD[j] && RR[j] == SS[j])
        { midc = c;
          midj = j;
        }
  for (j = N; j >= 0; j--)
    if ((c = DD[j] + SS[j] + q) > midc)
      { midc = c;
        midj = j;
        type = 2;
      }
}

/* Conquer: recursively around midpoint */

  if (type == 1)
    { diff(A,B,midi,midj,tb,q);
      diff(A+midi,B+midj,M-midi,N-midj,q,te);
    }
  else
    { diff(A,B,midi-1,midj,tb,zero);
      DEL(2);
      diff(A+midi+1,B+midj,M-midi-1,N-midj,zero,te);
    }
  return midc;
}

/* Alignment display routine */
/*
static char ALINE[51], BLINE[51], CLINE[51];

int display(A,B,M,N,S,AP,BP) char A[], B[]; int M, N; int S[], AP, BP;
{ register char *a, *b, *c;
  register int   i,  j, op;
           int   lines, ap, bp;

  i = j = op = lines = 0;
  ap = AP;
  bp = BP;
  a = ALINE;
  b = BLINE;
  c = CLINE;
  while (i < M || j < N)
    { if (op == 0 && *S == 0)
        { op = *S++;
          *a = sq[A[++i]];
          *b = sq[B[++j]];
          *c++ = (*a++ == *b++) ? '|' : ' ';
        }
      else
        { if (op == 0)
            op = *S++;
          if (op > 0)
            { *a++ = ' ';
              *b++ = sq[B[++j]];
              op--;
            }
          else
            { *a++ = sq[A[++i]];
              *b++ = ' ';
              op++;
            }
          *c++ = '-';
        }
      if (a >= ALINE+50 || i >= M && j >= N)
        { *a = *b = *c = '\0';
          printf("\n%5d ",50*lines++);
          for (b = ALINE+10; b <= a; b += 10)
            printf("    .    :");
          if (b <= a+5)
            printf("    .");
          printf("\n%5d %s\n      %s\n%5d %s\n",ap,ALINE,CLINE,bp,BLINE);
	  ap = AP + i;
	  bp = BP + j;
          a = ALINE;
          b = BLINE;
          c = CLINE;
        }
    }
}
*/

/* lib.c - library of C procedures. */

/* lfatal - print message and die */
lfatal(msg)
char *msg;
{
	fprintf(stderr, "%s\n", msg);
	exit(1);
}

/* lfatalf - format message, print it, and die */
lfatalf(msg, val)
char *msg, *val;
{
	fprintf(stderr, msg, val);
	putc('\n', stderr);
	exit(1);
}
	
/* ckopen - open file; check for success */
FILE *ckopen(name, mode)
char *name, *mode;
{
	FILE *fopen(), *fp;

	if ((fp = fopen(name, mode)) == NULL)
		lfatalf("Cannot open %s.", name);
	return(fp);
}

/* lckalloc - allocate space; check for success */
char *lckalloc(amount)
int amount;
{
	char *malloc(), *p;
	static long mtotal;

	mtotal += (long)amount;

	if ((p = malloc( (unsigned) amount)) == NULL) {
	fprintf(stderr,"Ran out of near memory: %d/%ld\n",amount,mtotal);
		exit(1);
	}
	return(p);
}

#ifdef FAR_PTR
#define FMALLOC farmalloc
#define MTYPE long
#define FFREE farfree

/* flckalloc - allocate space; check for success */
char far *flckalloc(amount)
	int amount;
{
	static long ftotal;
	static int nf;

	char far * FMALLOC(), far * p;

	ftotal += (long)amount;
	nf++;

	if ((p = FMALLOC( (MTYPE) amount)) == (char far *)NULL) {
	fprintf(stderr,"Ran out of far memory: %d/%ld (%d)\n",
		amount,ftotal,nf);
		exit(1);
	}
	return(p);
}
#endif
