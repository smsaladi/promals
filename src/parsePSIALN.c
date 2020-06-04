#include "header_cpp.h"
#include "amino.h"
#include "util.h"
#include "subalign.h"

int debug  = 1;
double low_id_thr = 0.25; // low threshold for removing sequences
double high_id_thr = 0.98; // high threshold for removing redundancies
int N = -1; // number of sequences selected
int S = 1; // shrink ali or not
int nonumber = 1;
double len_fraction = 0.5;

void purge_ali_by_mark(int n, int len, char **aseq, double low_id_thr, double high_id_thr, double len_fraction, int *marks);
void shrinkali(int n, int len, char **name, char **seq, int *start, char **seq1, int *len1);
void printali_marks(int chunk, int n, int len, char **name, char **seq, int *start, int *marks);
void printali_marks_nonumber(int chunk, int n, int len, char **name, char **seq, int *start, int *marks);

int  main(int argc, char **argv)
{
        int i, n;
        int *marks;
        int pn; // the number of printed sequences

        if (argc < 2) {
                fprintf(stderr, "Usage: readali file\n");
                fprintf(stderr, "-l: low threshold [0..1]\n");
                fprintf(stderr, "-h: high threhold [0..1]\n");
                fprintf(stderr, "-n: sequence number \n");
                fprintf(stderr, "-s: shrink alignment according the first sequence\n");
		fprintf(stderr, "-g: length fraction [0..1]\n");
                exit(1);
        }

        for(i=1;i<argc;i++) {
           if(strcmp(argv[i], "-l")==0) low_id_thr = atof(argv[i+1]);
           if(strcmp(argv[i], "-h")==0) high_id_thr = atof(argv[i+1]);
           if(strcmp(argv[i], "-n")==0) N = atoi(argv[i+1]);
           if(strcmp(argv[i], "-s")==0) S = atoi(argv[i+1]);
           if(strcmp(argv[i], "-no")==0) nonumber = atoi(argv[i+1]);
	   if(strcmp(argv[i], "-g")==0) len_fraction = atof(argv[i+1]);

           if( (low_id_thr>1.0) || (low_id_thr < 0) ) { fprintf(stdout, "low threshold must between 0 and 1\n"); exit(0); }
           if( (high_id_thr>1.0) || (high_id_thr < 0) ) { fprintf(stdout, "high threshold must between 0 and 1\n"); exit(0); }
        }

        if(debug>11)fprintf(stdout, "%f %f\n", low_id_thr, high_id_thr);

	subalign *a = new subalign(argv[1]);

        marks = ivector(a->nal); 
	

	purge_ali_by_mark(a->nal, a->alilen, a->aseq, low_id_thr, high_id_thr, len_fraction, marks);

        // selected sequence number
        //fprintf(stdout, "N: %d\n", N);
        if(N>0) {
            pn = 0;
            for(i=0;i<a->nal;i++) {
                if(marks[i]) { pn++;}
            }
            //fprintf(stdout, "pn: %d\n", pn);
            if(pn>N) {
                for(i=a->nal-1;i>0;i--) {
                    if(marks[i]) {marks[i]=0;pn--;}
                    if(pn==N) break;
                }
            }
        }

	int alilen1;
	char **seq1;
        if(S) {
            seq1 = (char **) malloc( (a->nal+1) * sizeof (char *) );
            for(i=0;i<=a->nal;i++) seq1[i] = (char *) malloc( (a->alilen+1) * sizeof( char ) );
            shrinkali(a->nal, a->alilen, a->aname, a->aseq, a->astart, seq1, &alilen1);
            if(!nonumber) printali_marks(80, a->nal, alilen1, a->aname, seq1, a->astart, marks);
            else printali_marks_nonumber(80, a->nal, alilen1, a->aname, seq1, a->astart, marks);
        }
        else {
           if(!nonumber) printali_marks(80, a->nal, alilen1, a->aname, seq1, a->astart, marks);
           else printali_marks_nonumber(80, a->nal, a->alilen, a->aname, a->aseq, a->astart, marks);
        }

        return (0);



}



#define MIN(x,y) ( (x) < (y) ? (x): (y) )

// for an alignment with sequence number n and alignment length len and amino acid sequences aseq,
// mark those sequences that have sequences with mutual sequence similarity below a high threshold (
// high_id_thr) and shows sequence identity to the first sequence above a low threshold (low_id_thr)
void purge_ali_by_mark(int n, int len, char **aseq, double low_id_thr, double high_id_thr, double len_fraction, int *marks)
{

	int i,j,k;
	int *aalens, count;
	int *hms;
	int querylength=0;

	// marks for the sequences
	//marks = (int *) malloc ( (n) * sizeof (int ) );

	// amino acid lenghs according to hms
	aalens = ivector(n);
	aalens = (int *) malloc ( (n) * sizeof( int ) );

	// positions with the first sequence occupied by amino acids
	// horizontal marks
	hms = ivector(len);
	hms = (int *) malloc( (len) * sizeof( int ) );

  	for(i=0;i<n;i++) {marks[i] = 1; aalens[i] = 0; }
	for(i=0;i<len;i++) hms[i] = 0;

	// mark the positions with hms
 	for(i=0;i<len;i++) {
	    if(aseq[0][i] == 'X') {aseq[0][i] = 'A'; }
	    if(aseq[0][i] != '-') {
		hms[i] = 1;
		querylength++;
	    }
	    if(debug>11) fprintf(stdout, "%d", hms[i]);
	}
	if(debug>11) fprintf(stdout, "\n");
	if(debug>11) fprintf(stdout, "querylength: %d\n", querylength);

	// aalens: 
	for(i=0;i<n;i++) {
	    for(j=0;j<len;j++) {
		if(hms[j])
	    	if(aseq[i][j] != '-') {
		    aalens[i]++;
		}
	    }
	    if(debug>11) fprintf(stdout, "aalens: %d %d \n", i, aalens[i]);
	}

	// mark by low identity threshold
	for(i=1;i<n;i++) {
	   count = 0;
	   for(k=0;k<len;k++) {
	      if(hms[k]) {
		if(aseq[0][k]==aseq[i][k]) count++;
	      }
	   }
	   if( count < low_id_thr * MIN(aalens[0], aalens[i]) ) {
		marks[i] = 0;
	   }
	}

	// mark by high identity threshold
	for(i=0;i<n-1;i++) {
	   if(marks[i]==0) continue;
	   for(j=i+1;j<n;j++) {
		if(marks[j]==0) continue;
		count = 0;
		for(k=0;k<len;k++) {
		    if(hms[k]) {
			if( (aseq[i][k]=='-') &&( aseq[j][k]=='-') ) continue;
			if(aseq[i][k]==aseq[j][k]) {
			   count++;
			}
		    }
		}
		if( (count / high_id_thr) > MIN (aalens[i], aalens[j]) ) {
		    marks[j] = 0;
		}
	    }
	}

	// mark by length fraction
	for(i=1;i<n;i++) {
		if(marks[i]==0) continue;
		if (aalens[i]*1.0/aalens[0] < len_fraction) marks[i] = 0;
	}

	delete [] aalens; delete [] hms;
		
}
	

// shrinkali according to the gap pattern of the first sequence
// if the position is a gap in the first sequence, it is deleted
void shrinkali(int n, int len, char **name, char **seq, int *start, char **seq1, int *len1)
{
        int i, j,k;
        int jj, l;
        char *sq;

        jj = 0;
        for(i=0;i<len;i++) {
             if(seq[0][i] != '-') {
                for(j=0;j<n;j++) {
                    seq1[j][jj] = seq[j][i];
                }
                jj++;
             }
        }

        *len1 = jj;
}

// print alignment consisting of sequences that are marked
void printali_marks(int chunk, int n, int len, char **name, char **seq, int *start, int *marks)
{
        int i, j, jj, mlen, sta, *pos;
        char *sq;

        pos = (int *) malloc(n * sizeof(pos[0]));
        memcpy(pos, start, n * sizeof(pos[0]));
        for (i=0; i < n && start[i] == 0; i++) ;
        sta = (i < n);
        for (i=1, mlen=strlen(name[0]); i < n; i++) {
                if (mlen < strlen(name[i])) {
                        mlen = strlen(name[i]);
                }
        }
        jj = 0;
        do {
                if (jj > 0) {
                        printf("\n");
                }
                for (i=0; i < n; i++) {
                    if(marks[i]) {
                        printf("%-*s ", mlen, name[i]);
                        if (sta) {
                                printf("%4d ", pos[i]);
                        }
                        sq = seq[i] + jj;
                        for (j=0; j+jj < len && j < chunk; j++) {
                                if (isalpha(sq[j])) {
                                        pos[i]++;
                                }
                                printf("%c", sq[j]);
                        }
                        if (sta) {
                                printf(" %4d", pos[i]-1);
                        }
                        printf("\n");
                    }
                }
                jj += j;
        } while (jj < len);
        free(pos);
}




void printali_marks_nonumber(int chunk, int n, int len, char **name, char **seq, int *start, int *marks)
{
        int i, j, jj, mlen, sta, *pos;
        char *sq;

        pos = (int *) malloc(n * sizeof(pos[0]));
        memcpy(pos, start, n * sizeof(pos[0]));
        for (i=0; i < n && start[i] == 0; i++) ;
        sta = (i < n);
        for (i=1, mlen=strlen(name[0]); i < n; i++) {
                if (mlen < strlen(name[i])) {
                        mlen = strlen(name[i]);
                }
        }
	cout << endl;
        jj = 0;
        do {
                if (jj > 0) {
                        printf("\n");
                }
                for (i=0; i < n; i++) {
                    if(marks[i]) {
                        printf("%-*s ", mlen, name[i]);
                        /* if (sta) {
                                printf("%4d ", pos[i]);
                        } */
                        sq = seq[i] + jj;
                        for (j=0; j+jj < len && j < chunk; j++) {
                                if (isalpha(sq[j])) {
                                        pos[i]++;
                                }
                                printf("%c", sq[j]);
                        }
                        /* if (sta) {
                                printf(" %4d", pos[i]-1);
                        } */
                        printf("\n");
                    }
                }
                jj += j;
        } while (jj < len);
	cout << endl;
        free(pos);
}


