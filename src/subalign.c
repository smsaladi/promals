#include "header_cpp.h"
#include "subalign.h"
#include "util.h"
#include "amino.h"

static int mydebug = 1;

// Initilize static data member at file scope
int subalign::subaligncount=0;

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
	beta = 10;

	subaligncount++;
}

void subalign::set_gapt(double gapt) {

	gapt_threshold = gapt;

}

// Copy constructor
subalign::subalign(const subalign &init) {
	int i, j;

	nal = init.nal;
	alilen = init.alilen;
	mnamelen = init.mnamelen;
	beta = init.beta;

	//cout << "nal: " << nal << " alilen: " << alilen << " mnamelen " << mnamelen << endl; 
	//cout << init.aname[0] << endl;
	//cout << init.aseq[0] << endl;
	//cout << init.nal << endl;
	//cout << init.alilen << endl;
	//cout << init.mnamelen << endl;
	
	alignment = imatrix(nal, alilen);
	aname = cmatrix(nal, mnamelen+1);	
	aseq = cmatrix(nal, alilen+1);
	
	for(i=1;i<=nal;i++) {
	   for(j=1;j<=alilen;j++) {
	        //cout << i << endl;
		alignment[i][j] = init.alignment[i][j];
		aseq[i-1][j-1] = init.aseq[i-1][j-1];
		//cout << alignment[i][j];
	   }
	   aseq[i-1][alilen] = '\0';
	   //cout << endl;
	}
	for(i=0;i<nal;i++) {
	   strcpy(aname[i], init.aname[i]);
	}

	
	subaligncount++;
	//cout << "alignment count: " << subaligncount << endl;

}

// Constructor reading an input alignment file
subalign::subalign(char *filename) {

	int i,j,k;

	readali(filename);

	// new
	reassign();

	alignment = imatrix(nal, alilen);

	for(i=1;i<=nal;i++) {
	   for(j=1;j<=alilen;j++) {
		alignment[i][j] = am2num(aseq[i-1][j-1]);
	   }
	}

	gap_threshold = 0.5;
	gapt_threshold = 1.0;
	done_profile = 0;
	done_prof = 0;
	beta = 10;

	subaligncount++;
}

// Constructor reading an input file name
subalign::subalign(string filename) {

        int i,j,k;

	char filename1[100];
	for(i=0;i<=filename.length();i++) {
	   filename1[i] = filename[i];
	}
        readali(filename1);

        // new
        reassign();

        alignment = imatrix(nal, alilen);

        for(i=1;i<=nal;i++) {
           for(j=1;j<=alilen;j++) {
                alignment[i][j] = am2num(aseq[i-1][j-1]);
           }
        }

        gap_threshold = 0.5;
        gapt_threshold = 1.0;
        done_profile = 0;
        done_prof = 0;
	beta = 10;

	subaligncount++;
}

// Destructor
subalign::~subalign() {
	int i;

	if(alignment) { for(i=0;i<=nal;i++) delete [] alignment[i]; delete [] alignment; }
	if(aseq) { for(i=0;i<=nal;i++) delete [] aseq[i]; delete [] aseq; }
	if(aname) { for(i=0;i<=nal;i++) delete [] aname[i]; delete [] aname; }
	subaligncount--;
	//cout << "alignment count decrease: " << subaligncount << endl;

	if(done_profile) {
		delete [] apos_filtr;
		for(i=0;i<=alilen;i++) delete [] n_effAa[i]; delete [] n_effAa;
		for(i=0;i<=alilen;i++) delete [] pseudoCnt[i]; delete [] pseudoCnt;
		delete [] sum_eff_let;
		delete [] maskgapRegion;
	}
	if(done_prof) {
		delete [] prof_pos;
		delete [] prof_sum_eff;
		delete [] prof_hwt_all;
		delete [] prof_gap_content;
		free_dmatrix(prof_effn, prof_len, 20);
		free_dmatrix(prof_freq, prof_len, 20);
	}
}

// Overload assignment operator
const subalign & subalign::operator=( const subalign &right) {

	int i,j,k;
	
	if(&right !=this) {
	    // delete the original arrays
	    if(alignment) {
		for(i=0;i<=nal;i++) {
		   delete [] alignment[i];
		}
		delete [] alignment;
	    }
	    if(aseq) {
		for(i=0;i<=nal;i++) {
		   delete [] aseq[i];
		}
		delete [] aseq;
	    }
	    if(aname) {
		for(i=0;i<=nal;i++) {
		   delete [] aname[i];
		}
	        delete [] aname;
	    }
	
	    // assignment of new values
	    nal = right.nal;
	    alilen = right.alilen;
	    mnamelen = right.mnamelen;
	    alignment = imatrix(nal, alilen);
	    aseq = cmatrix(nal, alilen);
	    aname = cmatrix(nal, mnamelen);
	    for(i=1;i<=nal;i++) {
		for(j=1;j<=alilen;j++) {
		    alignment[i][j] = right.alignment[i][j];
		    aseq[i-1][j-1] = right.aseq[i-1][j-1];
		}
		aseq[i-1][alilen] = '\0';
		strcpy(aname[i-1], right.aname[i-1]);
	    }
	    beta = right.beta;
	}
	return *this;
}

	    
// Return the number of subalign objects intantiated
int subalign::getsubaligncount() {
	return subaligncount;
}


/******************************************************************/
/*adapted from align.c, auxilary routines used for readali...*/
/*memory allocations follow the C language constoms */

void *subalign::mymalloc(int size)
{
        void *buf;

        if ((buf = malloc(size)) == NULL) {
                fprintf(stderr, "Not enough memory: %d\n", size);
                exit(1);
        }
        return buf;
}

char * subalign::strsave(char *str)
{
        char *buf;
        int l;

        l = strlen(str);
        buf = (char *) mymalloc((l + 1)*sizeof(str[0]));
        strcpy(buf, str);
        return buf;
}

char * subalign::strnsave(char *str, int l)
{
        char *buf;

        buf = (char *)mymalloc(l + 1);
        memcpy(buf, str, l);
        buf[l] = '\0';
        return buf;
}

char ** subalign::incbuf(int n, char **was)
{
        char **buf;

        buf = (char **)mymalloc((n+1) * sizeof(buf[0]));
        if (n > 0) {
                memcpy(buf, was, n * sizeof(was[0]));
                free(was);
        }
        buf[n] = NULL;
        return buf;
}

int * subalign::incibuf(int n, int *was)
{
        int *ibuf;

        ibuf = (int *) mymalloc((n+1) * sizeof(ibuf[0]));
        if (n > 0) {
                memcpy(ibuf, was, n * sizeof(was[0]));
                free(was);
        }
        ibuf[n] = 0;
        return ibuf;
}

/* end of auxilary routines from align.c ****************************/
/********************************************************************/

        
int subalign::getNal() 
{
        return (nal);
}

int subalign::getAlilen()
{
        return (alilen);
}

int **subalign::getAlignment() 
{
        /* int **ali;       
        int i, j;
        
        ali = imatrix(0, nal+1, 0, alilen+1);
        for(j=1;j<=nal;j++)
        for(i=1;i<=alilen;i++) {
                ali[j][i] = alignment[j][i];
        } */

        return alignment;
}

char **subalign::getAseq()
{
        /* char **seq;
        int i, j;
        
        seq = cmatrix(0, nal+1, 0, alilen+1);
        for(j=0;j<nal;j++)
        for(i=0;i<alilen;i++) {
                seq[j][i] = aseq[j][i];
        } */

        return aseq;
}

char **subalign::getAname()
{
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
        
        return(aname);
}

        

void subalign::setNal(int n)
{
        nal = n;
}

void subalign::setAlilen(int len)
{
        alilen = len;
}

void subalign::setAlignment(int **ali)
{
        int i, j;

        if(ali==NULL) {
                fprintf(stderr, "WARNING:: alignment is null\n");
                return;
        }

	alignment = new int * [nal+1];
	for(i=1;i<=nal;i++) alignment[i] = new int [alilen+1];
        //fprintf(stderr, "nal: %d\n", nal);
        for(i=1;i<=nal;i++){
        for(j=1;j<=alilen;j++){
                alignment[i][j] = ali[i][j];
//              fprintf(stderr,"%d",ali[i][j]);
        }
        }
        
        /*alignment = ali;*/

}

void subalign::setAseq(char **seq)
{
        int i, j;

        if(seq==NULL) {
                fprintf(stderr, "WARNING:: alignment sequences are null\n");
                return;
        }

	aseq = new char * [nal];
	for(i=0;i<nal;i++) aseq[i] = new char [alilen+1];
        for(i=0;i<nal;i++) {
	     for(j=0;j<alilen;j++) aseq[i][j] = seq[i][j];
	     aseq[i][alilen] = '\0';
	}

}

void subalign::setAname(char **name)
{
        int i, j;
        int mlen;
        
        if(name == NULL) {
                fprintf(stderr, "ALign name empty\n");
                return;
        }

        //fprintf(stderr, "nal: %d\n", nal);
        for (i=0, mlen=strlen(name[0]); i < nal; i++) {
                if (mlen < strlen(name[i])) {
                        mlen = strlen(name[i]);
                }
        }

	aname = new char * [ nal];
	for(i=0;i<nal;i++) aname[i] = new char [mlen+1];
        for(i=0;i<nal;i++) strcpy(aname[i], name[i]);
	mnamelen = mlen;
        
}

// Read an input alignment in ClustalW (with or without header) format
void subalign::readali(char *filename)
{
        FILE *fp;
        char *s, *ss, *seqbuf;
        int n, l, len, len0;
        char str[10001];
	int MAXSTR = 100001;

        if ((fp = fopen(filename, "r")) == NULL) {
                fprintf(stderr, "No such file: \"%s\"\n", filename);
                exit(1);
        }
        alilen = 0;
        nal = 0;
        n = 0;
        while (fgets(str, MAXSTR, fp) != NULL) {
		if(strncmp(str, "CLUSTAL ", 8)==0) continue;
                for (ss = str; isspace(*ss); ss++) ;
		if( (*ss=='.') || (*ss=='*') || (*ss==':') ) continue;
                if (*ss == '\0') {
                        if (n == 0) {
                                continue;
                        }
                        if (nal == 0) {
                                if (n == 0) {
                                        fprintf(stderr, "No alignments read\n");
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
                for (s = ss; *s != '\0' && !isspace(*s); s++) ;
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
                for (ss = s; isspace(*ss); ss++);
                for (s = ss; isdigit(*s); s++) ;
                if (isspace(*s)) {
                        if (nal == 0) {
                                astart[n] = atoi(ss);
                                /*fprintf(stderr, "n: %d  astart: %d\n", n, astart[n]);*/
                        }
                        for (ss = s; isspace(*ss); ss++);
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
                        seqbuf = (char *) mymalloc(alilen+1);
                        memcpy(seqbuf, aseq[n], alilen-len);
                        free(aseq[n]);
                        aseq[n] = seqbuf;
                        memcpy(seqbuf+alilen-len, ss, len);
                        seqbuf[alilen] = '\0';
                }
                n++;
		//cout << n << endl;
        }
        if (nal == 0) {
                if (n == 0) {
                        fprintf(stderr, "No alignments read\n");
                        exit(1);
                }
                nal = n;
        } else if (n != 0 && n != nal) {
                fprintf(stderr, "Wrong nal, was: %d, now: %d\n", nal, n);
                exit(1);
        }
	free(astart); free(alen);
        fclose(fp);
}

// convert the C-style malloc allocations to C++ style new allocations
void subalign::reassign() {
	int i,j,k;
	int mlen;
	char **tmpaseq;
	char **tmpaname;
	
	tmpaseq = cmatrix(nal, alilen);

	for (i=1, mlen=strlen(aname[0]); i < nal; i++) {
                if (mlen < strlen(aname[i])) {
                        mlen = strlen(aname[i]);
                }
        }
	mnamelen = mlen;
	tmpaname = cmatrix(nal, mlen+1);

	for(i=0;i<nal;i++) {
	    strcpy(tmpaseq[i], aseq[i]);
	    strcpy(tmpaname[i], aname[i]);
	    free(aseq[i]); free(aname[i]);
	}
	free(aseq); free(aname);
	aseq = cmatrix(nal, alilen);
	aname = cmatrix(nal, mlen+1);
	for(i=0;i<=nal;i++) {
	    strcpy(aseq[i], tmpaseq[i]);
	    strcpy(aname[i], tmpaname[i]);
	    aseq[i][alilen] = '\0';
	    delete [] tmpaseq[i];
	    delete [] tmpaname[i];
	}
	delete [] tmpaseq; delete [] tmpaname;
}

// Print alignment with a given block size
void subalign::printali(int blocksize, int header) {

	int i,j,k;
	
	if(header) {
	cout << endl << endl;
	cout << "blocksize: " << blocksize << endl;
	cout << "alilen: " << alilen << endl;
	cout << "nal: " << nal << endl;
	}
	int nblocks = (alilen-1)/blocksize + 1;
	if(header) {
	cout << "nblocks: " << nblocks << endl;
	}
	cout << endl;

	for(i=1;i<=nblocks;i++) {
	     for(j=0;j<nal;j++) {
		cout << setiosflags(ios::left) << setw(mnamelen+3) << aname[j];
	        for(k=1;k<=blocksize;k++) {
		     if((i-1)*blocksize+k-1>=alilen) break;
		     cout << aseq[j][(i-1)*blocksize+k-1];
		}
		cout << endl;
	     }
	     cout << endl << endl;
	}
}

// Print alignment with a given block size
void subalign::printali(int blocksize) {

	int i,j,k;
	
	cout << endl << endl;
	cout << "blocksize: " << blocksize << endl;
	cout << "alilen: " << alilen << endl;
	cout << "nal: " << nal << endl;
	int nblocks = (alilen-1)/blocksize + 1;
	cout << "nblocks: " << nblocks << endl;
	cout << endl;

	for(i=1;i<=nblocks;i++) {
	     for(j=0;j<nal;j++) {
		cout << setiosflags(ios::left) << setw(mnamelen+3) << aname[j];
	        for(k=1;k<=blocksize;k++) {
		     if((i-1)*blocksize+k-1>=alilen) break;
		     cout << aseq[j][(i-1)*blocksize+k-1];
		}
		cout << endl;
	     }
	     cout << endl << endl;
	}
}

// Print alignment to a file
void subalign::printali(char *filename, int blocksize) {
	int i,j,k;
	ofstream outfile(filename, ios::out);
	if(!outfile) {
	    cout << "cannot write the alignment to "<< filename<< endl;
	    exit(1);
	}
	outfile << "CLUSTAL W format multiple sequence alignment"<<endl;
	outfile << endl << endl;
        int nblocks = (alilen-1)/blocksize + 1;

        for(i=1;i<=nblocks;i++) {
             for(j=0;j<nal;j++) {
                outfile << setiosflags(ios::left) << setw(mnamelen+3) << aname[j];
                for(k=1;k<=blocksize;k++) {
                     if((i-1)*blocksize+k-1>=alilen) break;
                     outfile << aseq[j][(i-1)*blocksize+k-1];
                }
                outfile << endl;
             }
             outfile << endl << endl;
        }

	//outfile << "<!---mummals_now_finished--->" << endl;
	outfile.close();
}
	
	
/* from a given alignment with aa as numbers, computes effective aa counts (PSIC->our formula)
and marks the columns with EFFECTIVE content of gaps > threshold (effgapmax) */
void subalign::neffsForEachCol_maskGapReg(int **ali, int n, int len, double effgapmax, double effgapRegionMin, double **n_effAa, double *sum_eff_let, int *maskgapRegion, int *apos_filtr, int *len_lowgaps, double *nef)
{
        int i,j,k,l;
        int alilen_mat, nsymbols_col, nsymbols;
        double nef_loc;
        int ele;
        double effnu[21];
        double sum_let;
        int *mark;
        int flagmark;

	gap_content = new double [alilen+1];


	mark = new int [n+11];
        alilen_mat = 0;
        nsymbols = 0;
        for(j=1;j<=len;j++) {

                nsymbols_col = 0;
                sum_let=0;

                for(k=0;k<=20;++k){
/* Mark sequences that have amino acid  k (or gap, k=0) in this jth position */
                         flagmark =0;
                        for(i=1;i<=n;++i){
                                mark[i]=0;

                                ele=ali[i][j];
                                if(ele==k){mark[i]=1; flagmark =1;}
                                ele=ali[i][j]-25;
                                if(ele==k) {mark[i]=1; flagmark =1;}
                        }

/* If aa k (or gap) is present in this position call compute k-th effective count */
                      if (flagmark == 1) {

                                effnu[k]=effective_number_nogaps(ali,mark,n,1,len);
                                nsymbols_col++;

                        } else { effnu[k] = 0.0; }

                       if (k>0) sum_let += effnu[k];
                }

		// calculate raw gap fraction
		if(mydebug>100) {
		int gap_count = 0;
		for(i=1;i<=n;i++) {
		    if(ali[i][j]==0) gap_count++;
		}
		double gapfract = 1.0 * gap_count / n;
		//cout << j << "  " << 1.0*effnu[0]/(sum_let + effnu[0]) << "  " << gapfract << endl;

		}

		if(mydebug>100) cout << j << "  " << aseq[0][j-1] << " " << sum_let << "  " << effnu[0] << endl;

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

                }


        }

	//cout << nsymbols << "\t" << alilen_mat << endl;

        nef_loc = 1.0*nsymbols/alilen_mat;
        *nef = nef_loc;
        *len_lowgaps = alilen_mat;

	average_sum_eff_let = 0;
	for(i=1;i<=alilen_mat;i++) {
		average_sum_eff_let += sum_eff_let[i];
	}
	average_sum_eff_let /= alilen_mat;

        maskgapRegion[0] = maskgapRegion[alilen_mat+1] = 0;

        delete [] mark;

	/* for(i=1;i<=alilen_mat;i++) {
	    for(j=1;j<=20;j++) {
		cout << n_effAa[i][j] << " ";
	    }
	    cout << endl;
	} */

}


double subalign::effective_number_nogaps(int **ali, int *marks, int n, int start, int end){

/* from the alignment of n sequences ali[1..n][1..l]
calculates effective number of sequences that are marked by 1 in mark[1..n]
for the segment of positions ali[][start..end]
Neff=ln(1-0.05*N-of-different-letters-per-site)/ln(0.95)
*/

int i,k,a,flag;
int amco[21],lettercount=0,sitecount=0;
double letpersite=0,neff;

for(k=start;k<=end;++k){
/******************DUMP the condition "consider only positions without gaps in the marked seqs" ***
********/
/*****  flag=0;for(i=1;i<=n;++i)if(marks[i]==1 && ali[i][k]==0)flag=1;
        if(flag==1)continue;
*****/ 
        for(a=0;a<=20;++a)amco[a]=0;
        for(i=1;i<=n;++i)if(marks[i]==1)amco[ali[i][k]]++;
        flag=0;for(a=1;a<=20;++a)if(amco[a]>0){flag=1;lettercount++;}
        if(flag==1)sitecount++;
                               }
if(sitecount==0)letpersite=0;
else letpersite=1.0*lettercount/sitecount;

neff=-log(1.0-0.05*letpersite)/0.05129329438755;

return neff;
}


void subalign::pseudoCounts(double **matrix, double n_eff, int len, double **pseudoCnt)
{
        int i,j,k;
        double f[21], g[21];
        double sumN;
        double alpha;

        alpha = n_eff-1;
	//cout << n_eff << endl;
	//cout << alpha << endl;
	//cout << beta <<endl;

	if(beta==0) alpha = 1;

        for (i=1;i<=len;i++) {
                sumN = 0;
                for (j=1;j<=20;j++) sumN += matrix[i][j];
                for (j=1;j<=20;j++) {
                        f[j] = 1.0*matrix[i][j]/sumN;
                }
                for (j=1;j<=20;j++) {
                        g[j] = 0;
                        for (k=1;k<=20;k++) g[j]+= qmatrix[j][k]*f[k]/robinson_freq[k];
                        pseudoCnt[i][j]= (alpha*f[j] + beta*g[j])/(alpha+beta);
                }
        }

}

void subalign::profile() {

    	int i,j,k,m,n;

	apos_filtr = new int [alilen+1];
	n_effAa = new double * [alilen+1];
	pseudoCnt = new double * [alilen+1];
	double *tmp = new double [15];
	fflush(stdout);
	for(i=0;i<=alilen;i++) {
		n_effAa[i] = new double [20+1];
		pseudoCnt[i] = new double [20+1];
	}
	sum_eff_let = new double [alilen+1];
	maskgapRegion = new int [alilen+2];

	for(i=0;i<=20;i++) for(j=0;j<=20;j++) qmatrix[i][j] = q_blosum62[i][j];

	neffsForEachCol_maskGapReg(alignment, nal, alilen, gap_threshold, gapt_threshold, n_effAa, sum_eff_let, maskgapRegion, apos_filtr, &alilen_mat, &n_eff);

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
	double sum_hwt_all=0;

	for(j=1;j<=alilen_mat;j++) {
		position_num_letters = 0;
		for(i=1;i<=20;i++) tmp_mark[i] = 0;
		for(i=1;i<=nal;i++) {
			//cout << j << " " << i << " " << alignment[i][apos_filtr[j]] << endl;
			if(alignment[i][apos_filtr[j]]==0) continue; // ignore gaps
			if( tmp_mark[alignment[i][apos_filtr[j]]] ) {
				tmp_mark[alignment[i][apos_filtr[j]]]++;
			}
			else {
				position_num_letters+=1;
				tmp_mark[alignment[i][apos_filtr[j]]]++;
			}
		}
		//cout << "position_num_letters: " << position_num_letters << endl;
		//for(i=1;i<=20;i++) cout << "tmp_mark " << i << " " << tmp_mark[i] << endl;
		for(i=1;i<=nal;i++) {
			if(alignment[i][apos_filtr[j]]==0) continue; // ignore gaps
			hwt_all[i] += 1.0 / position_num_letters / tmp_mark[alignment[i][apos_filtr[j]]];
		}
	}
	for(i=1;i<=nal;i++) {sum_hwt_all += hwt_all[i];}
	for(i=1;i<=nal;i++) {hwt_all[i]/=sum_hwt_all;}
	//cout << "sum_hwt_all: " << sum_hwt_all << " alilen_mat: " << alilen_mat << endl;
}
	

void subalign::printProfile() {
	
	int i, j;
	cout << "nal: " << nal << endl;
	cout << "alilen: " << alilen << endl;
	cout << "gap_threshold: " << gap_threshold << endl;
	cout << "gapt_threshold: " << gapt_threshold << endl;
	cout << "alilen_mat: " << alilen_mat << endl;
	cout << "n_eff: " << n_eff << endl;

	for(i=1;i<=alilen;i++) {
		cout << "maskgapRegion apos_filtr sum_eff_let " << i << " " << maskgapRegion[i] << " " << apos_filtr[i] << " " << sum_eff_let[i] << endl;
	}
}
	

// Construct a sub alignment with sequences marked in array mark[]
subalign *subalign::sub2align(int *mark) {
	
	int i,j,k;
	int nal1=0, alilen1;
	char **aseq1, **aname1;
	int **alignment1;
	int mlen;
	int count;


	subalign *aln1 = new subalign();

	for(i=1;i<=nal;i++) {
	     if(mark[i]) nal1++;
	}
	aln1->nal = nal1;
	aln1->alilen = alilen;

	mlen = 0;
	for (i=0; i < nal; i++) {
	        if(mark[i+1]) if (mlen < strlen(aname[i])) {
                        mlen = strlen(aname[i]);
                }
        }
	aln1->mnamelen = mlen;
	
	aln1->aseq = cmatrix(nal1, alilen);
	aln1->aname = cmatrix(nal1, mnamelen);
	aln1->alignment = imatrix(nal1, alilen);
	
	count = 0;
	for(i=1;i<=nal;i++) 	{
	   if(mark[i]) {
		strcpy(aln1->aseq[count], aseq[i-1]);
		strcpy(aln1->aname[count], aname[i-1]);
		for(j=1;j<=alilen;j++) {
			aln1->alignment[count+1][j] = alignment[i][j];
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


// Construct a sub alignment with sequences marked in array mark[1..nal] mark1[1..alilen]
subalign *subalign::sub2align(int *mark, int *mark1) {
	
	int i,j,k;
	int nal1=0, alilen1=0;
	char **aseq1, **aname1;
	int **alignment1;
	int mlen;
	int count, countl;


	subalign *aln1 = new subalign();

	// sequence number
	for(i=1;i<=nal;i++) {
	     if(mark[i]) nal1++;
	}
	aln1->nal = nal1;
	// sequence length
	for(i=1;i<=alilen;i++) {
	     if(mark1[i]) alilen1++;
	}
	aln1->alilen = alilen1;
	// maximum name length
	mlen = 0;
        for (i=0; i < nal; i++) {
                if(mark[i+1]) if (mlen < strlen(aname[i])) {
                        mlen = strlen(aname[i]);
                }
        }
        aln1->mnamelen = mlen;	

	aln1->aseq = cmatrix(nal1, alilen1+1);
	aln1->aname = cmatrix(nal1, mlen);
	aln1->alignment = imatrix(nal1, alilen1);
	
	count = 0;
	for(i=1;i<=nal;i++) 	{
	   if(mark[i]) {
		countl = 0;
		for(j=1;j<=alilen;j++) {
		    if(mark1[j]) {
			aln1->aseq[count][countl] = aseq[i-1][j-1];
			//cout  <<"|"<< aseq1[count][countl] <<"|";
			aln1->alignment[count+1][countl+1] = alignment[i][j];
		        countl++;
		    }
		}
		//cout << endl;
		//cout << alilen1 << "\t" << countl << endl;
		aln1->aseq[count][aln1->alilen] = '\0';
		//if(mark1[j-1]) aseq1[count][countl]='\0';
		//else aseq1[count][countl+1]='\0';
		//cout << aname[i-1] << "\t" << aseq1[count] << "|" << endl;
		strcpy(aln1->aname[count], aname[i-1]);
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

// convert the alignment in letters to alignment in numbers
void subalign::convertAseq2Alignment() {

	int i,j;
	
	if(!alignment) alignment = imatrix(nal, alilen);

	for(i=1;i<=nal;i++) {
	   for(j=1;j<=alilen;j++) {
		alignment[i][j] = am2num(aseq[i-1][j-1]);
	   }
	}
	
}

void subalign::add_sequence(char *name, char *seq) {

	int i,j;

	if(nal == 0) {
		nal++;
		alilen = strlen(seq);
		mnamelen = strlen(name)+1;

		aname = cmatrix(nal, mnamelen);
		aseq = cmatrix(nal, alilen);
		strcpy(aname[nal-1], name);
		strcpy(aseq[nal-1], seq);
		alignment = imatrix(nal, alilen);
		for(i=1;i<=nal;i++) {
	  	    for(j=1;j<=alilen;j++) {
			alignment[i][j] = am2num(aseq[i-1][j-1]);
		    }
		}
	}

	else {
		nal++;
		if(mnamelen < strlen(name)+1) mnamelen = strlen(name)+1;
		char **new_aname = cmatrix(nal, mnamelen);
		char **new_aseq = cmatrix(nal, alilen);
		for(i=0;i<nal-1;i++) {
			strcpy(new_aname[i], aname[i]);
			strncpy(new_aseq[i], aseq[i], alilen);
			new_aseq[i][alilen] = '\0';
			//cout << new_aname[i] << "\t"<<new_aseq[i]<<endl;
		}
		strcpy(new_aname[nal-1], name);
		strncpy(new_aseq[nal-1], seq, alilen);
		new_aseq[nal-1][alilen] = '\0';
		//cout << new_aname[nal-1] << "\t"<<new_aseq[nal-1]<<endl;
		
        	int **new_alignment = imatrix(nal, alilen);

	        for(i=1;i<=nal;i++) {
       	    		for(j=1;j<=alilen;j++) {
                		new_alignment[i][j] = am2num(new_aseq[i-1][j-1]);
           		}
        	}

		for(i=0;i<nal-1;i++) {
			delete [] aseq[i];
			delete [] aname[i];
		}
		for(i=0;i<=nal-1;i++) delete [] alignment[i];
		delete [] aseq;
		delete [] aname;
		delete [] alignment;

		aseq = new_aseq;	
		aname = new_aname;
		alignment = new_alignment;
		//cout << "===========" << endl;
	}
}

	

subalign * oneSeq2subalign(char *seq, char *name) {
	
	int i,j;

	subalign *a = new subalign();
	
	a->aseq = cmatrix(1, strlen(seq)+1);
	a->aname = cmatrix(1, strlen(name)+1);
	a->alignment = imatrix(1, strlen(seq) );
	for(i=0;i<=strlen(seq);i++) a->aseq[0][i] = seq[i];
	for(i=0;i<=strlen(name);i++) a->aname[0][i] = name[i];
	for(i=1;i<=strlen(seq);i++) a->alignment[1][i] = am2num( seq[i-1]);
	a->nal = 1;
	a->alilen = strlen(seq);
	a->mnamelen = strlen(name);

	//cout << a->alilen << endl;

	return a;
	
}

// purge an alignment by
// 1. removing highly similar sequences (id > high_id_thr)
// 2. removing highly divergent sequences to the first sequence (id < low_id_thr)
// 3. removing sequence fragment relative to the first sequence (gap_fraction > gap_fraction_thr)
// 4. selecting only a subset if remaining num > max_num_kept
// *** all positions are kept ***
subalign * subalign::purge_align(double low_id_thr, double high_id_thr, int max_num_kept, double gap_fraction_thr) {

	int i, j, k;
	int pairs, id_pairs;
	double gap_fraction;

	int *mark = ivector(nal);
	for(i=1;i<=nal;i++) mark[i] = 1;

	assert(low_id_thr<=1);
	assert(high_id_thr<=1);

	if(low_id_thr>0)
	for(i=1;i<nal;i++) {
	    pairs = 0; id_pairs = 0;
	    for(k=0;k<alilen;k++) {
		if( (aseq[i][k]!='-') && (aseq[0][k]!='-') ) {
		    pairs += 1;
		    if(aseq[i][k]==aseq[0][k]) id_pairs += 1;
		}
	    }
	    if(1.0*id_pairs/pairs < low_id_thr) mark[i+1] = 0;
	    if(mark[i+1]==0) continue;

	    pairs = 0; id_pairs = 0;
	    for(k=0;k<alilen;k++) {
		if(aseq[0][k]!='-') {
		    pairs += 1;
		    if(aseq[i][k]=='-') id_pairs += 1;
		}
	    }
	    gap_fraction = 1.0 *  id_pairs / pairs;
	    if(gap_fraction > gap_fraction_thr) mark[i+1] = 0;
	}

	for(i=0;i<nal;i++) {
	    if(mark[i+1]==0) continue;
	    for(j=i+1;j<nal;j++) {
		if (mark[j+1]==0) { continue; }
		pairs = 0; id_pairs = 0;
		for(k=0;k<alilen;k++) {
		    if( (aseq[i][k]!='-') && (aseq[j][k]!='-') ) {
			pairs += 1;
			if(aseq[i][k]==aseq[j][k]) id_pairs += 1;
		    }
		}
		if(1.0*id_pairs/pairs > high_id_thr) mark[j+1] = 0;
	    }
	}

	if(mydebug>100)for(i=1;i<=nal;i++) cout << i << "\t" << mark[i] << endl;

	subalign *xaln = sub2align(mark);
	delete [] mark;
	return xaln;

}
	
void subalign::prof_h_weight_all(double raw_gap_fraction_cutoff) {

	int i, j, k;
	int position_num_letters;
	int tmp_mark[21];
	prof_hwt_all = dvector(nal);
	for(i=1;i<=nal;i++) prof_hwt_all[i] = 1e-06;
	double sum_hwt_all=0;
	int raw_gap_counts;

	for(j=1;j<=alilen;j++) {
		raw_gap_counts = 0;
		for(i=1;i<=nal;i++) {
			if(!alignment[i][j]) raw_gap_counts++;
		}
		if(raw_gap_counts*1.0/nal > raw_gap_fraction_cutoff) continue;
		//cout << "j: " << j << endl;
		position_num_letters = 0;
		for(i=0;i<=20;i++) tmp_mark[i] = 0;
		// *** IMPORTANT ***
		// gap character is not considered as the 21 letter
		// this is to avoid giving fragments unreasonably high weights
		for(i=1;i<=nal;i++) {
			//cout << j << " " << i << " " << alignment[i][apos_filtr[j]] << endl;
			if(alignment[i][j]==0) continue; // ignore gaps
			if( tmp_mark[alignment[i][j]] ) {
				tmp_mark[alignment[i][j]]++;
			}
			else {
				position_num_letters+=1;
				tmp_mark[alignment[i][j]]++;
			}
		}
		//cout << "position_num_letters: " << position_num_letters << endl;
		//for(i=1;i<=20;i++) cout << "tmp_mark " << i << " " << tmp_mark[i] << endl;
		for(i=1;i<=nal;i++) {
			if(alignment[i][j]==0) continue; // ignore gaps
			prof_hwt_all[i] += 1.0 / position_num_letters / tmp_mark[alignment[i][j]];
		}
	}
	for(i=1;i<=nal;i++) {sum_hwt_all += prof_hwt_all[i];}
	for(i=1;i<=nal;i++) {prof_hwt_all[i]/=sum_hwt_all;}
	//cout << "sum_hwt_all: " << sum_hwt_all << " alilen_mat: " << alilen_mat << endl;

	if(mydebug>1) {
		cout << "henikoff weights: " << endl;
		for(i=1;i<=nal;i++) {
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
	for(i=1;i<=alilen;i++) {
		prof_gap_content[i] = 0;
		for(j=1;j<nal;j++) {
			if(alignment[j][i]==0) prof_gap_content[i] += prof_hwt_all[j];
		}
		if(prof_gap_content[i] < prof_gap_threshold_here) {
			prof_len++;
			prof_pos[prof_len] = i;
		}
	}

	if(prof_len == 0) {
		cout << "Warning: number of positions with weighted gap fraction less than " << prof_gap_threshold_here << " is zero" << endl;
		cout << "    - Use all the positions for profile" << endl;
		cout << "    - Each sequence has the same weight" << endl;
		cout << endl;
		for(i=1;i<alilen;i++) {
			prof_pos[i] = i;
		}
		prof_len = alilen;
	}
		
	if(mydebug>1) {
		cout << "prof positions: " << endl;
		cout << "prof length: " << prof_len << endl;
		for(i=1;i<=prof_len;i++) {
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
	
	prof_get_effn(alignment, nal, alilen, prof_len, prof_effn, prof_sum_eff, prof_pos, &prof_nef);

	pseudoCounts(prof_effn, prof_nef, prof_len, prof_freq);

	done_prof = 1;
}
	
/* from a given alignment with aa as numbers, computes effective aa counts (PSIC->our formula)
and marks the columns with EFFECTIVE content of gaps > threshold (effgapmax) */
void subalign::prof_get_effn(int **ali, int n, int len, int prof_len, double **n_effAa, double *sum_eff_let, int *prof_pos, double *nef)
{
        int i,j,k,l;
        int nsymbols_col, nsymbols;
        double nef_loc;
        int ele;
        double effnu[21];
        double sum_let;
        int *mark;
        int flagmark;

	int p;

	gap_content = new double [alilen+1];

	cout << "prof_len: " << prof_len << endl;

	mark = new int [n+11];
        nsymbols = 0;
        for(j=1;j<=prof_len;j++) {

		p = prof_pos[j];

                nsymbols_col = 0;
                sum_let=0;
		if(mydebug>100) cout << "p: " << p << endl;

                for(k=0;k<=20;++k){
/* Mark sequences that have amino acid  k (or gap, k=0) in this jth position */
                         flagmark =0;
                        for(i=1;i<=n;++i){
                                mark[i]=0;

                                ele=ali[i][p];
                                if(ele==k){mark[i]=1; flagmark =1;}
                                ele=ali[i][p]-25;
                                if(ele==k) {mark[i]=1; flagmark =1;}
                        }

/* If aa k (or gap) is present in this position call compute k-th effective count */
                      if (flagmark == 1) {

                                effnu[k]=effective_number_nogaps(ali,mark,n,1,len);
                                nsymbols_col++;

                        } else { effnu[k] = 0.0; }

                       if (k>0) sum_let += effnu[k];
                }

		// calculate raw gap fraction
		if(mydebug>100) {
			int gap_count = 0;
			for(i=1;i<=n;i++) {
			 if(ali[i][p]==0) gap_count++;
			}
			//double gapfract = 1.0 * gap_count / n;
			//cout << j << "  " << 1.0*effnu[0]/(sum_let + effnu[0]) << "  " << gapfract << endl;
		}

		if(mydebug>100) cout << j << "  " << aseq[0][j-1] << " " << sum_let << "  " << effnu[0] << endl;
		for(k=0; k<=20; k++) {
			n_effAa[j][k] = effnu[k];
			nsymbols += nsymbols_col;
		}
		if(mydebug>100) cout << "nsymbols: " << nsymbols << endl;

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

	//cout << nsymbols << "\t" << alilen_mat << endl;

        nef_loc = 1.0*nsymbols/prof_len;
        *nef = nef_loc;

	prof_average_sum_eff_let = 0;
	for(i=1;i<=prof_len;i++) {
		prof_average_sum_eff_let += sum_eff_let[i];
	}
	average_sum_eff_let /= prof_len;

        delete [] mark;

	/* for(i=1;i<=alilen_mat;i++) {
	    for(j=1;j<=20;j++) {
		cout << n_effAa[i][j] << " ";
	    }
	    cout << endl;
	} */

}


