#include "header_cpp.h"
#include "amino.h"

static int debug = 1;

#include <cctype>
int a3let2num(char *let);
int am2num_c(int c);
int am2num(int c);
int am2num_dayhoff(int c);
int am2numBZX(int c);
int am2nm(int c);
double dayhoff_freq[]={0,0.010,0.040,0.030,0.015,0.085,0.037,0.065,0.087,0.033,0.089,
0.051,0.058,0.070,0.040,0.038,0.047,0.050,0.034,0.041,0.081};
double robinson_freq[]={0.00000,0.01330,0.03856,0.03216,0.02243,0.09019,0.05142,0.06441,0.07805,0.01925,0.07377,0.05203,0.05841,0.07120,0.04487,0.04264,0.05364,0.06295,0.02199,0.05130,0.05744};
/*in dayhoff_freq[] amino acids are arranged as in *am */
int dayhoff_mutab[]={0,18,41,41,94,40,96,74,100,20,49,56,97,120,134,93,106,102,66,65,56};
/*in dayhoff_mutab[] amino acids are arranged as in *am */
double hydrophobicity[]={0,1.17,1.56,1.17,0.72,1.38,1.38,1.00,0.33,0.53,-0.04,0.63,-0.21,-0.50,-0.89,-0.84,-1.56,-1.51,-1.12,-1.65,-0.74};
/* double hydrophobicity[]={0,0.878,1.00,0.88,0.738,0.943,0.943,0.825,0.616,0.680,0.501,0.711,0.45,0.359,0.236,0.251,0.028,0.043,0.165,0.0,0.283}; */
/*in hydrophobicity[] amino acids are arranged as in *am */
double amino_charge[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.0,-1.0,0.5,1.0,1.0};
/*in amino_charge[] amino acids are arranged as in *am */

// logarithm robinsion frequencies
float log_robinson_freq[21];
void get_log_robinson_freq() {
	int i;
	log_robinson_freq[0] = -100000;
	for(i=1;i<=20;i++) log_robinson_freq[i] = log(robinson_freq[i]);
	if(debug>1) for(i=1;i<=20;i++) { cout << "i: " << i << "\t" << log_robinson_freq[i] << endl; }

}

const char *am="-WFYMLIVACGPTSNQDEHRKBZX*.wfymlivacgptsnqdehrkbzx";

const char *am_dayhoff="-ARNDCQEGHILKMFPSTWYV";
const char *am3[]={
"---",
"TRP",
"PHE",
"TYR",
"MET",
"LEU",
"ILE",
"VAL",
"ALA",
"CYS",
"GLY",
"PRO",
"THR",
"SER",
"ASN",
"GLN",
"ASP",
"GLU",
"HIS",
"ARG",
"LYS",
"ASX",
"GLX",
"UNK",
"***"
"...",
};



int am2num_c(int c)
{
switch (c) {
           	 case 'W':
                	c=1; break;
		 case 'w':
			c=26; break;
           	 case 'F': 
                	c=2; break;
                 case 'f':
			c=27; break;
           	 case 'Y': 
                	c=3; break;
		 case 'y':
                        c=28; break;
           	 case 'M': 
                	c=4; break;
		 case 'm':
                        c=29; break;
           	 case 'L': 
                	c=5; break;
		 case 'l':
                        c=30; break;
           	 case 'I': 
          		c=6; break;
		 case 'i':
                        c=31; break;
           	 case 'V': 
           		c=7; break;
		 case 'v':
                        c=32; break;
          	 case 'A': 
			c=8; break;
		 case 'a':
                        c=33; break; 
           	 case 'C': 
                	c=9; break;
		 case 'c':
                        c=34; break;
		 case 'G': 
			c=10; break;
		 case 'g':
                        c=35; break;
           	 case 'P': 
             	 	c=11; break;
		 case 'p':
                        c=36; break;
       		 case 'T': 
			c=12; break;
		 case 't':
                        c=37; break;
	         case 'S': 
			c=13; break;
		 case 's':
                        c=38; break;
           	 case 'N': 
                	c=14; break;
		 case 'n':
                        c=39; break;
           	 case 'Q': 
                	c=15; break;
		 case 'q':
                        c=40; break;
           	 case 'D': 
                	c=16; break;
		 case 'd':
                        c=41; break;
           	 case 'E': 
                	c=17; break;
		 case 'e':
                        c=42; break;
           	 case 'H': 
                	c=18; break;
		 case 'h':
                        c=43; break;
           	 case 'R': 
                	c=19; break;
		 case 'r':
                        c=44; break;
           	 case 'K': 
                	c=20; break;
		 case 'k':
                        c=45; break;
           	 default : 
			c=0; 
		}
return (c);
}


int am2num(int c)
{
switch (c) {
           	 case 'W': case 'w':
                	c=1; break;
           	 case 'F': case 'f':
                	c=2; break;
           	 case 'Y': case 'y':
                	c=3; break;
           	 case 'M': case 'm':
                	c=4; break;
           	 case 'L': case 'l':
                	c=5; break;
           	 case 'I': case 'i':
          		c=6; break;
           	 case 'V': case 'v':
           		c=7; break;
          	 case 'A': case 'a': 
			c=8; break;
           	 case 'C': case 'c':
                	c=9; break;
		 case 'G': case 'g':
			c=10; break;
           	 case 'P': case 'p':
             	 	c=11; break;
       		 case 'T': case 't':
			c=12; break;
	         case 'S': case 's':
			c=13; break;
           	 case 'N': case 'n':
                	c=14; break;
           	 case 'Q': case 'q':
                	c=15; break;
           	 case 'D': case 'd':
                	c=16; break;
           	 case 'E': case 'e':
                	c=17; break;
           	 case 'H': case 'h':
                	c=18; break;
           	 case 'R': case 'r':
                	c=19; break;
           	 case 'K': case 'k':
                	c=20; break;
		 // dealing with X, or ambiguous amino acids: treat as 'A'
		 case 'X': case 'x':
		 case 'B': case 'b':
		 case 'J': case 'j':
		 case 'O': case 'o':
		 case 'U': case 'u':
		 case 'Z': case 'z':
			c=8; break;
           	 default : 
			c=0; 
		}
return (c);
}


int am2num_dayhoff(int c)
{
switch (c) {
           	 case 'A':
                	c=1; break;
           	 case 'R':
                	c=2; break;
           	 case 'N':
                	c=3; break;
           	 case 'D':
                	c=4; break;
           	 case 'C': 
                	c=5; break;
           	 case 'Q':
          		c=6; break;
           	 case 'E':
           		c=7; break;
          	 case 'G':
			c=8; break;
           	 case 'H': 
                	c=9; break;
		 case 'I':
			c=10; break;
           	 case 'L':
             	 	c=11; break;
       		 case 'K':
			c=12; break;
	         case 'M':
			c=13; break;
           	 case 'F':
                	c=14; break;
           	 case 'P':
                	c=15; break;
           	 case 'S':
                	c=16; break;
           	 case 'T': 
                	c=17; break;
           	 case 'W': 
                	c=18; break;
           	 case 'Y':
                	c=19; break;
           	 case 'V':
                	c=20; break;
           	 default : 
			c=0; 
		}
return (c);
}

int am2numBZX(int c)
{
switch (c) {
           	 case 'W': case 'w':
                	c=1; break;
           	 case 'F': case 'f':
                	c=2; break;
           	 case 'Y': case 'y':
                	c=3; break;
           	 case 'M': case 'm':
                	c=4; break;
           	 case 'L': case 'l':
                	c=5; break;
           	 case 'I': case 'i':
          		c=6; break;
           	 case 'V': case 'v':
           		c=7; break;
          	 case 'A': case 'a': 
			c=8; break;
           	 case 'C': case 'c':
                	c=9; break;
		 case 'G': case 'g':
			c=10; break;
           	 case 'P': case 'p':
             	 	c=11; break;
       		 case 'T': case 't':
			c=12; break;
	         case 'S': case 's':
			c=13; break;
           	 case 'N': case 'n':
                	c=14; break;
           	 case 'Q': case 'q':
                	c=15; break;
           	 case 'D': case 'd':
                	c=16; break;
           	 case 'E': case 'e':
                	c=17; break;
           	 case 'H': case 'h':
                	c=18; break;
           	 case 'R': case 'r':
                	c=19; break;
           	 case 'K': case 'k':
                	c=20; break;
           	 case 'B': case 'b':
                	c=21; break;
           	 case 'Z': case 'z':
                	c=22; break;
           	 case 'X': case 'x':
                	c=23; break;
           	 case '*': 
                	c=24; break;
           	 default : 
			c=0; 
		}
return (c);
}

int am2nm(int c) /* am2numBZX_c */
{
switch (c) {
		 case 'W':
                        c=1; break;
                 case 'w':
                        c=26; break;
                 case 'F':
                        c=2; break;
                 case 'f':
                        c=27; break;
                 case 'Y':
                        c=3; break;
                 case 'y':
                        c=28; break;
                 case 'M':
                        c=4; break;
                 case 'm':
                        c=29; break;
                 case 'L':
                        c=5; break;
                 case 'l':
                        c=30; break;
                 case 'I':
                        c=6; break;
                 case 'i':
                        c=31; break;
                 case 'V':
                        c=7; break;
                 case 'v':
                        c=32; break;
                 case 'A':
                        c=8; break;
                 case 'a':
                        c=33; break;
                 case 'C':
                        c=9; break;
                 case 'c':
                        c=34; break;
                 case 'G':
                        c=10; break;
                 case 'g':
                        c=35; break;
                 case 'P':
                        c=11; break;
                 case 'p':
                        c=36; break;
                 case 'T':
                        c=12; break;
                 case 't':
                        c=37; break;
                 case 'S':
                        c=13; break;
                 case 's':
                        c=38; break;
                 case 'N':
                        c=14; break;
                 case 'n':
                        c=39; break;
                 case 'Q':
                        c=15; break;
                 case 'q':
                        c=40; break;
                 case 'D':
                        c=16; break;
                 case 'd':
                        c=41; break;
                 case 'E':
                        c=17; break;
                 case 'e':
                        c=42; break;
                 case 'H':
                        c=18; break;
                 case 'h':
                        c=43; break;
                 case 'R':
                        c=19; break;
                 case 'r':
                        c=44; break;
                 case 'K':
                        c=20; break;
                 case 'k':
                        c=45; break;
		 case 'B':
			c=21; break;
		 case 'b':
			c=46; break;
		 case 'Z':
			c=22; break;
		 case 'z':
			c=47; break;
		 case 'X':
			c=23; break;
		 case 'x':
			c=48; break;
                 default :
                        c=0;
                }
return (c);
}


int a3let2num(char *let){

char s[3];
s[0]=toupper(let[0]);
s[1]=toupper(let[1]);
s[2]=toupper(let[2]);

     if(s[0]=='T' && s[1]=='R' && s[2]=='P')return(1);
else if(s[0]=='P' && s[1]=='H' && s[2]=='E')return(2);
else if(s[0]=='T' && s[1]=='Y' && s[2]=='R')return(3);
else if(s[0]=='M' && s[1]=='E' && s[2]=='T')return(4);
else if(s[0]=='L' && s[1]=='E' && s[2]=='U')return(5);
else if(s[0]=='I' && s[1]=='L' && s[2]=='E')return(6);
else if(s[0]=='V' && s[1]=='A' && s[2]=='L')return(7);
else if(s[0]=='A' && s[1]=='L' && s[2]=='A')return(8);
else if(s[0]=='C' && s[1]=='Y' && s[2]=='S')return(9);
else if(s[0]=='G' && s[1]=='L' && s[2]=='Y')return(10);
else if(s[0]=='P' && s[1]=='R' && s[2]=='O')return(11);
else if(s[0]=='T' && s[1]=='H' && s[2]=='R')return(12);
else if(s[0]=='S' && s[1]=='E' && s[2]=='R')return(13);
else if(s[0]=='A' && s[1]=='S' && s[2]=='N')return(14);
else if(s[0]=='G' && s[1]=='L' && s[2]=='N')return(15);
else if(s[0]=='A' && s[1]=='S' && s[2]=='P')return(16);
else if(s[0]=='G' && s[1]=='L' && s[2]=='U')return(17);
else if(s[0]=='H' && s[1]=='I' && s[2]=='S')return(18);
else if(s[0]=='A' && s[1]=='R' && s[2]=='G')return(19);
else if(s[0]=='L' && s[1]=='Y' && s[2]=='S')return(20);
else if(s[0]=='A' && s[1]=='S' && s[2]=='X')return(21);
else if(s[0]=='G' && s[1]=='L' && s[2]=='X')return(22);
else if(s[0]=='U' && s[1]=='N' && s[2]=='K')return(23);
else if(s[0]=='*' && s[1]=='*' && s[2]=='*')return(24);
else if(s[0]=='.' && s[1]=='.' && s[2]=='.')return(25);
return(0);
}

/*
Matrix file format:

#  Matrix made by matblas from blosum45.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/3 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 45
#  Entropy =   0.3795, Expected =  -0.2789
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  5 -2 -1 -2 -1 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -2 -2  0 -1 -1  0 -5
R -2  7  0 -1 -3  1  0 -2  0 -3 -2  3 -1 -2 -2 -1 -1 -2 -1 -2 -1  0 -1 -5
N -1  0  6  2 -2  0  0  0  1 -2 -3  0 -2 -2 -2  1  0 -4 -2 -3  4  0 -1 -5
D -2 -1  2  7 -3  0  2 -1  0 -4 -3  0 -3 -4 -1  0 -1 -4 -2 -3  5  1 -1 -5
C -1 -3 -2 -3 12 -3 -3 -3 -3 -3 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -2 -3 -2 -5
Q -1  1  0  0 -3  6  2 -2  1 -2 -2  1  0 -4 -1  0 -1 -2 -1 -3  0  4 -1 -5
E -1  0  0  2 -3  2  6 -2  0 -3 -2  1 -2 -3  0  0 -1 -3 -2 -3  1  4 -1 -5
G  0 -2  0 -1 -3 -2 -2  7 -2 -4 -3 -2 -2 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -5
H -2  0  1  0 -3  1  0 -2 10 -3 -2 -1  0 -2 -2 -1 -2 -3  2 -3  0  0 -1 -5
I -1 -3 -2 -4 -3 -2 -3 -4 -3  5  2 -3  2  0 -2 -2 -1 -2  0  3 -3 -3 -1 -5
L -1 -2 -3 -3 -2 -2 -2 -3 -2  2  5 -3  2  1 -3 -3 -1 -2  0  1 -3 -2 -1 -5
K -1  3  0  0 -3  1  1 -2 -1 -3 -3  5 -1 -3 -1 -1 -1 -2 -1 -2  0  1 -1 -5
M -1 -1 -2 -3 -2  0 -2 -2  0  2  2 -1  6  0 -2 -2 -1 -2  0  1 -2 -1 -1 -5
F -2 -2 -2 -4 -2 -4 -3 -3 -2  0  1 -3  0  8 -3 -2 -1  1  3  0 -3 -3 -1 -5
P -1 -2 -2 -1 -4 -1  0 -2 -2 -2 -3 -1 -2 -3  9 -1 -1 -3 -3 -3 -2 -1 -1 -5
S  1 -1  1  0 -1  0  0  0 -1 -2 -3 -1 -2 -2 -1  4  2 -4 -2 -1  0  0  0 -5
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -1 -1  2  5 -3 -1  0  0 -1  0 -5
W -2 -2 -4 -4 -5 -2 -3 -2 -3 -2 -2 -2 -2  1 -3 -4 -3 15  3 -3 -4 -2 -2 -5
Y -2 -1 -2 -2 -3 -1 -2 -3  2  0  0 -1  0  3 -3 -2 -1  3  8 -1 -2 -2 -1 -5
V  0 -2 -3 -3 -1 -3 -3 -3 -3  3  1 -2  1  0 -3 -1  0 -3 -1  5 -3 -3 -1 -5
B -1 -1  4  5 -2  0  1 -1  0 -3 -3  0 -2 -3 -2  0  0 -4 -2 -3  4  2 -1 -5
Z -1  0  0  1 -3  4  4 -2  0 -3 -2  1 -1 -3 -1  0 -1 -2 -2 -3  2  4 -1 -5
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -2 -1 -1 -1 -1 -1 -5
* -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1
*/

void read_aa_matrix(char *filename, int **arr)
{
	int i;
	char str[100];
	char amino[25][10];
	
	ifstream FP(filename);
	if(!FP) {cout << "not readable matrix file" << endl; exit(0);}

	int indi = 0;
	while(FP) {
	  FP >> str;
	  if(strcmp(str, "#") == 0) {
		FP.getline(str, 255);
		continue;
	  }
	  
	  if(strcmp(str, "A") == 0 && indi == 0 ) {
		for(i=2;i<=24;i++) {
		    FP >> amino[i];
		}
		strcpy(amino[0], "-");
		strcpy(amino[1], "A");
		indi = 1;
		continue;
	  }

	  if(strcmp(str, "*")) {
	    for(i=1;i<24;i++) {
	       FP >> arr[am2numBZX(str[0])][am2numBZX(amino[i][0])];
	    }
	    FP >> arr[am2numBZX(str[0])][0];
	  }
	  else {
	    for(i=1;i<24;i++) {
		FP >> arr[0][am2numBZX(amino[i][0])];
	    }
	    FP >> arr[0][0];
	  }
	}
	    
}

double q_blosum62[21][21] = {
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0, 0.0065, 0.0008, 0.0009, 0.0002, 0.0007, 0.0004, 0.0004, 0.0004, 0.0001, 0.0004, 0.0001, 0.0003, 0.0003, 0.0002, 0.0002, 0.0002, 0.0003, 0.0002, 0.0003, 0.0003 },
{0, 0.0008, 0.0183, 0.0042, 0.0012, 0.0054, 0.0030, 0.0026, 0.0016, 0.0005, 0.0012, 0.0005, 0.0012, 0.0012, 0.0008, 0.0005, 0.0008, 0.0009, 0.0008, 0.0009, 0.0009 },
{0, 0.0009, 0.0042, 0.0102, 0.0006, 0.0022, 0.0014, 0.0015, 0.0013, 0.0003, 0.0008, 0.0005, 0.0009, 0.0010, 0.0007, 0.0007, 0.0006, 0.0009, 0.0015, 0.0009, 0.0010 },
{0, 0.0002, 0.0012, 0.0006, 0.0040, 0.0049, 0.0025, 0.0023, 0.0013, 0.0004, 0.0007, 0.0004, 0.0010, 0.0009, 0.0005, 0.0007, 0.0005, 0.0007, 0.0004, 0.0008, 0.0009 },
{0, 0.0007, 0.0054, 0.0022, 0.0049, 0.0371, 0.0114, 0.0095, 0.0044, 0.0016, 0.0021, 0.0014, 0.0033, 0.0024, 0.0014, 0.0016, 0.0015, 0.0020, 0.0010, 0.0024, 0.0025 },
{0, 0.0004, 0.0030, 0.0014, 0.0025, 0.0114, 0.0184, 0.0120, 0.0032, 0.0011, 0.0014, 0.0010, 0.0027, 0.0017, 0.0010, 0.0009, 0.0012, 0.0012, 0.0006, 0.0012, 0.0016 },
{0, 0.0004, 0.0026, 0.0015, 0.0023, 0.0095, 0.0120, 0.0196, 0.0051, 0.0014, 0.0018, 0.0012, 0.0036, 0.0024, 0.0012, 0.0012, 0.0013, 0.0017, 0.0006, 0.0016, 0.0019 },
{0, 0.0004, 0.0016, 0.0013, 0.0013, 0.0044, 0.0032, 0.0051, 0.0215, 0.0016, 0.0058, 0.0022, 0.0037, 0.0063, 0.0019, 0.0019, 0.0022, 0.0030, 0.0011, 0.0023, 0.0033 },
{0, 0.0001, 0.0005, 0.0003, 0.0004, 0.0016, 0.0011, 0.0014, 0.0016, 0.0119, 0.0008, 0.0004, 0.0009, 0.0010, 0.0004, 0.0003, 0.0004, 0.0004, 0.0002, 0.0004, 0.0005 },
{0, 0.0004, 0.0012, 0.0008, 0.0007, 0.0021, 0.0014, 0.0018, 0.0058, 0.0008, 0.0378, 0.0014, 0.0022, 0.0038, 0.0029, 0.0014, 0.0025, 0.0019, 0.0010, 0.0017, 0.0025 },
{0, 0.0001, 0.0005, 0.0005, 0.0004, 0.0014, 0.0010, 0.0012, 0.0022, 0.0004, 0.0014, 0.0191, 0.0014, 0.0017, 0.0009, 0.0008, 0.0012, 0.0014, 0.0005, 0.0010, 0.0016 },
{0, 0.0003, 0.0012, 0.0009, 0.0010, 0.0033, 0.0027, 0.0036, 0.0037, 0.0009, 0.0022, 0.0014, 0.0125, 0.0047, 0.0022, 0.0014, 0.0019, 0.0020, 0.0007, 0.0018, 0.0023 },
{0, 0.0003, 0.0012, 0.0010, 0.0009, 0.0024, 0.0017, 0.0024, 0.0063, 0.0010, 0.0038, 0.0017, 0.0047, 0.0126, 0.0031, 0.0019, 0.0028, 0.0030, 0.0011, 0.0023, 0.0031 },
{0, 0.0002, 0.0008, 0.0007, 0.0005, 0.0014, 0.0010, 0.0012, 0.0019, 0.0004, 0.0029, 0.0009, 0.0022, 0.0031, 0.0141, 0.0015, 0.0037, 0.0022, 0.0014, 0.0020, 0.0024 },
{0, 0.0002, 0.0005, 0.0007, 0.0007, 0.0016, 0.0009, 0.0012, 0.0019, 0.0003, 0.0014, 0.0008, 0.0014, 0.0019, 0.0015, 0.0073, 0.0016, 0.0035, 0.0010, 0.0025, 0.0031 },
{0, 0.0002, 0.0008, 0.0006, 0.0005, 0.0015, 0.0012, 0.0013, 0.0022, 0.0004, 0.0025, 0.0012, 0.0019, 0.0028, 0.0037, 0.0016, 0.0213, 0.0049, 0.0010, 0.0016, 0.0024 },
{0, 0.0003, 0.0009, 0.0009, 0.0007, 0.0020, 0.0012, 0.0017, 0.0030, 0.0004, 0.0019, 0.0014, 0.0020, 0.0030, 0.0022, 0.0035, 0.0049, 0.0161, 0.0014, 0.0027, 0.0041 },
{0, 0.0002, 0.0008, 0.0015, 0.0004, 0.0010, 0.0006, 0.0006, 0.0011, 0.0002, 0.0010, 0.0005, 0.0007, 0.0011, 0.0014, 0.0010, 0.0010, 0.0014, 0.0093, 0.0012, 0.0012 },
{0, 0.0003, 0.0009, 0.0009, 0.0008, 0.0024, 0.0012, 0.0016, 0.0023, 0.0004, 0.0017, 0.0010, 0.0018, 0.0023, 0.0020, 0.0025, 0.0016, 0.0027, 0.0012, 0.0178, 0.0062 },
{0, 0.0003, 0.0009, 0.0010, 0.0009, 0.0025, 0.0016, 0.0019, 0.0033, 0.0005, 0.0025, 0.0016, 0.0023, 0.0031, 0.0024, 0.0031, 0.0024, 0.0041, 0.0012, 0.0062, 0.0161 }
};

float log_q_blosum62_ratio[21][21]; 

void get_log_q_blosum62_ratio() {

	int i,j;

	for(i=0;i<=20;i++) {
	     for(j=0;j<=20;j++) {
	        if(!q_blosum62[i][j]) log_q_blosum62_ratio[i][j] = -10000;
		else log_q_blosum62_ratio[i][j] = log(q_blosum62[i][j]/robinson_freq[i]/robinson_freq[j] );
		if(debug>1) cout << log_q_blosum62_ratio[i][j] << " ";
	     }
	     if(debug>1) cout << endl;
	}
}

float log_q_blosum62[21][21]; 

void get_log_q_blosum62() {

	int i,j;

	for(i=0;i<=20;i++) {
	     for(j=0;j<=20;j++) {
	        if(!q_blosum62[i][j]) log_q_blosum62[i][j] = -10000;
		else log_q_blosum62[i][j] = log(q_blosum62[i][j]);
	     }
	}
}

// Below are blosum62 frequencies from probcons
string alphabetDefault = "ARNDCQEGHILKMFPSTWYV";
float emitSingleDefault[20] = {
  0.07831005f, 0.05246024f, 0.04433257f, 0.05130349f, 0.02189704f,
  0.03585766f, 0.05615771f, 0.07783433f, 0.02601093f, 0.06511648f,
  0.09716489f, 0.05877077f, 0.02438117f, 0.04463228f, 0.03940142f,
  0.05849916f, 0.05115306f, 0.01203523f, 0.03124726f, 0.07343426f
};

float emitPairsDefault[20][20] = {
  {0.02373072f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00244502f, 0.01775118f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00210228f, 0.00207782f, 0.01281864f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00223549f, 0.00161657f, 0.00353540f, 0.01911178f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00145515f, 0.00044701f, 0.00042479f, 0.00036798f, 0.01013470f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00219102f, 0.00253532f, 0.00158223f, 0.00176784f, 0.00032102f, 0.00756604f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00332218f, 0.00268865f, 0.00224738f, 0.00496800f, 0.00037956f, 0.00345128f, 0.01676565f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00597898f, 0.00194865f, 0.00288882f, 0.00235249f, 0.00071206f, 0.00142432f, 0.00214860f, 0.04062876f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00114353f, 0.00132105f, 0.00141205f, 0.00097077f, 0.00026421f, 0.00113901f, 0.00131767f, 0.00103704f, 0.00867996f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00318853f, 0.00138145f, 0.00104273f, 0.00105355f, 0.00094040f, 0.00100883f, 0.00124207f, 0.00142520f, 0.00059716f, 0.01778263f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00449576f, 0.00246811f, 0.00160275f, 0.00161966f, 0.00138494f, 0.00180553f, 0.00222063f, 0.00212853f, 0.00111754f, 0.01071834f, 0.03583921f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00331693f, 0.00595650f, 0.00257310f, 0.00252518f, 0.00046951f, 0.00312308f, 0.00428420f, 0.00259311f, 0.00121376f, 0.00157852f, 0.00259626f, 0.01612228f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00148878f, 0.00076734f, 0.00063401f, 0.00047808f, 0.00037421f, 0.00075546f, 0.00076105f, 0.00066504f, 0.00042237f, 0.00224097f, 0.00461939f, 0.00096120f, 0.00409522f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00165004f, 0.00090768f, 0.00084658f, 0.00069041f, 0.00052274f, 0.00059248f, 0.00078814f, 0.00115204f, 0.00072545f, 0.00279948f, 0.00533369f, 0.00087222f, 0.00116111f, 0.01661038f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00230618f, 0.00106268f, 0.00100282f, 0.00125381f, 0.00034766f, 0.00090111f, 0.00151550f, 0.00155601f, 0.00049078f, 0.00103767f, 0.00157310f, 0.00154836f, 0.00046718f, 0.00060701f, 0.01846071f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00631752f, 0.00224540f, 0.00301397f, 0.00285226f, 0.00094867f, 0.00191155f, 0.00293898f, 0.00381962f, 0.00116422f, 0.00173565f, 0.00250962f, 0.00312633f, 0.00087787f, 0.00119036f, 0.00180037f, 0.01346609f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00389995f, 0.00186053f, 0.00220144f, 0.00180488f, 0.00073798f, 0.00154526f, 0.00216760f, 0.00214841f, 0.00077747f, 0.00248968f, 0.00302273f, 0.00250862f, 0.00093371f, 0.00107595f, 0.00147982f, 0.00487295f, 0.01299436f, 0.0f, 0.0f, 0.0f},
  {0.00039119f, 0.00029139f, 0.00021006f, 0.00016015f, 0.00010666f, 0.00020592f, 0.00023815f, 0.00038786f, 0.00019097f, 0.00039549f, 0.00076736f, 0.00028448f, 0.00016253f, 0.00085751f, 0.00015674f, 0.00026525f, 0.00024961f, 0.00563625f, 0.0f, 0.0f},
  {0.00131840f, 0.00099430f, 0.00074960f, 0.00066005f, 0.00036626f, 0.00070192f, 0.00092548f, 0.00089301f, 0.00131038f, 0.00127857f, 0.00219713f, 0.00100817f, 0.00054105f, 0.00368739f, 0.00047608f, 0.00102648f, 0.00094759f, 0.00069226f, 0.00999315f, 0.0f},
  {0.00533241f, 0.00169359f, 0.00136609f, 0.00127915f, 0.00119152f, 0.00132844f, 0.00178697f, 0.00194579f, 0.00071553f, 0.01117956f, 0.00914460f, 0.00210897f, 0.00197461f, 0.00256159f, 0.00135781f, 0.00241601f, 0.00343452f, 0.00038538f, 0.00148001f, 0.02075171f}
};

// change the robinson frequencies to the probcons background frequencies
// calculate the log_q_blosum62[][] and log_q_blosum62_ratio[][] using the 
// probcons frequencies
void useProbconsFrequencies() {

	int i,j;
	char c, d;
	for(i=0;i<20;i++) {
		c = alphabetDefault[i];
		robinson_freq[am2num(c)] = emitSingleDefault[i]; 
	}
	if(debug>1) cout << "Frequencies: " << endl;
	if(debug>1) for(i=1;i<=20;i++) {
		cout << am[i] << "\t" << robinson_freq[i] << endl;
	}

	for(i=0;i<20;i++) {
	    for(j=0;j<=i;j++) {
		c = alphabetDefault[i];
		d = alphabetDefault[j];
		q_blosum62[am2num(c)][am2num(d)] = q_blosum62[am2num(d)][am2num(c)] = emitPairsDefault[i][j];
	    }
	}
	for(i=0;i<=20;i++) {
	    q_blosum62[i][0] = 0;
	    q_blosum62[0][i] = 0;
	}

	if(debug>1) cout << "Blosum62 matrix: " << endl;
	if(debug>1) for(i=1;i<=20;i++) {
	    cout << am[i] << " ";
	    for(j=1;j<=20;j++) {
		cout << q_blosum62[i][j] << " ";
	    }
	    cout << endl;
	}

	// get the logarithms
	if(debug>1) cout << "===============" << endl;
        get_log_robinson_freq();
        get_log_q_blosum62_ratio();
        get_log_q_blosum62();
	if(debug>1) cout << "===============" << endl;

}
