#include "param.h"

double ave_grp_thr=1;
double minProb = 0.01;
string outFile = "";
int probconsBLOSUM = 1; // probcons blosum62
int useLocal = 3;
float weightG = 1;
int ss = 1;
int solv = 1;
int unaligned = 0;
char parameter_file[200];
char parameter_file1[200];
char parameter_file2[200];
int relax_number = 2;
int reverse_align_order = 0;
double id_thr = 0.6;


void getParameter(int argc, char **argv, int prog) {

	int i,j;
	string argStr;
	
	//cout << argc << endl;
	parameter_file[0] = '\0';

	if(argc<=1) {
	    printHelp(prog);
	    exit(0);
	}    

	if(argv[1][0] == '-') {
	    printHelp(prog);
	    exit(0);
	}

	for(i=1;i<argc;i++) {
		argStr = argv[i];
		//if(argStr == "-a") { ave_grp_thr = atof(argv[i+1]); continue; }
		if(argStr == "-minprob") {
			minProb = atof(argv[i+1]);
			continue;
		}
		if(argStr == "-outfile") {
			outFile = argv[i+1];
			continue;
		}
		if(argStr == "-pb") {
			probconsBLOSUM = atoi(argv[i+1]);
			continue;
		}
		if(argStr == "-local") {
			useLocal = atoi(argv[i+1]);
		}
		if(argStr == "-wg") {
			weightG = atof(argv[i+1]);
			assert( weightG>=0 );
			assert( weightG<=1 );
		}
		if(argStr == "-ss") {
			ss = atoi(argv[i+1]);
			try {
				if( (ss!=1) && (ss!=3) ) throw 3;
			}
			catch(int ae) {
				cout << endl;
				cout << "Error: with option -ss" << endl;
				cout << "       number of secondary structure types must be 1 or 3" << endl;
				exit(1);
			}
		}
		if(argStr == "-solv") {
			solv = atoi(argv[i+1]);
			try {
				if( (solv<=0) || (solv>3) ) throw 3;
			}
			catch (int ae) {
				cout << endl;
				cout << "Error: with option -solv" << endl;
				cout << "       number of solvent accessibility categories must be 1,2 or 3" << endl;
				exit(1);
			}
		}
		if(argStr == "-unaligned") {
			unaligned = atoi(argv[i+1]);
			try {
				if( (unaligned<0) || (unaligned>1) ) throw 3;
			}
			catch (int ae) {
				cout << endl;
				cout << "Error: with option -unaligned" << endl;
				cout << "       number of unaligned match state must be 0 or 1" << endl;
				exit(1);
			}
		}
		if(argStr == "-param") {
			strcpy(parameter_file, argv[i+1]);
		}
		if(argStr == "-param1") {
			strcpy(parameter_file1, argv[i+1]);
		}
		if(argStr == "-param2") {
			strcpy(parameter_file2, argv[i+1]);
		}
		if(argStr == "-r") {
		        relax_number = atoi(argv[i+1]);
		}
		if(argStr == "-reverse") {
			reverse_align_order = atoi(argv[i+1]);
		}
		if(argStr == "-id_thr") {
			id_thr = atof(argv[i+1]);
		}
	}
	//cout << "Here" << endl;
}

void printParameters() {

	cout << endl << "MUMMALS - MUltiple alignment with Multiple MAtch state models of Local Structure" << endl;
	cout << endl << "List of parameters: " << endl;
	//cout << "  ave_grp_thr: " << ave_grp_thr << endl;
	//cout << "  minprob: " << minProb << endl;
	//cout << "  outfile name: " << outFile << endl;
	//cout << "  using probscons blosum62: " << probconsBLOSUM << endl;
	//cout << "  using local alignment option: " << useLocal;
	/*switch (useLocal) {
		case 1: cout << "  local only" << endl; break;
		case 0: cout << "  global only" << endl; break;
		case 2: cout << "  both local and global" << endl; 
			cout << "  weight of global: " << weightG << endl; break;
		default: cout << "   global only" << endl;
	}*/
	cout << "  Number of secondary structure types (ss): " << ss << endl;
	cout << "  Number of solvent accessibility categories (solv): " << solv << endl; 
	cout << "  Number of separate match state for unaligned regions (unaligned): " << unaligned << endl;
	cout << "  Identity threshold (id_thr): " << id_thr << endl;
	cout << "  Parameter file (param): " << parameter_file << endl;
	//cout << "  relax_number: " << relax_number << endl;
	//cout << "  reverse align order: " << reverse_align_order << endl;
}

void printHelp(int prog) {

    
	if(prog == 1)  // mummals
	{

	cout << endl << " MUMMALS - MUltiple alignment with Multiple MAtch state models of Local Structure" << endl << endl;

	cout << " Usage:" << endl;
	cout << " \t mummals input_fasta [options]" << endl;
	cout << endl;

	cout << " options: " << endl << endl;

	cout << " -ss         Number of secondary structural types, 1 or 3" << endl;
	cout << " -solv       Number of solvent accessibility categories, 1, 2 or 3" << endl;
	cout << " -unaligned  Number of additional match states for unaligned regions, 0 or 1" << endl;
	cout << " -param      Input parameter file for hidden Markov model" << endl;
	cout << " -id_thr     Indentity threshold above which neighboring groups are aligned in a fast way, between 0 and 1, default: 0.6" << endl;
	cout << " -outfile    Output file name" << endl;
	cout << endl;

	cout << " Example: using model HMM_1_3_1" << endl << endl;
	cout << " mummals tmp1.fa -ss 3 -solv 1 -unaligned 1 -param hmm_parameters/dataset_0.20_0.40_0.60_abcd.dali.solv1_ss3.mat -outfile tmp1.mummals.aln" << endl;

	cout << endl;
	}

	else if (prog == 0) // meta_align
	{
	    cout << endl << "meta_align - merge several multiple alignments into one alignment based on consistency";
	    cout << endl << endl;

	    cout << " Input should be a file that contains a list of alignment file names." << endl;
	    cout << " Each alignment file should be in fasta format." << endl;
	    cout << endl;
	    cout << " Usage:" << endl;
	    cout << " \t meta_align input_file_list [options]" << endl;
	    cout << endl;

	    cout << " options: " << endl << endl;
	    cout << " -outfile    Output file name " << endl;
	    cout << "	   if not specified, it will be input_file_list.meta_aln.aln" << endl;

	    cout << endl;
	}
}
		
