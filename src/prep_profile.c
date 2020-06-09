#include "all.h"

// -q: queryfile containing the fasta record of the query
// -o: blast output file im m6 format
// -n: query name

int main(int argc, char **argv) {

        int i, j, k;

        string argstr;
        char queryfile[500];
        char blastoutfile[500];
        char queryname[500];
        for(i=1;i<argc;i++) {

                argstr = argv[i];
                if(argstr == "-q") {
                        strcpy(queryfile, argv[i+1]);
                        continue;
                }
                if(argstr == "-o") {
                        strcpy(blastoutfile, argv[i+1]);
                        continue;
                }
                if(argstr == "-n") {
                        strcpy(queryname, argv[i+1]);
                        continue;
                }
                //if(argstr == "-j") { iteration_number = atoi(argv[i+1]); continue; }
        }
                
        ifstream fp(queryfile, ios::in);
        char tmpline[10000];
        char queryseq[10000];
        queryseq[0] = '\0';
        while(fp.getline(tmpline, 10000)){
                if ( !fp.good() ) break;
                fp.getline(tmpline, 10000);
                if(tmpline[0]=='>') continue;
                strcat(queryseq, tmpline);
        }
        fp.close();

        //cout << queryseq << endl;

        subalign *x1 = read_blastpgp_result(blastoutfile, queryname, queryseq);
	subalign *x = x1->purge_align(low_id_thr, 0.9, max_num_sequences, 0.9);
        x->printali( (string(queryname)+string(".aln") ).c_str(), 80); 
        
}
        
        
