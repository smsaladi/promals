#include "all.h"

int main(int argc, char **argv) {

        int i, j;
        int header = 0;
        int blocksize = 80;
        float gap_thr = 0.5;

        if(argc<3) {
                cout << "Usage: elimiategap alnfile -g gap_thr [-b blocksize] [-h header]" << endl;
                exit(0);
        }
        string argstr = "";
        for(i=1;i<argc;i++) {
                argstr = argv[i];
                if(argstr == "-b") blocksize = atoi(argv[i+1]);
                if(argstr == "-h") header = atoi(argv[i+1]);
                if(argstr == "-g") gap_thr = atof(argv[i+1]);
        }

        subalign *a = new subalign(argv[1]);

        int *mark = new int[a->alilen+1];
        int *mark1 = new int[a->nal+1];
        for(i=1;i<=a->nal;i++) mark1[i] = 1;
        
        int gap_count;
        for(i=1;i<=a->alilen;i++) {
                mark[i] = 1;
                gap_count = 0;
                for(j=1;j<=a->nal;j++) {
                        if(a->aseq[j-1][i-1]=='-') gap_count++;
                }
                if(gap_count*1.0/a->nal>=gap_thr) mark[i] = 0;
        }

        subalign *b = a->sub2align(mark1, mark);
        b->printali(blocksize, header);
}
