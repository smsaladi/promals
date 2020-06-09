#include  "all.h"

int main(int argc, char **argv) {

        int i, j;

        char header[500];

        strcpy(header, "CLUSTAL format multiple sequence alignment by PROMALS3D\n");
                       
        subalign myalign(argv[1]);
        
        cout << header << endl;

        int blocksize = atoi(argv[2]);
        if(blocksize <10) {
                blocksize = 70;
        }
        myalign.printali(blocksize, 0);


        return 0;
}

                        

                
