#include "all.h"

int main(int argc, char **argv) {

        int i;

        subalign *a = new subalign(argv[1]);
        int start = atoi(argv[2]);
        int end = atoi(argv[3]);
                

        int *mark = new int[a->alilen+1];
        int *mark1 = new int[a->nal+1];
        for(i=1;i<=a->nal;i++) {
                mark1[i] = 1;
                if(i>1000) mark1[i] = 0;
        }

        for(i=1;i<=a->alilen;i++) {
                if( (i>=start) && (i<=end) ) {
                        mark[i] = 1;
                }
                else mark[i] = 0;
        }
        /*
        for(i=1;i<=a->nal;i++) {
                a->aname[i-1] = new char [300];
                sprintf(a->aname[i-1], "%d", i);
        }*/

        subalign *b = a->sub2align(mark1, mark);
        b->printali(80, 0);
}
