#include "all.h"

int main(int argc, char **argv) {

        int i, j;

        subalign *a = new subalign(argv[1]);

        int gap_count;

        for(i=1;i<=a->alilen;i++) {
                gap_count=0;
                for(j=1;j<=a->nal;j++) {
                        if(a->alignment[j][i]==0)gap_count++;
                }
                cout << i << " " << 1.0 * gap_count/a->nal << endl;
        }
}
