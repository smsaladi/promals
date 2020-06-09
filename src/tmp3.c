#include "header_cpp.h"

int main() {

        ifstream fp("test.fa", ios::in);

        char tmpstr[1500];
        
        while(fp){
               fp.getline(tmpstr, 3000000);
               if(fp) cout << tmpstr << endl;
        }

        fp.close();
        fp.clear();

        fp.open("test.fa", ios::in);
        while(fp){
                fp.getline(tmpstr, 3000000);
                if(fp) cout << tmpstr << endl;
        }


        exit(0);

        long int a1 = time(NULL);
        cout << getpid() << endl;
        
        srand(a1);
        int myrand = rand();

        char tmpstring[100];
        sprintf(tmpstring, "%d_%d", getpid(), myrand);
        cout << tmpstring << endl;
        

        cout << myrand << endl;
        cout << a1 << endl;
}
