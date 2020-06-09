#include "header_cpp.h"

int main() {

        cout << sizeof(int) << endl;

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
