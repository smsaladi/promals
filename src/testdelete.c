#include "header_cpp.h"

int main() {

	char **aseq = new char * [20];

	char **bseq = new char * [2];

	delete [] aseq;

	aseq = bseq;

	cout << aseq << endl;

}
