#ifndef __sequences_
#define __sequences_

#include "all.h"
#include <string>
#include "header_cpp.h"
#include <map>
#include <vector>


class sequences {
	
     public: 
	vector <string> name;
	vector <string> seq;

	map<string, string> name2seq;
	void get_map();

	map<string, int> name2index;
	
	int nseqs;

	int isAlign;

	//constructors
	sequences(); // {nseqs = 0; distMat = 0;}
	sequences(char *inputFastaName, int zap) {
		readFasta(inputFastaName, zap);
		distMat = 0;
	}
	// copy constructors
	sequences(const sequences &s);
	sequences(const subalign &s);

	// reading and processing the sequences
	void readFasta(char *inputFastaName, int zap);
	void zapLetters(string &str);
	void toUpper(string &str);
	void toDayhoff6();
	void printSeqs();

	// generate and process maps (hashes) of K-mers
	vector <dayhoff6Table> d6t;
	dayhoff6Table seqToD6t(string str, int K);
	void generateD6t(int K);
	int diffCountD2t(dayhoff6Table t1, dayhoff6Table t2);
	int commonCountD2t(dayhoff6Table t1, dayhoff6Table t2);

	// distance matrix
	double **distMat;
	void d6t2DistMat(int K); //use dayhoff6 table to get distMat
	void printDistMat();
	void seqIdentity2DistMat(); // dist = 1 - (percentage seq. ident.)

};

double seqIdentity(const string &c1, const string &c2, int length);

#endif
