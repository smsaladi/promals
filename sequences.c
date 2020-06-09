#include "sequences.h"

sequences::sequences() {
	nseqs = 0;
}

void sequences::get_map() {

	int i;
	
	for(i=1;i<=nseqs;i++) {
		name2seq.insert(pair<string, string>(name[i], seq[i]));
		name2index.insert(pair<string,int>(name[i], i) );
	}
	
	//for(i=1;i<=nseqs;i++) { cout << i << " " << name2index.find(name[i])->first << " " << name2index.find(name[i])->second << endl; }

}

// copy constructor
sequences::sequences(const sequences &s) {
	
	int i;
	string str;
	
	nseqs = s.nseqs;
	isAlign = s.isAlign;
	name.erase(name.begin(), name.end());
	seq.erase(seq.begin(), seq.end() );
	// unlike subalign, index begins from 1 and not 0
	name.push_back("");
	seq.push_back("");
	for(i=1;i<=nseqs;i++) {
		str = s.name[i];
		name.push_back(str);
		//name.push_back(s.name[i]);
		str = s.seq[i];
		seq.push_back(str);
		//seq.push_back(s.seq[i]);
	}

	distMat = 0;
}

// copy constructor, input is a subalign
sequences::sequences(const subalign &s) {
	
	int i;
	
	nseqs = s.nal;
	isAlign = 1;
	name.erase(name.begin(), name.end());
	seq.erase(seq.begin(), seq.end() );
	// unlike subalign, index begins from 1 and not 0
	name.push_back("");
	seq.push_back("");
	for(i=0;i<s.nal;i++) {
		name.push_back(s.aname[i]);
		seq.push_back(s.aseq[i]);
	}

	distMat = 0;
}

// read input from a fasta format file; remove non-alphabetic letters
// read input from a fasta format file; remove non-alphabetic letters
void sequences::readFasta(char *inputFastaName, int zap) {

	int i,j,k;
	string tmpstring;
	string tmpstring1;
	char tmp1[10000];
	
	nseqs = 0;
	isAlign = 0;
	name.push_back("");
	seq.push_back("");

	ifstream fastaFile(inputFastaName);
	if(!fastaFile) {
		cout << endl;
		cout << "Error: failed to read the fasta file: " << inputFastaName << endl;
		exit(1);
	}

	while(fastaFile) {
		//fastaFile >> tmpstring;
		fastaFile.getline(tmp1, 10000);
		tmpstring = tmp1;
		//cout << tmpstring << "===" << endl;
		if(tmpstring[0]=='>') {
			nseqs++;
			//cout << nseqs << endl;
			tmpstring = tmpstring.substr(1); 
			tmpstring1 = tmpstring;
			for(i=0;i<tmpstring.length();i++) {
			    if( (tmpstring[i]!=' ') && (tmpstring[i]!='\t') &&(tmpstring[i]!='\n') ) {
				break;
			    }
			    tmpstring1.erase(0,1);
			}
			if(tmpstring1.length()==0) {
				cout << endl;
			    cout << "Error: one or more sequence names of input fasta are empty" << endl;
			    exit(1);
			}
			tmpstring = tmpstring1;
			tmpstring1  = "";
			for(i=0;i<tmpstring.length();i++) {
				if( (tmpstring[i]==' ') || (tmpstring[i]=='\t') || (tmpstring[i]=='\n') ) {
					break;
				}
				tmpstring1 += tmpstring[i];
			}
			tmpstring = tmpstring1;
			string::size_type loc = tmpstring.find(" ", 0);
			if( loc != string::npos ) {
			    tmpstring = tmpstring.substr(0, loc);
			}
			name.push_back(tmpstring);
			seq.push_back("");
			//cout << tmpstring << "===" << endl;
		}
		else {
			seq[nseqs] += tmpstring;
		}
		tmpstring = "";
	}	

	if(nseqs==0) {
		cout << endl;
		cout << "Error: wrong input format - should be in fasta format" << endl;
		exit(1);
	}
	if(nseqs==1) {
		cout << endl;
		cout << "Error: only one sequence is read" << endl;
		cout << "       - there should be at least two sequences for alignment" << endl;
		exit(1);
	}

	for(i=1;i<=nseqs;i++) {
		if(zap) zapLetters(seq[i]);
		toUpper(seq[i]);
		//cout << name[i] << endl;
		//cout << seq[i] << endl;
	}

	for(i=1;i<=nseqs;i++) {
		if(seq[i].length()<=0) {
			cout << endl;
			cout << "Error: sequence length is 0 for " <<  name[i] << endl;
			exit(1);
		}
	}

}

// zapLetters: remove letters that are not alphabets such as gaps
void sequences::zapLetters(string &str) {

	int i;
	int num=0;
	string tmpstr;
	for(i=0;i<str.length();i++) {
		if( (str[i] >= 'A') && (str[i] <='Z') ) {
			tmpstr += str[i];
			continue;
		}
		if( (str[i] >= 'a') && (str[i] <='z') ) {
			tmpstr += str[i];
		}
	}
	//cout << tmpstr << endl;
	str = tmpstr;
}

// toUpper: change the letters to upper case
void sequences::toUpper(string &str) {
	
	int i;
	for(i=0;i<str.length();i++) {
		str[i] = toupper(str[i]);
	}
}

// print out the sequences
void sequences::printSeqs() {
	int i;
	cout << "number of sequences: "<< nseqs << endl;
	for(i=1;i<=nseqs;i++) {
		// zapLetters(seq[i]);
		cout << ">" << name[i] << endl;
		cout << seq[i] << endl;
	}
}

// convert to Dayhoff 6-alphabet
void sequences::toDayhoff6() {
	
	int i,j;
	for(i=1;i<=nseqs;i++) {
		for(j=0;j<seq[i].length();j++) {
			seq[i][j] = am2d6(seq[i][j]);
		}
	}
}
	
// derive a K-mer map(table, hash) for a sequence(str)
dayhoff6Table sequences::seqToD6t(string str, int K) {
	int i,j;

	dayhoff6Table tmpd6t;

	if (str.length() < K) return tmpd6t;
	
	for(i=0;i<str.length()-K+1;i++) {
		tmpd6t[str.substr(i,K)]++;
	}

	return tmpd6t;
}

// generate the vector of K-mer maps
void sequences::generateD6t(int K) {
	int i;
	dayhoff6Table d;

	d6t.erase(d6t.begin(), d6t.end());
	d6t.push_back(d);
	for(i=1;i<=nseqs;i++) {
		d6t.push_back( seqToD6t(seq[i], K) );
	}
}

// the difference between two dayhoff6tables
int sequences::diffCountD2t(dayhoff6Table t1, dayhoff6Table t2) {
        int i,j;
	int count = 0;

	dayhoff6Table::iterator iter1 = t1.begin();

	while(iter1 !=t1.end() ) {
		//cout << iter1->first << "\t" << iter1->second << endl;
		count += abs(iter1->second - t2[iter1->first] );
		iter1++;
	}

	dayhoff6Table::iterator iter2 = t2.begin();
	while(iter2 != t2.end() ) {
		if(t1[iter2->first] == 0) { count += iter2->second; }
		iter2++;
	}
	//cout << "Count: " << count << endl;
	return count;
}


// the common counts between two dayhoff6tables
int sequences::commonCountD2t(dayhoff6Table t1, dayhoff6Table t2) {
        int i,j;
	int count = 0;

	dayhoff6Table::iterator iter1 = t1.begin();

	while(iter1 !=t1.end() ) {
		//cout << iter1->first << "\t" << iter1->second << endl;
		if(iter1->second> t2[iter1->first]) {
			count += t2[iter1->first];
		} else count += iter1->second;
		iter1++;
	}

	//dayhoff6Table::iterator iter2 = t2.begin();
	//while(iter2 != t2.end() ) {
	//	if(t1[iter2->first] == 0) { count += iter2->second; }
	//	iter2++;
	//}
	//cout << "Count: " << count << endl;
	return count;
}


void sequences::d6t2DistMat(int K) {

	int i,j;
	int minLength;

	if(!distMat) {
		distMat = dmatrix(nseqs, nseqs);
	}
	
	for(i=1;i<=nseqs;i++) {
		for(j=1;j<=nseqs;j++) {
			if( (seq[i].length()<K) || (seq[j].length()<K) ) {
				distMat[i][j] = 1;
				if(i==j) distMat[i][j] = 0;
				continue;
			}	
			if(seq[i].length() > seq[j].length() ) {
				minLength = seq[j].length();	
			} else minLength = seq[i].length();
			// cout << i << "\t" << j << "\t" << minLength << endl;
			// distMat[i][j] = 1; // - 1.0 * commonCountD2t(d6t[i],d6t[j])/(minLength-K+1);
			//distMat[i][j] = 1 - 1.0 * commonCountD2t(d6t[i],d6t[j])/(minLength-K+1);
			distMat[i][j] = log(0.1 + 1.0 * commonCountD2t(d6t[i],d6t[j])/(minLength-K+1) )/log(10.0);
			if(distMat[i][j]>=0) distMat[i][j] = 0;
			else distMat[i][j] = 0-distMat[i][j];
		}
	}
}

void sequences::printDistMat() {

	int i,j;
	
	cout << "  " << nseqs << endl;
	for(i=1;i<=nseqs;i++) {
		cout << left << setw(10)  << name[i].substr(0,9) << " ";
		for(j=1;j<=nseqs;j++) {
			cout << left << setw(10) << setprecision(6) << distMat[i][j];
		}
		cout << endl;
	}
}

void sequences::seqIdentity2DistMat() { // dist = 1 - (percentage seq. ident.)

	int i,j;
	
	if(!isAlign) {
		cout <<"warning:sequences may not be aligned"<<endl;
	}
	if(!distMat) {
		distMat = dmatrix(nseqs, nseqs);
	}
	for(i=1;i<=nseqs;i++){
		//cout << seq[i].length() << endl;
		distMat[i][i] = 0;
		for(j=i+1;j<=nseqs;j++) {
			distMat[i][j] = distMat[j][i] = 1.0;
			distMat[i][j] = distMat[j][i] = 1-seqIdentity(seq[i], seq[j], seq[i].length());

		}
	}
}

double seqIdentity(const string &c1, const string &c2, int length) {

	int i,j;
	int ip = 0;
	int len1=0, len2=0;
	for(i=0;i<length;i++) {
		if( (c1[i]==c2[i])&& isalpha(c1[i])  && isalpha(c2[i])  ) {
			ip++;
		}
		if( isalpha(c1[i]) ) len1++;
		if( isalpha(c2[i]) ) len2++;
	}
	return (1.0*ip/(len1?len2:(len1<len2) ) );
}

