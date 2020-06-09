#ifndef __kmer_dist_
#define __kmer_dist_

#include <map>
#include <string>
#include "amino.h"
using namespace std;

// map a k-mer string to its number of occurances
typedef std::map<string, int, std::less<string> > dayhoff6Table;

char am2d6(char x);

#endif
