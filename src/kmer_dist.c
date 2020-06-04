#include "kmer_dist.h"


char am2d6(char x) {
	switch (x) {
		case 'A':
		case 'G':
		case 'P':
		case 'S':
		case 'T':
		case 'a':
		case 'g':
		case 'p':
		case 's':
		case 't': return 'A';
		case 'c':
		case 'C': return 'B';
		case 'D':
		case 'E':
		case 'N':
		case 'Q':
		case 'd':
		case 'e':
		case 'n':
		case 'q': return 'C';
		case 'F':
		case 'W':
		case 'Y':
		case 'f':
		case 'y':
		case 'w': return 'D';
		case 'H':
		case 'K':
		case 'R':
		case 'h':
		case 'r':
		case 'k': return 'E';
		case 'I':
		case 'L':
		case 'V':
		case 'M':
		case 'i':
		case 'l':
		case 'v':
		case 'm': return 'F';
		// if not the standard 20 amino acids
		default: return 'A';
	}
}
		
		

