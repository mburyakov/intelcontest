#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <cstdlib>
//#include "/home/mburyakov/bin/intel/ipp/include/ipp.h"
#include "../include/cyclichash.h"

using namespace std;

pair<string,string> getSequence(ifstream& s){
	string sequence;
	string name;

	if(s.get()=='>')
		getline(s,name);

	while (!s.eof()){
		int c = s.get(); //get char from reference Sequence file
		if(c == '>'){
			s.unget();
			break;
		}
		else if (c == 'G' || c=='T' || c=='A' || c=='C'){
			sequence.push_back(c);
		}
		// ignore every other chars (newline, etc)
	}

	return make_pair(name,sequence);
}

int bruteforce(int argc, char* argv[]){

	// first command line argument : number of worker threads to use
	// not used in this program at this time

	// first command line argument : window size
	size_t minMatchLength = atol(argv[2]);

	// second command line argument : reference sequence file
	ifstream refSeqFile(argv[3],ios::in);
	pair<string,string> ref = getSequence(refSeqFile);
	string refSeqName = ref.first;
	string refSeq = ref.second;

	// following command line arguments : other files containing sequences
	// result is stored in an associative array ordered by title
	map<string,string> otherSequences;
	for(int i=4; i<argc; i++){ // iterate over command arguments
		ifstream seqFile(argv[i],ios::in);
		while(!seqFile.eof()){
			pair<string,string> other = getSequence(seqFile);
			otherSequences[other.first] = other.second;
		}
	}

	// compare other sequences to reference sequence
	// iterate over other sequences
	for(map<string,string>::iterator sequencesIter = otherSequences.begin(); sequencesIter!=otherSequences.end(); sequencesIter++){
		// output sequence name
		cout << sequencesIter->first << "\n";
		string otherSeq = sequencesIter->second;

		// L[i][j] will contain length of the longest substring
		// ending by positions i in refSeq and j in otherSeq
		size_t **L = new size_t*[refSeq.length()];
		for(size_t i=0; i<refSeq.length();++i)
			L[i] = new size_t[otherSeq.length()];

		// iteration over the characters of the reference sequence
		for(size_t i=0; i<refSeq.length();i++){
			// iteration over the characters of the sequence to compare
			for(size_t j=0; j<otherSeq.length();j++){
				// if the characters are the same,
				// increase the consecutive matching score from the previous cell
				if(refSeq[i]==otherSeq[j]){
					if(i==0 || j==0)
						L[i][j]=1;
					else
						L[i][j] = L[i-1][j-1] + 1;
				}
				// or reset the matching score to 0
				else
					L[i][j]=0;
			}
		}

		// output the matches for this sequence
		// length must be at least minMatchLength
		// and the longest possible.
		for(size_t i=0; i<refSeq.length();i++){
			for(size_t j=0; j<otherSeq.length();j++){

				if(L[i][j]>=minMatchLength) {
					//this match can be shifted on i and j
					if(i+1<refSeq.length() && j+1<otherSeq.length() && L[i][j]<=L[i+1][j+1])
						continue;
					//this match can be shifted on j
					if(i<refSeq.length() && j+1<otherSeq.length() && L[i][j]<=L[i][j+1])
						continue;
					//this match can be shifted on i
					if(i+1<refSeq.length() && j<otherSeq.length() && L[i][j]<=L[i+1][j])
						continue;
					cout << i-L[i][j]+2 << " " << i+1 << " " << j-L[i][j]+2 << " " << j+1 << "\n";

					// output the matching sequences for debugging :
					//cout << refSeq.substr(i-L[i][j]+1,L[i][j]) << "\n";
					//cout << otherSeq.substr(j-L[i][j]+1,L[i][j]) << "\n";
				}
			}
		}

		for(size_t i=0; i<refSeq.length();++i)
			delete[] L[i];
		delete[] L;
	}
	return 0;
}

#define hashKernel 2147483647


//reads A -> 00. C -> 01, T -> 10, G -> 11
//input:  16 ascii letters
//output: 32 bit
void readDataBlock(unsigned long letters[], unsigned int* dest) {
    *dest = 0;
    for (int i=0;i<2;i++) {
        unsigned long tmp  = ((letters[i] & 0x0606060606060606l) << 5);
        tmp = tmp | (tmp >> 10);
        tmp = tmp | (tmp >> 20);
        tmp = tmp & 0x000000ff000000ffl;
        tmp = tmp | (tmp >> 24);
        tmp = tmp & 0x000000000000ffffl;
        *dest = *dest | (tmp << (i*16));
    }
}

/*unsigned long calcSingleHash(unsigned long const data[], size_t len, size_t shift) {
    unsigned long accum = 0;
    for (size_t i=0;i<len/8;i++) {
        accum *= hashKernel;
        accum += (data[i]);
    }
    return accum;
} */

/*unsigned long calcSingleAlignedHash(unsigned long const data[], size_t len) {
    unsigned long accum = 0;
    for (size_t i=0;i<len/8;i++) {
        accum *= hashKernel;
        accum += (data[i]);
    }
    return accum;
} */

void testReadDataBlock() {
    unsigned int dest;
    char input[17] = "AAAAAAAAAAAAACTG";
    readDataBlock((unsigned long *) &input,&dest);
	cout << dest << endl;
}

/* unsigned long * initHashKernelPows(size_t hashSize) {
    unsigned long *ans = new unsigned long[hashSize/8];
    unsigned long kpow = 1;
    for (size_t i=0;i<hashSize/8;i++) {
        ans[i] = kpow;
        kpow *= hashKernel;
    }
    return ans;
}  */

int main(int argc, char* argv[]) {

    //return bruteforce(argc, argv);
    size_t numerOfThreads = atol(argv[1]);
	size_t minMatchCharLen = atol(argv[2]);

	ifstream refSeqFile(argv[3],ios::in);
	pair<string,string> ref = getSequence(refSeqFile);
	string refSeqName = ref.first;
	string refSeq = ref.second;

	// following command line arguments : other files containing sequences
	// result is stored in an associative array ordered by title
	map<string,string> otherSequences;
	for(int i=4; i<argc; i++){ // iterate over command arguments
		ifstream seqFile(argv[i],ios::in);
		while(!seqFile.eof()){
			pair<string,string> other = getSequence(seqFile);
			otherSequences[other.first] = other.second;
		}
	}

	//char const* refCharSeq = refSeq.c_str();
	//size_t refCharLen = refSeq.length();

	//TODO: remove following 3 lines
	char const* refCharSeq = "AAAAAAAAAAAAACTGAAAAAAAAAAAAACTGAAAAAAAAAAAAACTGAAAAAAAAAAAAACTGAAAAAAAAAAAAACTGAAAAAAAAAAAAACTGAAAAAAAAAAAAACTGAAAAAAAAAAAAACTG";
	size_t refCharLen = 128;
	minMatchCharLen = 32; //TODO: if it is less than 32 then use other algorithm
    size_t minMatchBinLen = minMatchCharLen/4;

    size_t refBinLen = (refCharLen+127)/128*32;
    char *refBinSeq = new char[refBinLen];
    for (size_t i=0; i+16<=refCharLen; i+=16) {
        readDataBlock((unsigned long  *)(refCharSeq+i),(unsigned  *)(refBinSeq+i/4));
        cout << *((unsigned int *) (refBinSeq+i/4)) << endl;
    }

    size_t hashLen = minMatchBinLen/8*4;
    cout << "hashlen = " << hashLen << endl;

    CyclicHash hasher(hashLen, hashLen);

    unsigned int h0;
    h0 = hasher.singleHash((unsigned short *)refBinSeq);
    cout << "hash[0] = " << h0 << endl;
    cout << "hash[4] = " << hasher.singleHash((unsigned short *)(refBinSeq+4)) << endl;
    hasher.moveRight(h0,(unsigned short *)refBinSeq);
    cout << "hash[0+4] = " << h0 << endl;

    /*unsigned long *hashKernelPows = initHashKernelPows(hashSize);

	for (int i=0;i+hashSize<=refBinLen;i+=hashSize) {
	    cout << "i = " << (unsigned long const *)(refBinSeq+i) << endl;
        cout << calcSingleAlignedHash((unsigned long const *)(refBinSeq+i),hashSize) << endl;
	} */

    //testReadDataBlock();
    cout << endl << endl;
	return 0;
}
