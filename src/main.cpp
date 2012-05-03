#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <cstdlib>

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

int main(int argc, char* argv[]){

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
