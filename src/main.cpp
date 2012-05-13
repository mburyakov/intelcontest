#define DEBUG

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <cstdlib>
#include <omp.h>
#include <ipp.h>
#include <cmath>
//#include <sparsehash/dense_hash_map>
#include "../include/cyclichash.h"
#include "../include/reader.h"

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
        } else if (c == 'G' || c=='T' || c=='A' || c=='C'){
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


void testHash(CyclicHash hasher, char *refBinSeq) {
    unsigned int h0;
    cout << "testHash" << endl;
    h0 = hasher.singleHash((unsigned short *)refBinSeq);
    cout << "hash[0] = " << h0 << endl;
    cout << "hash[4] = " << hasher.singleHash((unsigned short *)(refBinSeq+4)) << endl;
    hasher.moveRight(&h0,(unsigned short *)refBinSeq);
    cout << "hash[0+4] = " << h0 << endl;
    h0 = hasher.singleHash((unsigned short *)refBinSeq,2);
    cout << "truesrc = " << *((unsigned int*) refBinSeq)<< endl;
    cout << "hash[0:2] = " << h0 << endl;

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

inline unsigned long index(unsigned long byte, unsigned char shift) {
    return (byte<<8)+shift;
}

inline void insert(map<unsigned int,unsigned long> &map,unsigned int hash, unsigned long byte, unsigned char shift) {
    map.insert(pair<unsigned int, unsigned long>(hash,index(byte,shift)));
}


int main(int argc, char* argv[]) {

    //return bruteforce(argc, argv);
    size_t numerOfThreads = atol(argv[1]);
    size_t numberOfWorkingThreads = numerOfThreads - 1;

    omp_set_dynamic(0);      // запретить библиотеке openmp менять число потоков во время исполнения
    omp_set_num_threads(numberOfWorkingThreads); // установить число потоков
    ippStaticInit();


    size_t minMatchCharLen = atol(argv[2]);

    //ifstream refSeqFile(argv[3],ios::in);
    Reader reader;
    Reader::mapped refFile = Reader::mapfile(argv[3]);
    
    
    pair<string,char *> namepair = reader.readName(refFile.pnt);
    
    size_t refBinLen = ((int)(refFile.pnt-namepair.second+refFile.len))/128*32+32;
    char *refBinSeq = (char *) ippsMalloc_16u(refBinLen);
    Reader::mapped refcharpair = reader.readData((unsigned char *) namepair.second,(unsigned char *)refBinSeq,((int)(refFile.pnt-namepair.second+refFile.len)));
    long refCharLen = refcharpair.len;
    assert(refBinLen >= (refCharLen)/128*32+32);
    
    
    //pair<string,string> ref = getSequence(refSeqFile);
    //string refSeqName = ref.first;
    //string refSeq = ref.second;

    //TODO: read data with mmap
    //char const* refCharSeq = refSeq.c_str();
    //size_t refCharLen = refSeq.length();
    //TODO: remove following 3 lines
    //char const* refCharSeq = "CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTAAAAACCTGAAAAAAAAAAAAACTGAAAAAAAAAAAAACTGAAAAAAAAAAAAACTGAAAAAAAAAAAAACTGAAAAAAAAAAAAACTG";
    //size_t refCharLen = 128;
    minMatchCharLen = 32; //TODO: if it is less than 32 then use other algorithm

    size_t minMatchBinLen = minMatchCharLen/4;
    size_t refBitLen = refCharLen*2;
    //char *refBinSeq = new char[refBinLen];//TODO: ippMalloc
    //for (size_t i=0; i+16<=refCharLen; i+=16) {
    //    readDataBlock((unsigned long  *)(refCharSeq+i),(unsigned  *)(refBinSeq+i/4));
    //    //cout << *((unsigned int *) (refBinSeq+i/4)) << endl;
    //}

    size_t hashLen = minMatchBinLen/8*4;
    cout << "hashlen = " << hashLen << endl;

    CyclicHash hasher(hashLen, hashLen);

    testHash(hasher,refBinSeq);
    
    map<unsigned int,unsigned long> stripes;

    for (int shift=0; shift<hasher.charLen*8; shift+=2) {
        int n = (refBitLen-shift)/hasher.charLen-hasher.wordLen;
        if (n>=0) {
            unsigned int hash = hasher.singleHash((basechartype*)refBinSeq,shift);
            insert(stripes,hash,0,shift);
        
            for (int i = 0; i<n; i+=hasher.charLen) {
                hasher.moveRight(&hash,(basechartype*) (refBinSeq+i), shift);
                insert(stripes,hash,i,shift);
            }
        }
    }

    // following command line arguments : other files containing sequences
    // result is stored in an associative array ordered by title
    map<string,string> otherSequences;
    long threadLen = 0;
    for (int i=4; i<argc; i++) {
        threadLen += Reader::lookLen(argv[i]);
    }
    cout << "threadLen = " << threadLen << endl;
    threadLen /= numberOfWorkingThreads;
    threadLen++;
    
    cout << "threadLen = " << threadLen << endl;
    
    vector<Reader::mapped> *files = new vector<Reader::mapped>[numberOfWorkingThreads];
    vector<Reader::mapped> *bins = new vector<Reader::mapped>[numberOfWorkingThreads];
    vector<string> *names = new vector<string>[numberOfWorkingThreads];
    string *firstnames = new string[numberOfWorkingThreads];
    long *whereToStart = new long[numberOfWorkingThreads+10];
    long *whereToStop = new long[numberOfWorkingThreads+10];
    {
        int j=0;
        long thisThreadLen = 0;
        whereToStart[0] = 0;
        Reader::mapped refFile = Reader::mapfile(argv[4]);
        for (int i=4; i<argc; i++) {
            files[j].push_back(refFile);
            if (thisThreadLen+refFile.len>=threadLen) {
                i--;                
                whereToStop[j] = threadLen - thisThreadLen;
                whereToStart[j+1] = whereToStop[j];
                thisThreadLen = -whereToStop[j];
                
                cout << "j = " << j << endl;
                cout << "i = " << i << endl;
                cout << "thisThreadLen = " << thisThreadLen << endl;
                cout << " whereToStop[j] = " <<  whereToStop[j] << endl;
                cout << " whereToStart[j+1] = " <<  whereToStart[j+1] << endl;
                j++;
            } else {
                thisThreadLen += refFile.len;
                if (i+1<argc) {
                    refFile = Reader::mapfile(argv[i+1]);
                }
                cout << "i = " << i << endl;
            }
        }
        whereToStop[j] = refFile.len;
        cout << " whereToStop[j] = " <<  whereToStop[j] << endl;
    }
    
    //#pragma omp parallel for
    for (int thread = 0; thread<numberOfWorkingThreads; thread++) {
        for (int file = 0; file<files[thread].size(); file++) {
            long start = file==0 ? whereToStart[thread] : 0;
            long stop = file==files[thread].size()-1 ? whereToStop[thread] : files[thread][file].len;
            stop = min(stop+hashLen,(unsigned long)refFile.len);      
            cout << "thread = " << thread << endl;
            cout << "start = " << start << endl;
            cout << "stop = " << stop << endl;
            Reader::mapped refFile = files[thread][file];
            size_t refBinLen = (stop-start)/32*32+32;
            char *refBinSeq = (char *) ippsMalloc_16u(refBinLen);
            Reader::mapped bin;
            Reader::mapped bin1;
            bin.pnt = refBinSeq;
            while (stop>start) {
                cout << "start = " << start << endl;
                bin1 = reader.readData((unsigned char *) refFile.pnt+start,(unsigned char *)refBinSeq,stop-start);                
                bin.len = bin1.len;
                //cout << " text = " << refFile.pnt+start << endl;
                //cout << "bin.pnt = " << bin.pnt << endl;
                bins[thread].push_back(bin);
                refBinSeq += bin1.len/4/32*32;
                start = (bin1.pnt-refFile.pnt);
                if (stop>start) {
                    pair<string,char *> namepair = reader.readName(bin1.pnt);
                    names[thread].push_back(namepair.first);
                    firstnames[thread+1] = namepair.first;
                    start = namepair.second-refFile.pnt;
                    cout << "strt = " << start << endl;
                }
            }         
                      
            
            
            //long refCharLen = reader.readData((unsigned char *) namepair.second,(unsigned char *)refBinSeq,((int)(refFile.pnt-namepair.second+refFile.len)));
            //assert(refBinLen >= (refCharLen)/128*32+32);
            //for (long len = start; len < stop; len++){
            //    putchar (files[thread][file].pnt[len]);   
            //}
        }
    }
    
    for (int thread = 0; thread<numberOfWorkingThreads; thread++) {
        cout << "thread = " << thread << endl;
        for (int bin = 0; bin<names[thread].size(); bin++) {
            cout << names[thread][bin] << endl;
        }
    }

    /*int *whichSeqStart = new int[numberOfWorkingThreads];
    int *whichSeqStop = new int[numberOfWorkingThreads];
    int *whichPosStart = new int[numberOfWorkingThreads];
    int *whichPosStop = new int[numberOfWorkingThreads];

    long threadLen = otherSequencesSumLen/numberOfWorkingThreads;
    for ((long i=0),(int j=0); i<otherSequencesNum && j<numberOfWorkingThreads; j++) {
        long thisThreadLen = 0;
        whichSeqStart[j] = i;
        if (j>0) {
            whichPosStart[j] = whichPosStop[j-1];
        } else {
            whichPosStart[j] = 0;
        }
        while (i<otherSequenceNum && otherSequencesLen[i]+2*thisThreadLen-2*threadLen<0) {
             thisThreadLen += otherSequencesLen[i];
             i++;
        }
        whichSeqStop[j] = i;
        whichPosStop[j] = threadLen - thisThreadLen;
    }
*/
    /*unsigned long *hashKernelPows = initHashKernelPows(hashSize);

    for (int i=0;i+hashSize<=refBinLen;i+=hashSize) {
        cout << "i = " << (unsigned long const *)(refBinSeq+i) << endl;
        cout << calcSingleAlignedHash((unsigned long const *)(refBinSeq+i),hashSize) << endl;
    } */

    //testReadDataBlock();
    cout << endl << endl;
    return 0;
}
