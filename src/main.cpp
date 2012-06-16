#define DEBUG

#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
//#include <hash_map>
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


void testHash(CyclicHash &hasher, char *refBinSeq) {
    unsigned int h0;
    cout << "testHash" << endl;
    h0 = hasher.singleHash((unsigned short *)refBinSeq);
    cout << "hash[0] = " << h0 << endl;
    cout << "hash[4] = " << hasher.singleHash((unsigned short *)(refBinSeq+4)) << endl;
    hasher.moveRight(&h0,(unsigned short *)refBinSeq);
    cout << "hash[0+4] = " << h0 << endl;
    h0 = hasher.singleHash((unsigned short *)refBinSeq,2);
    //cout << "truesrc = " << *((unsigned int*) refBinSeq)<< endl;
    h0 = hasher.singleHash((unsigned short *)refBinSeq,2);
    cout << "hash[0:2] = " << h0 << endl;
    h0 = hasher.singleHash((unsigned short *)refBinSeq,4);
    cout << "hash[0:4] = " << h0 << endl;
    h0 = hasher.singleHash((unsigned short *)refBinSeq,6);
    cout << "hash[0:6] = " << h0 << endl;
    h0 = hasher.singleHash((unsigned short *)refBinSeq,8);
    cout << "hash[0:8] = " << h0 << endl;

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

inline void insert(multimap<unsigned int,unsigned long> &map,unsigned int hash, unsigned long byte, unsigned char shift) {
    map.insert(pair<unsigned int, unsigned long>(hash,index(byte,shift)));
}

struct foundation {
    long seqNum;
    long seqStartPos;
    long seqEndPos;
    long refStartPos;
};

struct anspos {
    unsigned short seqNum;
    long refStartPos;
    long seqStartPos;
};

union ansunion {
    anspos ansstruct;
    long something;
};

bool compareans(foundation const &a, foundation const &b) {
    if (a.seqNum < b.seqNum) {
        //cout << "<";
        return true;
    }
    if (a.seqNum > b.seqNum) {
        //cout << "<";
        return false;
    }
    if (a.refStartPos+a.seqEndPos-a.seqStartPos < b.refStartPos+b.seqEndPos-b.seqStartPos) {
        //cout << "<";
        return true;
    }
    if (a.refStartPos+a.seqEndPos-a.seqStartPos > b.refStartPos+b.seqEndPos-b.seqStartPos) {
        //cout << "<";
        return false;
    }
    if (a.seqEndPos < b.seqEndPos) {
        return true;
    }
    return false;
}

bool compareprev(foundation const &a, foundation const &b) {
    if (a.seqNum < b.seqNum) {
        //cout << "<";
        return true;
    }
    if (a.seqNum > b.seqNum) {
        //cout << "<";
        return false;
    }
    if (a.refStartPos-a.seqStartPos < b.refStartPos-b.seqStartPos) {
        return true;
    }
    if (a.refStartPos-a.seqStartPos > b.refStartPos-b.seqStartPos) {
        return false;
    }
    if (a.refStartPos < b.refStartPos) {
        return true;
    }
    if (a.refStartPos > b.refStartPos) {
        return false;
    }
    if (a.seqStartPos < b.seqStartPos) {
        return true;
    }
    return false;
}

bool operator<(anspos const &a, anspos const &b) {
    if (a.seqNum < b.seqNum) {
        return true;
    }
    if (a.refStartPos < b.refStartPos) {
        return true;
    }
    if (a.seqStartPos < b.seqStartPos) {
        return true;
    }
    return false;
}

bool operator==(anspos const &a, anspos const &b) {
    if (a.seqNum != b.seqNum) {
        return false;
    }
    if (a.refStartPos != b.refStartPos) {
        return false;
    }
    if (a.seqStartPos != b.seqStartPos) {
        return false;
    }
    return true;
}

bool operator<=(anspos const &a, anspos const &b) {
    return (a<b) || (a==b);
}

bool operator>(anspos const &a, anspos const &b) {
    return !(a<=b);
}

bool operator>=(anspos const &a, anspos const &b) {
    return !(a<b);
}



/*inline bool tryCandidate(long seqPos, long refPos, int seqNum, char* pnt, char* refBinSeq, bool goBack, CyclicHash &hasher) {
    char *refpnt = refBinSeq + refPos>>8;
    if (goBack) {
        hasher.deShiftChar(pnt, basechartype dest[],size_t shift) {
    }
}*/


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
    Reader::mapped refcharpair = reader.readData((unsigned char *) namepair.second,(unsigned char *)refBinSeq,((int)(refFile.pnt-namepair.second+refFile.len))).first;
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
    //minMatchCharLen = 32; //TODO: if it is less than 32 then use other algorithm

    if (minMatchCharLen < 32) {
        return bruteforce(argc, argv);
    }
    return 0;
    
    size_t minMatchBinLen = minMatchCharLen/4;
    size_t refBitLen = refCharLen*2;
    //char *refBinSeq = new char[refBinLen];//TODO: ippMalloc
    //for (size_t i=0; i+16<=refCharLen; i+=16) {
    //    readDataBlock((unsigned long  *)(refCharSeq+i),(unsigned  *)(refBinSeq+i/4));
    //    //cout << *((unsigned int *) (refBinSeq+i/4)) << endl;
    //}

    size_t hashLen;
    cout << "hashlen = " << hashLen << endl;
    
    CyclicHash *hasher;
    
    if (minMatchCharLen < 31) {
        hasher = new TinyHash[numberOfWorkingThreads];
        hashLen = 2;
    } else {
        hashLen = minMatchBinLen/8*4;
        hasher = new CyclicHash[numberOfWorkingThreads](hashLen,hashLen);        
    }    
    
    cout << "hashlen = " << hashLen << endl;
    
    
    //testHash(hasher[0],refBinSeq);
    //return 0;
    
    multimap<unsigned int,unsigned long> stripes;

    //cout << "refBitLen = " << refBitLen << endl;
    //cout << "refBinSeq = " << (int)*refBinSeq << endl;
    for (int shift=0; shift<hasher[0].charLen*8; shift+=2) {
        //cout << (shift&0x0e) << "= shift";
        int n = (refBitLen-shift)/hasher[0].charLen/8-hasher[0].charsInWord;
        //cout << " n = " << n << endl;
        if (n>=0) {
            unsigned int hash = hasher[0].singleHash(((basechartype*)refBinSeq)+(shift>>4),shift&0x0e);
            insert(stripes,hash,(shift>>4)*2,shift&0x0e);
            //cout << "hash = " << hash << "   byte = " << 0 << "    shift = " << shift << endl;
        
            for (int i = 0; i<n*hasher[0].charLen; i+=hasher[0].charLen) {
                hasher[0].moveRight(&hash,(((basechartype*) (refBinSeq+i))+(shift>>4)), shift&0x0e);
                insert(stripes,hash,i+hasher[0].charLen+(shift>>4)*2,shift&0x0e);
                cout << "hash = " << hash << "   byte = " << i+hasher[0].charLen+(shift>>4)*2 << "    shift = " << (shift&0x0e) << endl;
            }
        }
    }
    
    //for (map<unsigned int,unsigned long>::iterator it = stripes.begin(); it!=stripes.end(); it++) {
        //cout << "stripe:" << (*it).first << " " << (*it).second << endl;   
    //}

    // following command line arguments : other files containing sequences
    // result is stored in an associative array ordered by title
    //map<string,string> otherSequences;
    long threadLen = 0;
    for (int i=4; i<argc; i++) {
        threadLen += Reader::lookLen(argv[i]);
    }
    //cout << "threadLen = " << threadLen << endl;
    threadLen /= numberOfWorkingThreads;
    threadLen++;
    
    //cout << "threadLen = " << threadLen << endl;
    
    vector<Reader::mapped> *files = new vector<Reader::mapped>[numberOfWorkingThreads];
    vector<Reader::mapped> *bins = new vector<Reader::mapped>[numberOfWorkingThreads];
    vector<string> *names = new vector<string>[numberOfWorkingThreads];
    vector<foundation> *found = new vector<foundation>[numberOfWorkingThreads];
    vector<int> *intnames = new vector<int>[numberOfWorkingThreads];
    string *firstnames = new string[numberOfWorkingThreads+1];   
    long *prefixes = new long[numberOfWorkingThreads+1];
    bool *dosumprefix = new bool[numberOfWorkingThreads+1];
    long *whereToStart = new long[numberOfWorkingThreads+1];
    long *whereToStop = new long[numberOfWorkingThreads+1];
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
                
                //cout << "j = " << j << endl;
                //cout << "i = " << i << endl;
                //cout << "thisThreadLen = " << thisThreadLen << endl;
                //cout << " whereToStop[j] = " <<  whereToStop[j] << endl;
                //cout << " whereToStart[j+1] = " <<  whereToStart[j+1] << endl;
                j++;
            } else {
                thisThreadLen += refFile.len;
                if (i+1<argc) {
                    refFile = Reader::mapfile(argv[i+1]);
                }
                //cout << "i = " << i << endl;
            }
        }
        whereToStop[j] = refFile.len;
        //cout << " whereToStop[j] = " <<  whereToStop[j] << endl;
    }
    
    #pragma omp parallel for
    for (int thread = 0; thread<numberOfWorkingThreads; thread++) {
        names[thread].push_back("firstname");
        for (int file = 0; file<files[thread].size(); file++) {
            Reader::mapped refFile = files[thread][file];
            long start = file==0 ? whereToStart[thread] : 0;
            long stop = file==files[thread].size()-1 ? whereToStop[thread] : refFile.len;
            long newstop = min(stop+4*hashLen*4,(unsigned long)refFile.len);
            //cout << "thread = " << thread << endl;
            //cout << "start = " << start << endl;
            //cout << "stop = " << stop << endl;            
            size_t refBinLen = (stop-start)/32*32+1024;
            char *refBinSeq = hasher[thread].charLen + (char *) ippsMalloc_16u(refBinLen/2+hasher[thread].baseCharsInChar);
            pair<Reader::mapped,long> bin1;
            Reader::mapped bin;
            bin.pnt = refBinSeq;
            dosumprefix[thread+1] = true;
            while (stop>start) {
                //cout << "start = " << start << endl;
                bin1 = reader.readData((unsigned char *) refFile.pnt+start,(unsigned char *)refBinSeq,newstop-start,stop-start);                
                bin.len = bin1.first.len;
                //cout << " bin.len = " << bin.len << endl;
                //cout << "bin.pnt = " << (int*)refBinSeq << endl;
                bin.pnt = refBinSeq;
                bins[thread].push_back(bin);
                refBinSeq += bin1.first.len/4/32*32+32;
                start = (bin1.first.pnt-refFile.pnt);
                //cout << "start1 = " << start << endl;
                if (stop>start) {
                    pair<string,char *> namepair = reader.readName(bin1.first.pnt);
                    names[thread].push_back(namepair.first);
                    firstnames[thread+1] = namepair.first;
                    start = namepair.second-refFile.pnt;
                    dosumprefix[thread+1] = false;                    
                } else {
                    prefixes[thread+1] = bin1.second;
                }
                //cout << "start2 = " << start << endl;
            }         
                      
            
            
            //long refCharLen = reader.readData((unsigned char *) namepair.second,(unsigned char *)refBinSeq,((int)(refFile.pnt-namepair.second+refFile.len)));
            //assert(refBinLen >= (refCharLen)/128*32+32);
            //for (long len = start; len < stop; len++){
            //    putchar (files[thread][file].pnt[len]);   
            //}
        }
    }

    //do not parallel this cycle
    for (int thread = 1; thread<numberOfWorkingThreads; thread++) {
        //cout << "thread = " << thread << endl;
        //cout << "firstname = " << firstnames[thread] << endl;
        //cout << "prefix = " << prefixes[thread] << endl;
        //cout << "dosumprefix = " << dosumprefix[thread] << endl;
        if (dosumprefix[thread]) {
            prefixes[thread]+=prefixes[thread-1];
            firstnames[thread] = firstnames[thread-1];
            //cout << "new prefix = " << prefixes[thread] << endl;
        }
    }
    
    
    map<string,int> seqNames;
       
    //do not parallel this cycle
    for (int thread = 0; thread<numberOfWorkingThreads; thread++) {
        //cout << firstnames[thread] << endl;
        names[thread][0] = firstnames[thread];
        for (int name=0; name<names[thread].size(); name++) {
            seqNames.insert(pair<string,int>(names[thread][name],0));
        }
    }
    
    //do not parallel this cycle
    {
        int i = 0;
        for (map<string,int>::iterator it = seqNames.begin(); it!=seqNames.end(); it++,i++) {
            seqNames[(*it).first] = i;
        }
    }
    
    #pragma omp parallel for
    for (int thread=0; thread<numberOfWorkingThreads; thread++) {
        for (int i=0; i<names[thread].size(); i++) {
            intnames[thread].push_back(seqNames[names[thread][i]]);
        }
    }
    
    vector<pair<string,int> > seqVector(seqNames.begin(), seqNames.end());
        
    #pragma omp parallel for
    for (int thread = 0; thread<numberOfWorkingThreads; thread++) {
        //cout << "thread = " << thread << endl;
        for (int ibin = 0; ibin<bins[thread].size(); ibin++) {
            Reader::mapped bin = bins[thread][ibin];
            unsigned int hash;
            for (long i=0; i<bin.len/4; i+=hashLen) {
                //cout << "ibin = " << ibin << endl;
                //cout << "pos = " <<  i*4+(ibin==0?prefixes[thread]:0) << endl;
                //cout << "bin.len = " << bin.len << endl;
                //cout << (int)(unsigned char)bin.pnt[i] << endl;
                hash = hasher[thread].singleHash((basechartype *) (bin.pnt+i));
                cout << "hash = " << hash << endl;
                //multimap<unsigned int, unsigned long>::iterator it = stripes.find(hash);
                pair<multimap<unsigned int, unsigned long>::iterator,multimap<unsigned int, unsigned long>::iterator> range = stripes.equal_range(hash);
                for (multimap<unsigned int, unsigned long>::iterator it=range.first; it!=range.second; it++) {
                    long thispos = i*4+(ibin==0?prefixes[thread]:0);
                    char* thispnt = bin.pnt+i;
                    long thatpos = (*it).second;
                    long realthatpos = (thatpos>>6)+((thatpos&0xff)>>1);
                    cout << "found at  " << thispos << " :" << realthatpos << " | " << seqVector[intnames[thread][ibin]].first << endl;
                    cout << "that = " << (short)(unsigned char)(*(refBinSeq+(thatpos>>8))) << endl;
                    cout << "this = " << (short)(unsigned char)(*(thispnt)) << endl;
                    long forwlen = min((bin.len-i*4)+0ul,(refBitLen)/2-realthatpos);
                    long backwlen = min(i*4,realthatpos);
                    long endpos = hasher[thread].compareForward((basechartype *) thispnt, (basechartype *)(refBinSeq+(thatpos>>8)), thatpos&0xff, forwlen);
                    long startpos = hasher[thread].compareBackward((basechartype *) thispnt, (basechartype *)(refBinSeq+(thatpos>>8)), thatpos&0xff, backwlen);
                    cout << "thispos = " << thispos << endl;
                    cout << "startpos = " << startpos << endl;                    
                    cout << "endpos = " << endpos << endl;
                    endpos = min(endpos,forwlen);
                    startpos = max(startpos,-backwlen);
                    cout << "startpos = " << startpos << endl;                    
                    cout << "endpos = " << endpos << endl;
                    endpos += thispos;
                    startpos += thispos;
                    cout << "startpos = " << startpos << endl;                    
                    cout << "endpos = " << endpos << endl;
                    if ((endpos-startpos)/4>=hashLen) {
                        foundation f;
                        f.seqNum = intnames[thread][ibin];
                        f.seqStartPos = startpos;
                        f.seqEndPos = endpos;
                        f.refStartPos = (thatpos>>6)+((thatpos&0xff)>>1)-thispos+startpos;
                        found[thread].push_back(f);
                        //i+=((endpos-thispos)/4/hashLen-1)*hashLen;
                    }
                }
            }
        }
    }
    
    //old dummy code
    /*{
    int thr_1 = 0;        
    for (int thread = 1; thread<numberOfWorkingThreads; thread++) {
        if (found[thr_1].size()>0) {
            if (found[thread].size()>0) {                
                foundation &fp = found[thr_1][found[thr_1].size()-1];
                foundation &ft = found[thread][0];
                
                if (fp.seqNum == ft.seqNum) {
                    cout << "thread = " << thread << endl;
                    if (fp.seqEndPos>=ft.seqStartPos) {
                        if (fp.seqStartPos-fp.refStartPos==ft.seqStartPos-ft.refStartPos) {
                            
                            cout << "seqStartPos = " << fp.seqStartPos << endl;
                            cout << "found[thread].size() = " << found[thread].size() << endl;
                            
                            ft.seqStartPos = fp.seqStartPos;
                            ft.refStartPos = fp.refStartPos;
                            found[thr_1].pop_back();                            
                        }
                    }
                }
                thr_1 = thread;
            }            
        } else if (found[thread].size()>0) {
            thr_1 = thread;
        }
    }
    }*/
    
    
    vector<foundation> answer;
    
    for (int thread = 0; thread<numberOfWorkingThreads; thread++) {
        answer.insert(answer.end(),found[thread].begin(),found[thread].end());       
    }
    
    sort(answer.begin(),answer.end(),compareprev);
    
    for (int ia = 1; ia<answer.size(); ia++) {
        if (answer[ia-1].seqNum==answer[ia].seqNum) {
            if (answer[ia-1].seqStartPos-answer[ia-1].refStartPos==answer[ia].seqStartPos-answer[ia].refStartPos) {
                if (answer[ia-1].seqEndPos>=answer[ia].seqStartPos) {
                    answer[ia].seqStartPos = answer[ia-1].seqStartPos;
                    answer[ia].refStartPos = answer[ia-1].refStartPos;
                    answer[ia-1].seqNum = -1;
                }
            }
        }
    }
    
    sort(answer.begin(),answer.end(),compareans);
    
    long name1 = 0, name2;
    for (vector<foundation>::iterator it = answer.begin(); it!=answer.end(); it++) {
        foundation itt = (*it);
        name2 = itt.seqNum;
        if (name2!=-1 && itt.seqEndPos-itt.seqStartPos>=minMatchCharLen) {
            if (name2!=name1) {
                for (int i=name1+1; i<=name2; i++) {
                    cout << seqVector[i].first << endl;
                }
                name1 = name2;                            
            }
            //cout << ((*(++it)).first < (*(it)).first) << endl;
            //it--;
            cout << itt.refStartPos+1 << " " << itt.refStartPos+itt.seqEndPos-itt.seqStartPos << " " << itt.seqStartPos+1 << " " << itt.seqEndPos << endl;
        }
    }
    for (int i=name1+1; i<seqVector.size(); i++) {
        cout << seqVector[i].first << endl;
    }
    
    
    
    //cout << endl << endl;
    return 0;
}
