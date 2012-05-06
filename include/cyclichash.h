#ifndef CYCLICHASH
#define CYCLICHASH

#include <cmath>
#include <assert.h>
#include <stdlib.h>
#include "/home/mburyakov/bin/intel/ipp/include/ipp.h"
#include "/home/mburyakov/bin/intel/ipp/include/ipps.h"

using namespace std;


typedef unsigned int hashanstype;
typedef Ipp16u basechartype;

class CyclicHash {

    //all Lens are in bytes

    size_t wordLen;
    const static size_t baseCharLen = 2;
    size_t charLen;
    size_t baseCharsInWord;
    size_t baseCharsInChar;
    size_t charsInWord;
    size_t charSpace; // pow(2,(charLen*charLen))



    //unsigned char charPermute[];

    public:
        CyclicHash(size_t wordLen, size_t charLen) {
            this->wordLen = wordLen;
            this->charLen = charLen;
            assert(wordLen%charLen==0);
            this->charsInWord = wordLen/charLen;
            this->baseCharsInWord = wordLen/baseCharLen;
            this->baseCharsInWord = charLen/baseCharLen;
            charSpace == 1<<(charLen*8);
            //charPermute = new unsigned char[charSpace];
            //srand(time(NULL));
            //for (int i=0; i<charSpace; i++) {
            //    charPermute[i] = (unsigned char) rand();
            //}
        }

        hashanstype singleHash(basechartype src[]) {
            Ipp16u *pnt1 = ippsMalloc_16u(wordLen);
            Ipp16u *pnt2 = ippsMalloc_16u(wordLen);
            Ipp16u *pnt3 = ippsMalloc_16u(wordLen);

            //cyclic shifting
            ippsCopy_16s((Ipp16s *)src, (Ipp16s *)pnt1, wordLen);
            ippsCopy_16s((Ipp16s *)src, (Ipp16s *)pnt2, wordLen);
            for (size_t i=0; i<charsInWord; i++) {
                ippsLShiftC_16s_I(i,(Ipp16s *)pnt1+i*baseCharsInChar,baseCharsInChar);
                ippsRShiftC_16u_I(16-i,pnt2+i*baseCharsInChar,baseCharsInChar);

            }

            ippsOr_8u((Ipp8u *)pnt1,(Ipp8u *)pnt2,(Ipp8u *)pnt3,wordLen);

            //TODO:test and calc hash
        }

    /*CyclicHash(int myn, int mywordsize=19) : hashvalue(0),
    n(myn), wordsize(mywordsize),
      hasher( ( 1<<wordsize ) - 1),
      mask1((static_cast<hashvaluetype>(1)<<(wordsize-1)) -1),
      myr(n%wordsize),
      maskn((static_cast<uint32>(1)<<(wordsize-myr)) -1 )
      {
       if(static_cast<uint>(wordsize) > 8*sizeof(hashvaluetype)) {
       cerr<<"Can't create "<<wordsize<<"-bit hash values"<<endl;
       throw "abord";
      }
    }

    inline void fastleftshiftn(hashvaluetype & x) const {
        x = ((x & maskn) << myr ) | (x >> (wordsize-myr)) ; }

    inline void fastleftshift(hashvaluetype & x, int r) const {
        r = r % wordsize;
        const uint32 mask = (static_cast<uint32>(1)<<(wordsize-r)) -1 ;
        x = ((x & mask) << r ) | (x >> (wordsize-r)) ;
    }


    inline void fastleftshift1(hashvaluetype & x) const {
        x = ((x & mask1) << 1 ) | (x >> (wordsize-1)) ;
    }


    template<class container>
    hashvaluetype hash(container & c) {
     assert(c.size()==static_cast<uint>(n));
     hashvaluetype answer(0);
     for(uint k = 0; k<c.size();++k) {
     fastleftshift(answer, 1) ;
     answer ^= hasher.hashvalues[c[k]];
     }
     return answer;
    }


    inline void update(chartype outchar, chartype inchar) {
      hashvaluetype z (hasher.hashvalues[outchar]);
      fastleftshiftn(z);
      hashvalue = ( ((hashvalue & mask1) << 1 ) | (hashvalue >> (wordsize-1)) )
      ^ z
      ^ hasher.hashvalues[inchar];
    }



    void eat(chartype inchar) {
      fastleftshift1(hashvalue);
      hashvalue ^= hasher.hashvalues[inchar];
    }

    uint32 hashvalue;
    const int n, wordsize;
    CharacterHash hasher;
    const hashvaluetype mask1;
    const int myr;
    const hashvaluetype maskn;    */

};



#endif