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

    const basechartype tmp[];
    const basechartype tmp1[];
    const basechartype tmp2[];
    const basechartype tmp3[];
    const basechartype tmp4[];
    const basechartype accum[];

    //unsigned char charPermute[];

    public:
        CyclicHash(size_t wordLen, size_t charLen) {
            this->wordLen = wordLen;
            this->charLen = charLen;
            assert(wordLen%charLen==0);
            assert(charLen%baseCharLen==0);
            assert(charLen>=sizeof(hasnanstype));
            this->charsInWord = wordLen/charLen;
            this->baseCharsInWord = wordLen/baseCharLen;
            this->baseCharsInChar = charLen/baseCharLen;
            this->charsInWord = wordLen/charLen;
            charSpace == 1<<(charLen*8);
            tmp = ippsMalloc_16u(baseCharsInChar+1);
            tmp1 = ippsMalloc_16u(baseCharsInChar+1);
            tmp2 = ippsMalloc_16u(baseCharsInChar+1);
            tmp3 = ippsMalloc_16u(baseCharsInChar+1);
            tmp4 = ippsMalloc_16u(baseCharsInChar+1);
            accum = ippsMalloc_16u(baseCharsInChar+1);
            //charPermute = new unsigned char[charSpace];
            //srand(time(NULL));
            //for (int i=0; i<charSpace; i++) {
            //    charPermute[i] = (unsigned char) rand();
            //}
        }

        inline static void cyclicShiftChar(basechartype const src[],basechartype dest[],size_t i) {
            ippsCopy_16s((Ipp16s *)src, (Ipp16s *)tmp1, baseCharsInChar);
            ippsCopy_16s((Ipp16s *)src, (Ipp16s *)tmp2, baseCharsInChar);
            ippsLShiftC_16s_I(i,(Ipp16s *)tmp1,baseCharsInChar);
            ippsRShiftC_16u_I(16-i,tmp2,baseCharsInChar);
            ippsOr_8u((Ipp8u *)tmp1,(Ipp8u *)tmp2,(Ipp8u *)dest,CharLen);
        }

        inline static void cyclicShiftChar(basechartype const src[],basechartype dest[],size_t i,size_t shift) {
            ippsCopy_16s((Ipp16s *)src, (Ipp16s *)tmp1, baseCharsInChar+1);
            ippsCopy_16s((Ipp16s *)src, (Ipp16s *)tmp2, baseCharsInChar+1);
            if (shift+i < 16) {
                ippsLShiftC_16s_I(i+shift,(Ipp16s *)tmp1,baseCharsInChar);
                ippsRShiftC_16u_I(16-i-shift,tmp2,baseCharsInChar);
                ippsCopy_16s((Ipp16s *)tmp2+1, (Ipp16s *)tmp3, baseCharsInChar);
                Ipp16u mask = ((1<<shift)-1)<<i;
                ippsAndC_16u_I(mask,tmp3,baseCharsInChar)  //TODO: use avx function
                mask = ~mask;
                ippsAndC_16u_I(mask,tmp2,baseCharsInChar)  //TODO: use avx function
                ippsOr_8u((Ipp8u *)tmp1,(Ipp8u *)tmp2,(Ipp8u *)dest,CharLen);
                ippsOr_8u_I((Ipp8u *)tmp3,(Ipp8u *)dest,CharLen);  //TODO: use avx function
            } else {
                mask = (1<<i)-1;
                ippsLShiftC_16s_I(i+shift-16,(Ipp16s *)tmp1,baseCharsInChar);
                ippsAndC_16u_I(mask,tmp1,baseCharsInChar)  //TODO: use avx function
                ippsRShiftC_16u_I(32-i-shift,tmp2,baseCharsInChar);
                ippsCopy_16s((Ipp16s *)tmp2+1, (Ipp16s *)tmp3, baseCharsInChar);
                ippsOr_8u((Ipp8u *)tmp3,(Ipp8u *)tmp2,(Ipp8u *)dest,CharLen);
                ippsCopy_16s((Ipp16s *)tmp1+1, (Ipp16s *)tmp3, baseCharsInChar);
                mask = ~mask;
                ippsAndC_16u_I(mask,tmp3,baseCharsInChar)  //TODO: use avx function
                ippsOr_8u_I((Ipp8u *)tmp3,(Ipp8u *)dest,CharLen); //TODO: use avx function
            }
        }

        inline static void deShiftChar(basechartype const src[],basechartype dest[],size_t shift) {
            ippsCopy_16s((Ipp16s *)src, (Ipp16s *)tmp1, baseCharsInChar+1);
            ippsCopy_16s((Ipp16s *)src, (Ipp16s *)tmp2, baseCharsInChar+1);
            ippsLShiftC_16s_I(shift,(Ipp16s *)tmp1,baseCharsInChar);
            ippsRShiftC_16u_I(16-shift,tmp2,baseCharsInChar);
            ippsOr_8u((Ipp8u *)tmp1,(Ipp8u *)tmp2,(Ipp8u *)dest,CharLen);
        }


        hashanstype singleHash(basechartype const src[]) {
            //cyclic shifting
            ippsZero_16s((Ipp16s *)accum,baseCharsInChar);
            for (size_t i=0, size_t ibcc=0; i<charsInWord; i++, ibcc+=baseCharsInChar) {
                // TODO: some random permutation
                cyclicShiftChar(src+ibcc,tmp,i);
                ippsXor_16u_I(tmp,accum,baseCharsInChar); //TODO: use avx function
            }
            //TODO:test and calc hash
        }

        hashanstype collapseChar(basechartype const src[]) {
            return *((hashanstype *) src); //TODO: use some linear collapse
        }

        hashanstype singleHash(basechartype const src[], size_t shift) {
            //cyclic shifting
            ippsZero_16s((Ipp16s *)accum,baseCharsInChar);
            for (size_t i=0, size_t ibcc=0; i<charsInWord; i++, ibcc+=baseCharsInChar) {
                // TODO: some random permutation
                cyclicShiftChar(src+ibcc,tmp,i,shift);
                ippsXor_16u_I(tmp,accum,baseCharsInChar); //TODO: use avx function
            }
            return collapseChar(accum);
            //TODO:test and calc hash
        }

        void moveRight(hashanstype &prev, basechartype const src[]) {
            // TODO: some random permutation
            cyclicShiftChar(src,tmp,charsInWord);
            ippXor_16u_I(src+baseCharsInWord, tmp, baseCharsInChar); //TODO: use avx function
            prev ^= collapse(tmp);
        }

        void moveRight(hashanstype &prev, basechartype const src[], size_t shift) {
            // TODO: some random permutation
            cyclicShiftChar(src,tmp,charsInWord,shift);
            deShiftChar(src+baseCharsInWord,accum,shift);
            ippXor_16u_I(accum, tmp, baseCharsInChar); //TODO: use avx function
            prev ^= collapse(tmp);
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