#ifndef CYCLICHASH
#define CYCLICHASH

#include <cmath>
#include <assert.h>
#include <stdlib.h>
#include <assert.h>
#include <ipp.h>
//#include "/home/mburyakov/bin/intel/ipp/include/ipps.h"

using namespace std;


typedef unsigned int hashanstype;
typedef Ipp16u basechartype;

class CyclicHash {
    public:
    //all Lens are in bytes

    size_t wordLen;
    const static size_t baseCharLen = 2;
    size_t charLen;
    size_t baseCharsInWord;
    size_t baseCharsInChar;
    size_t baseCharsInAns;
    size_t ansLen;
    size_t charsInWord;
    size_t charSpace; // pow(2,(charLen*charLen))

    basechartype *tmp;
    basechartype *tmp1;
    basechartype *tmp2;
    basechartype *tmp3;
    basechartype *tmp4;
    basechartype *accum;

    //unsigned char charPermute[];


        CyclicHash(size_t wordLen, size_t charLen) {
            this->wordLen = wordLen;
            this->charLen = charLen;
            assert(wordLen%charLen==0);
            assert(charLen%baseCharLen==0);
            //assert(charLen>=sizeof(hashanstype));
            this->charsInWord = wordLen/charLen;
            this->baseCharsInWord = wordLen/baseCharLen;
            this->baseCharsInChar = charLen/baseCharLen;
            this->ansLen = sizeof(hashanstype);
            assert(ansLen%baseCharLen==0);
            this->baseCharsInAns = ansLen/baseCharLen;
            //assert(baseCharsInChar%2==0);
            this->charsInWord = wordLen/charLen;
            charSpace = 1<<(charLen*8);
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

        inline void cyclicShiftChar(basechartype const src[],basechartype dest[],size_t i) {
            ippsSwapBytes_16u((Ipp16u *)src, (Ipp16u *)tmp1, baseCharsInChar);
            ippsCopy_16s((Ipp16s *)tmp1, (Ipp16s *)tmp2, baseCharsInChar);
            ippsLShiftC_16s_I(i,(Ipp16s *)tmp1,baseCharsInChar);
            ippsRShiftC_16u_I(16-i,tmp2,baseCharsInChar);
            ippsOr_8u((Ipp8u *)tmp1,(Ipp8u *)tmp2,(Ipp8u *)dest,charLen);
        }

        inline void cyclicShiftAns_I(basechartype dest[],size_t i) {
            ippsCopy_16s((Ipp16s *)dest, (Ipp16s *)tmp1, baseCharsInAns);
            ippsCopy_16s((Ipp16s *)dest, (Ipp16s *)tmp2, baseCharsInAns);
            ippsLShiftC_16s_I(i,(Ipp16s *)tmp1,baseCharsInAns);
            ippsRShiftC_16u_I(16-i,tmp2,baseCharsInAns);
            ippsOr_8u((Ipp8u *)tmp1,(Ipp8u *)tmp2,(Ipp8u *)dest,ansLen);
        }

        inline void cyclicShiftChar(basechartype const src[],basechartype dest[],size_t i,size_t shift) {
            ippsSwapBytes_16u((Ipp16u *)src, (Ipp16u *)tmp1, baseCharsInChar+1);
            ippsCopy_16s((Ipp16s *)tmp1, (Ipp16s *)tmp2, baseCharsInChar+1);
            if (shift+i < 16) {
                ippsLShiftC_16s_I(i+shift,(Ipp16s *)tmp1,baseCharsInChar);
                //cout << "src + 1 = "<< *((unsigned int*) (src+1)) << endl;
                //cout << "src = "<< *((unsigned int*) (src)) << endl;
                //cout << "tmp2 = "<< *((unsigned int*) (tmp2)) << endl;
                //cout << "16-i-shift" << 16-i-shift;
                ippsRShiftC_16u_I(16-i-shift,tmp2,baseCharsInChar+1);
                //cout << "tmp2 = "<< *((unsigned int*) (tmp2)) << endl;
                //cout << "tmp2 + 1 = "<< *((unsigned int*) (tmp2 + 1)) << endl;
                ippsCopy_16s((Ipp16s *)tmp2+1, (Ipp16s *)tmp3, baseCharsInChar);
                Ipp16u mask = ((1<<shift)-1)<<i;
                //cout << "tmp2 = "<< *((unsigned int*) (tmp2)) << endl;
                //cout << "tmp3 = "<< *((unsigned *) tmp3) << endl;
                ippsAndC_16u_I(mask,tmp3,baseCharsInChar);  //TODO: use avx function
                mask = ~mask;
                ippsAndC_16u_I(mask,tmp2,baseCharsInChar);  //TODO: use avx function
                //cout << "tmp1 = "<< *((unsigned*) tmp1) << endl;
                //cout << "tmp2 = "<< *((unsigned *) tmp2) << endl;
                //cout << "tmp3 = "<< *((unsigned *) tmp3) << endl;
                ippsOr_8u((Ipp8u *)tmp1,(Ipp8u *)tmp2,(Ipp8u *)dest,charLen);
                ippsOr_8u_I((Ipp8u *)tmp3,(Ipp8u *)dest,charLen);  //TODO: use avx function
            } else {
                Ipp16u mask = (1<<i)-1;
                ippsLShiftC_16s_I(i+shift-16,(Ipp16s *)tmp1,baseCharsInChar+1);
                ippsAndC_16u_I(mask,tmp1,baseCharsInChar+1);  //TODO: use avx function
                ippsRShiftC_16u_I(32-i-shift,tmp2,baseCharsInChar+1);
                ippsCopy_16s((Ipp16s *)tmp2+1, (Ipp16s *)tmp3, baseCharsInChar);
                ippsOr_8u((Ipp8u *)tmp3,(Ipp8u *)tmp2,(Ipp8u *)dest,charLen);
                ippsCopy_16s((Ipp16s *)tmp1+1, (Ipp16s *)tmp3, baseCharsInChar);
                mask = ~mask;
                ippsAndC_16u_I(mask,tmp3,baseCharsInChar);  //TODO: use avx function
                ippsOr_8u_I((Ipp8u *)tmp3,(Ipp8u *)dest,charLen); //TODO: use avx function
            }
        }
        
        inline long compareForward(basechartype const *src1, basechartype const *src2, size_t shift, long len) {
            //cout << "src2+0 = " << src2[0] << endl;
            //cout << "src1+0 = " << src1[0] << endl;
            for (int i=0; i<len/4+1; i+=baseCharsInChar) {
                //cout << "len = " << len << endl;
                //cout << "kyky0" << endl;
                //cout << "tmp+0 = " << tmp[i] << endl;
                deShiftChar(src2+i,tmp,shift);
                //cout << "tmp+0 = " << tmp[i] << endl;
                //cout << "src2+i = " << src2[i] << endl;
                ippsSwapBytes_16u(src1+i,accum,baseCharsInChar);
                ippsXor_16u_I(accum,tmp,baseCharsInChar);
                //cout << "shift = " << shift << endl;
                //cout << "tmp+0 = " << tmp[0] << endl;
                //cout << "src2+0 = " << src2[0] << endl;
                //cout << "src1+0 = " << src1[0] << endl;
                for (int j=0; j<baseCharsInChar; j++) {
                    if (tmp[j]!=0) {
                        
                        int stop;
                        //cout << "shift = " << shift << endl;
                        //cout << "tmp+0 = " << tmp[0] << endl;
                        //cout << "src1+0 = " << src1[0] << endl;
                        for (stop=0; tmp[j]!=0; stop+=1,tmp[j]>>=2);
                        //cout << "kuku" << endl;
                        stop = 8-stop+4*baseCharLen*j+4*baseCharLen*i;
                        return stop;
                    }
                }
            }
            return len;
        }
        
        inline long compareBackward(basechartype const *src1, basechartype const *src2, size_t shift, long len) {
            cout << "len = " << len << endl;
            for (int i=-baseCharsInChar; -i<len/4+1; i-=baseCharsInChar) {
                
                cout << "kyky0" << endl;
                cout << "tmp+0 = " << tmp[0] << endl;
                deShiftChar(src2+i,tmp,shift);
                cout << "tmp+0 = " << tmp[0] << endl;
                cout << "src2+i = " << src2[i] << endl;
                ippsSwapBytes_16u(src1+i,accum,baseCharsInChar);
                ippsXor_16u_I(accum,tmp,baseCharsInChar);
                cout << "shift = " << shift << endl;
                cout << "tmp+0 = " << tmp[0] << endl;
                cout << "src2+i = " << src2[i] << endl;
                cout << "src1+i = " << src1[i] << endl;
                for (int j=baseCharsInChar-1; j>=0; j--) {
                    if (tmp[j]!=0) {
                        
                        int stop;
                        cout << "shift = " << shift << endl;
                        cout << "tmp+i = " << tmp[i] << endl;
                        cout << "src1+i = " << src1[i] << endl;
                        for (stop=0; tmp[j]!=0; stop+=1,tmp[j]<<=2);
                        cout << "kuku" << "  j = " << j << "   i=" << i << endl;
                        stop = stop+4*baseCharLen*j+4*baseCharLen*i;
                        return stop;
                    }
                }
            }
            return -len;
        }        

        inline void deShiftChar(basechartype const src[],basechartype dest[],size_t shift) {
            if (shift == 0) {
                ippsSwapBytes_16u((Ipp16u *)src, (Ipp16u *)dest, baseCharsInChar);
                //cout << "src[0] = " << src[0];
                //cout << "dest[0] = " << dest[0];
            } else {
                ippsSwapBytes_16u((Ipp16u *)src, (Ipp16u *)tmp1, baseCharsInChar+1);
                ippsCopy_16s((Ipp16s *)tmp1, (Ipp16s *)tmp2, baseCharsInChar+1);
                ippsLShiftC_16s_I(shift,(Ipp16s *)tmp1,baseCharsInChar);
                ippsRShiftC_16u_I(16-shift,tmp2,baseCharsInChar+1);
                //cout << "deshiftR = " << tmp2[1] << endl;
                ippsOr_8u((Ipp8u *)tmp1,(Ipp8u *)(tmp2+1),(Ipp8u *)dest,charLen);
            }
        }

        hashanstype collapseChar(basechartype const src[]) {
            hashanstype ans = 0;  // TODO: use vectorization
            for (int i=0; i<baseCharsInChar/2; i+=2) {
                ippsXor_32u_I((Ipp32u *)(src+i),&ans,1);
            }
            return ans;
        }


        virtual hashanstype singleHash(basechartype const src[]) {
            //cyclic shifting
            ippsZero_16s((Ipp16s *)accum,baseCharsInChar);
            {
                size_t ibcc;
                size_t i;
                for ((i=0), ibcc=0; i<charsInWord; i++, ibcc+=baseCharsInChar) {
                    // TODO: some random permutation
                    cyclicShiftChar(src+ibcc,tmp,i);
                    ippsXor_16u_I(tmp,accum,baseCharsInChar); //TODO: use avx function
                }
            }
            return collapseChar(accum);
        }



        virtual hashanstype singleHash(basechartype const src[], size_t shift) {
            //cyclic shifting
            ippsZero_16s((Ipp16s *)accum,baseCharsInChar);
            {
                size_t ibcc;
                size_t i;
                for ((i=0), ibcc=0; i<charsInWord; i++, ibcc+=baseCharsInChar) {
                    // TODO: some random permutation
                    cyclicShiftChar(src+ibcc,tmp,i,shift);
                    ippsXor_16u_I(tmp,accum,baseCharsInChar); //TODO: use avx function
                }
            }
            return collapseChar(accum);
        }

        virtual void moveRight(hashanstype *prev, basechartype const src[]) {
            // TODO: some random permutation
            cyclicShiftChar(src+baseCharsInWord,tmp,charsInWord);
            ippsXor_16u_I(src, tmp, baseCharsInChar); //TODO: use avx function
            *prev ^= collapseChar(tmp);
            cyclicShiftAns_I((basechartype *)prev,baseCharLen*8-1);
            //cout << "charsInWord = " << charsInWord << endl;

        }

        virtual void moveRight(hashanstype *prev, basechartype const src[], size_t shift) {
            // TODO: some random permutation
            cyclicShiftChar(src+baseCharsInWord,tmp,charsInWord,shift);
            deShiftChar(src,accum,shift);
            ippsXor_16u_I(accum, tmp, baseCharsInChar); //TODO: use avx function
            *prev ^= collapseChar(tmp);
            cyclicShiftAns_I((basechartype *)prev,baseCharLen*8-1);
        }


};



class TinyHash:public CyclicHash {
    public:

        
        TinyHash():CyclicHash(baseCharLen,baseCharLen) {
            
        }
    
        virtual hashanstype singleHash(basechartype const src[], size_t shift) {
            ippsSwapBytes_16u((Ipp16u *)src, (Ipp16u *)tmp1, baseCharsInChar+1);
            return (*((unsigned int *)tmp1)<<shift)>>16;
        }
        
        virtual hashanstype singleHash(basechartype const src[]) {
            cout << "kuku" << endl;
            ippsSwapBytes_16u((Ipp16u *)src, (Ipp16u *)tmp1, baseCharsInChar+1);
            return *tmp1;
        }
        virtual void moveRight(hashanstype *prev, basechartype const src[]) {            
            *prev=singleHash(src+1);
        }
        virtual void moveRight(hashanstype *prev, basechartype const src[], size_t shift) {            
            *prev=singleHash(src+1, shift);
        }
};



#endif