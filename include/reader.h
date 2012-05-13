#ifndef CYCLICHASH
#define CYCLICHASH

#include <cmath>
#include <assert.h>
#include <stdlib.h>
#include <assert.h>
#include <ipp.h>

class Reader {

        const A2 = "AA";
        const A16 = "AAAAAAAAAAAAAAAA";
        const A32 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        const A96 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        const C16 = "CCCCCCCCCCCCCCCC";
        const C32 = "CCCCCCCCCCCCCCCC";

        Ipp16u *tmp1;
        Ipp16u *tmp2;
        const size_t baseCharsInBlock = 16;

        Reader() {
            tmp1 = new ippsMalloc_16u(baseCharsInBlock*2);
            tmp2 = tmp1 + baseCharsInBlock;
        }

        //returns Charlength of letters
        int readData(char letters[], char dest[]) {
            int j = 0;
            for (int i=0; letters[i]!='>'; i++) {
                if (letters[i]=='A'||letters[i]=='C'||letters[i]=='T'||letters[i]=='G') {
                    dest[j>>3] |= (letters[i]&0x06)<<(5-(shift&0x06));
                    j+=2;
                }
            }
            return j/2;
        }



}

#endif