#ifndef READER
#define READER

#include <cmath>
#include <assert.h>
#include <stdlib.h>
#include <assert.h>
#include <ipp.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

        const char* A2 = "AA";
        const char* A16 = "AAAAAAAAAAAAAAAA";
        const char* A32 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        const char* A96 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        const char* C16 = "CCCCCCCCCCCCCCCC";
        const char* C32 = "CCCCCCCCCCCCCCCC";

class Reader {
    
public:



        Ipp16u *tmp1;
        Ipp16u *tmp2;
        static const size_t baseCharsInBlock = 16;

        Reader() {
            tmp1 = ippsMalloc_16u(baseCharsInBlock*2);
            tmp2 = tmp1 + baseCharsInBlock;
        }

struct mapped {
    char *pnt;
    long len;
};        
        
  static mapped mapfile(char *filename) {
          struct stat sb;
          char *p;
          int fd;

          fd = open (filename, O_RDONLY);
          if (fd == -1) {
                  perror ("open");
                  throw "error mapping file";
          }

          if (fstat (fd, &sb) == -1) {
                  perror ("fstat");
                  throw "error mapping file";
          }

          if (!S_ISREG (sb.st_mode)) {
                  fprintf (stderr, "%s is not a file\n", filename);
                  throw "error mapping file";
          }

          p = (char *) mmap (0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
          if (p == MAP_FAILED) {
                  perror ("mmap");
                  throw "error mapping file";
          }

          if (close (fd) == -1) {
                  perror ("close");
                  throw "error mapping file";
          }

          /*for (len = 0; len < sb.st_size; len++)
                  putchar (p[len]);

          if (munmap (p, sb.st_size) == -1) {
                  perror ("munmap");
                  throw "error mapping file";
          }*/
          mapped ans;
          ans.pnt = p;
          ans.len = sb.st_size;
          return ans;
  }
  
  static long lookLen(char *filename) {
          struct stat sb;
          char *p;
          int fd;

          fd = open (filename, O_RDONLY);
          if (fd == -1) {
                  perror ("open");
                  throw "error mapping file";
          }

          if (fstat (fd, &sb) == -1) {
                  perror ("fstat");
                  throw "error mapping file";
          }

          if (!S_ISREG (sb.st_mode)) {
                  fprintf (stderr, "%s is not a file\n", filename);
                  throw "error mapping file";
          }

          p = (char *) mmap (0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
          if (p == MAP_FAILED) {
                  perror ("mmap");
                  throw "error mapping file";
          }

          if (close (fd) == -1) {
                  perror ("close");
                  throw "error mapping file";
          }
          return sb.st_size;
  }
  
  




        //returns Charlength of letters
        pair<string,char * > readName(char letters[]) {
            string ans;
            int i;
            for (i=0; letters[i]!='\n'; i++) {
                if (letters[i]!='>') {
                    ans.push_back(letters[i]);
                }
            }
            return pair<string,char * >(ans,letters+i);
        }
        
        mapped readData(unsigned char letters[],unsigned char dest[], long len) {
            long j = 0;
            long i;
            for (i=0; letters[i]!='>' && i<len; i++) {
                if (letters[i]=='A'||letters[i]=='C'||letters[i]=='T'||letters[i]=='G') {
                    dest[j>>3] |= ((letters[i]&0x06)<<5)>>(j&0x06);
                    j+=2;
                }
            }
            mapped ans;
            ans.pnt = (char *) letters+i;
            ans.len = j/2;
            return ans;
        }



};

#endif