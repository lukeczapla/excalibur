#ifndef _MISC_H
#define _MISC_H

#include "../simple.h"


// matches list with ref names to N set of names in forward or reverse order
integer SymmetricMatch(char **names, char **ref, integer *indices, integer N);
integer SpecificSymmetricMatch(char **names, char **ref, integer *indices, integer N);
integer SpecificMatch(char **names, char **ref, integer *indices, integer N);

integer StringToInteger(char *str, int base);
real StringToReal(char *str);

integer IsKeyword(const char *str, integer N, const char **keywords, integer length);
integer IsDelim(char test, char *delimiters, int ndelim);
integer TokenizeString(char *str, char *delimiters, char **pointers, int maxptrs);

real ran1(long *idum);
void UniformSpherePoint(real *rvec, long *seed);


#endif

