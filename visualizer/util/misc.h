#ifndef _MISC_H
#define _MISC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define THROWERROR(msg) {  printf("ERROR: %s(): %s line %d: \"%s\"\n", __func__, __FILE__, __LINE__, msg); exit(-1);  }


// Convert null-terminated character strings to various data types. Better error checking than standard routines!
int StringToInt(char *str, int base, char *file, int lineno);
double StringToDouble(char *str, char *file, int lineno);


// Some handy string tokenisation routines to make things more pleasant than using the c standard library functions.
int IsDelim(char test, char *delimiters, int ndelim);
int TokenizeString(char *str, char *delimiters, char **pointers, int maxptrs);


// The ran1() function from Numerical Recipes in C, and  routine to get a unit vector uniformly distributed on a sphere
double ran1(long *idum);
void UniformSpherePoint(double *rvec, long *seed);


#endif

