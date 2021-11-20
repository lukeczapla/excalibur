#ifndef _SIMPLE_H
#define _SIMPLE_H

/* includes for MPI and OpenMP */

#ifdef _USE_MPI
#include <mpi.h>
#endif
#ifdef _USE_OPENMP
#include <omp.h>
#endif

/* simple type definitions */
typedef double real;
typedef int integer;

/* maximum sizes for MD and GCMC */
// maximum dihedral potentials per dihedral defined
#define MAXDIHEDRAL 10

// simple error functions
#define THROWERROR(msg) {  printf("ERROR: %s(): %s line %d: \"%s\"\n", __func__, __FILE__, __LINE__, msg); exit(-1);  }
#define THROWERROR2(msg,msg2) {  printf("ERROR: %s(): %s line %d: \"%s : %s\"\n", __func__, __FILE__, __LINE__, msg, msg2); exit(-1);  }

#endif

