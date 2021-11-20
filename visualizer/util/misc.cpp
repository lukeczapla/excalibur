
#include "misc.h"


int StringToInt(char *str, int base, char *file, int lineno) {
	char *endptr;
	char message[1024];
	int returned;
	
	returned = strtol(str, &endptr, 10);
	if(*endptr != '\0')
	{
		sprintf( message, "Unable to convert '%s' (%s, line %d) into an integer of base %d", str, file, lineno, base );
		return 0;
	}
	return returned;
}


double StringToDouble(char *str, char *file, int lineno) {
	char *endptr;
	char message[1024];
	double returned;

	returned = strtod(str, &endptr);
	if(*endptr != '\0')
	{
		sprintf(message, "Unable to convert '%s' (%s, line %d) into a double", str, file, lineno);
		return 0.0;
	}
	return returned;
}


int IsDelim(char test, char *delimiters, int ndelim) {
	int i;
	for (i = 0; i < ndelim; i++)
		if (test == delimiters[i]) return 1;
	return 0;
}


int TokenizeString(char *str, char *delimiters, char **pointers, int maxptrs) {
	int ntoks, len, ndelim, i;
	
	len = strlen(str);
	ndelim = strlen(delimiters);
	
	ntoks = 0;
	i = 0;
	while (1)
	{
		/* skip and zero any leading delimiters - should only matter on first iteration. */
		while (IsDelim(str[i],delimiters,ndelim) && i < len) { 
			str[i] = '\0'; 
			i++; 
		}
		if (i >= len) return ntoks;
		/* store start of string */
		pointers[ntoks] = &str[i];
		/* skip until next delimiter or end of string found */
		while (!IsDelim(str[i],delimiters,ndelim) && i < len) i++;
		ntoks++;
		if (ntoks >= maxptrs) return ntoks;
		if (i >= len) return ntoks;
	}
}



#define IA 16807 
#define IM 2147483647 
#define AM (1.0/IM) 
#define IQ 127773 
#define IR 2836 
#define NTAB 32 
#define NDIV (1+(IM-1)/NTAB) 
#define EPS 1.2e-7 
#define RNMX (1.0-EPS)

double ran1(long *idum) {
	int j; 
	long k; 
	static long iy = 0; 
	static long iv[NTAB]; 
	double temp;
	
	if( *idum <= 0 || !iy )
	{
		/* initialize */
		if (-(*idum) < 1) *idum=1; /* prevent idum == 0 */
		else *idum =-(*idum);

		for( j = NTAB+7; j >= 0; j--)
		{
			/* Load the shuffle table ( after 8 warm-ups ) */
			k = (*idum) / IQ; 
			*idum = IA*(*idum-k*IQ)-IR*k; 
			if (*idum < 0) *idum += IM; 
			if (j < NTAB) iv[j] = *idum; 
		} 
		iy = iv[0]; 
	}
	
	k = (*idum) / IQ; /* Start here when not initializing */
	*idum = IA*(*idum-k*IQ)-IR*k; /* Compute idum = (IA*idum) % IM without overflows by Schrage's method */
	if (*idum < 0) *idum += IM; 
	j = iy / NDIV; /* Will be in the range 0..NTAB-1. */
	iy = iv[j]; /* Output previously stored value and refill the shuffle table */
	iv[j] = *idum; 
	if ((temp = AM*iy) > RNMX) return RNMX; /* Because users don’t expect end point values. */
	else return temp; 
}

#undef IA 
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


void UniformSpherePoint(double *rvec, long *seed) {
	double u, theta;
	
	u = (ran1(seed)-0.5) * 2.0;		/* u \elem [-1,1] */
	u = sqrt(1.0 - u*u);
	theta = ran1(seed) * 2.0*M_PI;	/* theta \elem [0,2pi) */
	rvec[0] = u * cos(theta);
	rvec[1] = u * sin(theta);
	rvec[2] = u;
}


