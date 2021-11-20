#ifndef _DUMPDCD_H
#define _DUMPDCD_H

#include <stdio.h>
#include <stdlib.h>
#include "../simple.h"

class dcdtool {

	FILE *fp;
	char *DCDfilename;
	bool assigned;
	int nframes;
	int ntimestep;
	int natoms;
	
public:

	dcdtool(integer nat);
	~dcdtool();

	void openfile(const char *filename);

	void writeframe(real *x, real *y, real *z);
	void writeDCDheader(const char *remarks, integer ntimestep, float dt);

};


#endif

