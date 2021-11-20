#include <cstdint>
#include <time.h>
#include <string.h>
#include <stdint.h>

#include "dcdtool.h"

#define NFILE_POS 8L
#define NSTEP_POS 20L

// necessary to set SEEK params b/c MPI-2 messes with these settings

#ifndef SEEK_SET
#define SEEK_SET        0
#define SEEK_CUR        1
#define SEEK_END        2
#endif


static inline void fwrite_int32(FILE* fd, uint32_t i) {
	fwrite(&i,sizeof(uint32_t),1,fd);
}


dcdtool::dcdtool(integer nat) {
	nframes = 0;
	natoms = nat;
	assigned = false;
}


dcdtool::~dcdtool() {
	if (assigned) {
		delete [] DCDfilename;
		fclose(fp);
	}
}


void dcdtool::openfile(const char *filename) {
	DCDfilename = new char[128];
	strcpy(DCDfilename, filename);
	fp = fopen(filename,"wb");
	if (fp == NULL) THROWERROR2("Could not open DCD file", filename);
	assigned = true;
}


void dcdtool::writeframe(real *x, real *y, real *z) {
	// vector to save coordinates
	float *xf = new float[natoms];
	float *yf = new float[natoms];
	float *zf = new float[natoms];
	
	for (int i = 0; i < natoms; i++) {
		xf[i] = x[i];
		yf[i] = y[i];
		zf[i] = z[i];
	}

	double dim[6];
	for (int i = 0; i < 6; i++) dim[i] = 90.0;
	uint32_t out_integer = 48;
	fwrite_int32(fp,out_integer);
	fwrite(dim,out_integer,1,fp);
	fwrite_int32(fp,out_integer);

	out_integer = 4*natoms;
	fwrite_int32(fp,out_integer);
	fwrite(xf,out_integer,1,fp);
	fwrite_int32(fp,out_integer);
	fwrite_int32(fp,out_integer);
	fwrite(yf,out_integer,1,fp);
	fwrite_int32(fp,out_integer);
	fwrite_int32(fp,out_integer);
	fwrite(zf,out_integer,1,fp);
	fwrite_int32(fp,out_integer);

	// update NFILE and NSTEP fields in DCD header

	nframes++;
	out_integer = nframes;
	fseek(fp,NFILE_POS,SEEK_SET);
	fwrite_int32(fp,out_integer);
	out_integer = (nframes-1)*ntimestep;
	fseek(fp,NSTEP_POS,SEEK_SET);
	fwrite_int32(fp,out_integer);
	fseek(fp,0,SEEK_END);

	delete [] xf;
	delete [] yf;
	delete [] zf;
	
}



void dcdtool::writeDCDheader(const char *remarks, integer nts, float dt) {
	uint32_t out_integer;
	float out_float;
	char title_string[200];
	time_t cur_time;
	struct tm *tmbuf;
	ntimestep = nts;
	nframes = 0;
	ntimestep = 1;

	out_integer = 84;
	fwrite_int32(fp,out_integer);
	strcpy(title_string,"CORD");
	fwrite(title_string,4,1,fp);
	fwrite_int32(fp,0);					// NFILE = # of snapshots in file
	fwrite_int32(fp,0);					// START = timestep of first snapshot
	fwrite_int32(fp,ntimestep);			// SKIP = interval between snapshots
	fwrite_int32(fp,nframes*ntimestep);	// NSTEP = timestep of last snapshot
	fwrite_int32(fp,0);					// NAMD writes NSTEP or ISTART
	fwrite_int32(fp,0);
	fwrite_int32(fp,0);
	fwrite_int32(fp,0);
	fwrite_int32(fp,0);
	out_float = dt;
	fwrite(&out_float,sizeof(float),1,fp);
	fwrite_int32(fp,1);
	fwrite_int32(fp,0);
	fwrite_int32(fp,0);
	fwrite_int32(fp,0);
	fwrite_int32(fp,0);
	fwrite_int32(fp,0);
	fwrite_int32(fp,0);
	fwrite_int32(fp,0);
	fwrite_int32(fp,0);
	fwrite_int32(fp,27);					// pretend to be Charmm version 27
	fwrite_int32(fp,84);
	fwrite_int32(fp,164);
	fwrite_int32(fp,2);
	strncpy(title_string,remarks,80);
	title_string[79] = '\0';
	fwrite(title_string,80,1,fp);
	cur_time=time(NULL);
	tmbuf=localtime(&cur_time);
	memset(title_string,' ',81);
	strftime(title_string,80,"REMARKS Created by OpenMM interface %d %B,%Y at %H:%M",tmbuf);
	fwrite(title_string,80,1,fp);
	fwrite_int32(fp,164);
	fwrite_int32(fp,4);
	fwrite_int32(fp,natoms);				// number of atoms in each snapshot
	fwrite_int32(fp,4);
	fflush(fp);

}

