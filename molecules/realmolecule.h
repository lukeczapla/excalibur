#ifndef _REALMOLECULE_H
#define _REALMOLECULE_H

#include <stdio.h>
#include "../simple.h"

namespace molecules {

typedef int cmapindex[8];

struct Vdihedral {
	integer N;
	real Kchi[MAXDIHEDRAL];
	integer nD[MAXDIHEDRAL];
	real delta[MAXDIHEDRAL];
	real psi0[MAXDIHEDRAL];
	real Kpsi[MAXDIHEDRAL];
};

typedef char groupname[8];

// the parameter set is generally derived for CHARMM-style values

class realmolecule {

public:

	groupname *resname;
	groupname *atomtype;
	groupname *atomname;
	integer *resid;
	
	real scaling14es;
	real scaling14lj;
	real scale;

	real totalmass;
	real mu;  // chemical potential

	integer Natoms;
	integer Natomsbase;
	integer Nbonds;
	integer Nangles;
	integer NuniqueD;
	integer Ndihedrals;

	real *x, *y, *z;
	bool *coordinateassigned;
	real **r;

	real *bond;
	real *angle;
	real *dihedral;

	bool *bondassigned;
	bool *angleassigned;
	bool *dihedralassigned;
	bool *vdwassigned;

	// mass in A.M.U.
	real *mass;
	// partial charge on the atom
	real *charge;
	// sigma is the effective diameter with energy of zero
	real *sigma;
	// epsilon_ij = sqrt(eqsilon_i*epsilon_j)
	real *epsilon;
	real *radii;
	real *gamma, *gamma2;

	integer *bondindex1;
	integer *bondindex2;
	real *b0;
	real *Kb;
	
	integer *angleindex1;
	integer *angleindex2;
	integer *angleindex3;
	real *theta0;
	real *Ktheta;
	real *Kub;
	real *S0;

	integer *dihedralindex1;
	integer *dihedralindex2;
	integer *dihedralindex3;
	integer *dihedralindex4;
	
	Vdihedral *dihedralV;


	bool hascmap;
	integer Ncmap;

	// these are CMAPMAXN arrays
	real ***cmap;
	integer ncmaps;
	integer *cmaplength;
	
	// this is for each phi/psi with a CMAP term
	integer *cmapnumber;
	bool *cmapassigned;

	cmapindex *cindex;

	
	// Monte Carlo sampling options
	integer Nsamplers;
	integer NsampleT;
	integer NsampleR;
	integer NsampleD;
	integer *Tstart;
	integer *Tend;
	integer *Rstart;
	integer *Rend;
	integer *dsdirection;
	
	// Monte Carlo sampling options
	integer *dsindex1;
	integer *dsindex2;
	integer *dsindex3;
	integer *dsindex4;
	
	bool visualize;
	realmolecule();
	realmolecule(const realmolecule& m);
	~realmolecule();
	
	realmolecule& operator+=(const realmolecule& m);
	realmolecule operator+(realmolecule& m1);
	
	void expandcoordinates(int size);
	
	void adddihedral(integer index, real Kch, integer nD, real delta, real psi0, real Kpsi);
	integer checkparameters();

	void drawmolecule(const char *hostname, integer port);
	void drawmolecule(const char *hostname, integer port, real lx, real ly, real lz, bool minimize = false);

	void loadPDBcoordinates(const char *pdbname, bool verbose = false);
	void loadPDBmolecule(const char *pdbname, int N, int startindex = 0, bool verbose = false);

	void loadgamma(FILE *f, bool verbose = false);
	void loadgamma2(FILE *f, bool verbose = false);
	void loadGBradii(FILE *f, bool verbose = false);
	
	void writePDBmolecule(const char *pdbname, int N, bool verbose = false);
	
};

}

#endif

