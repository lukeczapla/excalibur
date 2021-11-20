#include "realmolecule.h"
#include "../util/misc.h"
#include "../util/socket.h"
#include "../util/dcdtool.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


using namespace molecules;


realmolecule::realmolecule() {
	visualize = false;
	Natoms = 0;
	Nbonds = 0;
	Nangles = 0;
	NuniqueD = 0;
	Ndihedrals = 0;
	Ncmap = 0;
	ncmaps = 0;
	totalmass = 0.0;
	cmap = NULL; hascmap = false;
	r = new real*[3];
}


realmolecule::realmolecule(const realmolecule& m) : Natomsbase(m.Natomsbase), Natoms(m.Natoms), Nbonds(m.Nbonds), Nangles(m.Nangles),
	Ndihedrals(m.Ndihedrals), NuniqueD(m.NuniqueD), Ncmap(m.Ncmap), hascmap(m.hascmap), ncmaps(m.ncmaps),
	totalmass(m.totalmass), mu(m.mu), scaling14es(m.scaling14es), scaling14lj(m.scaling14lj), scale(m.scale) {

	x = new real[Natoms];
	y = new real[Natoms];
	z = new real[Natoms];
	r = new real*[3];
	r[0] = x;
	r[1] = y;
	r[2] = z;
	charge = new real[Natoms];
	sigma = new real[Natoms];
	epsilon = new real[Natoms];
	mass = new real[Natoms];
	coordinateassigned = new bool[Natoms];
	vdwassigned = new bool[Natoms];

	for (int i = 0; i < Natoms; i++) {
		x[i] = m.x[i];
		y[i] = m.y[i];
		z[i] = m.z[i];
		charge[i] = m.charge[i];
		sigma[i] = m.sigma[i];
		epsilon[i] = m.epsilon[i];
		mass[i] = m.mass[i];
		coordinateassigned[i] = m.coordinateassigned[i];
		vdwassigned[i] = m.vdwassigned[i];
	}
	
	bond = new real[Nbonds];
	bondindex1 = new integer[Nbonds];
	bondindex2 = new integer[Nbonds];
	b0 = new real[Nbonds];
	Kb = new real[Nbonds];
	bondassigned = new bool[Nbonds];
	for (int i = 0; i < Nbonds; i++) {
		bond[i] = m.bond[i];
		bondindex1[i] = m.bondindex1[i];
		bondindex2[i] = m.bondindex2[i];
		b0[i] = m.b0[i];
		Kb[i] = m.Kb[i];
		bondassigned[i] = m.bondassigned[i];
	}
	
	angle = new real[Nangles];
	angleindex1 = new integer[Nangles];
	angleindex2 = new integer[Nangles];
	angleindex3 = new integer[Nangles];
	Ktheta = new real[Nangles];
	theta0 = new real[Nangles];
	Kub = new real[Nangles];
	S0 = new real[Nangles];
	angleassigned = new bool[Nangles];
	for (int i = 0; i < Nangles; i++) {
		angle[i] = m.angle[i];
		angleindex1[i] = m.angleindex1[i];
		angleindex2[i] = m.angleindex2[i];
		angleindex3[i] = m.angleindex3[i];
		Ktheta[i] = m.Ktheta[i];
		theta0[i] = m.theta0[i];
		Kub[i] = m.Kub[i];
		S0[i] = m.S0[i];
		angleassigned[i] = m.angleassigned[i];
	}

	dihedral = new real[NuniqueD];
	dihedralindex1 = new integer[NuniqueD];
	dihedralindex2 = new integer[NuniqueD];
	dihedralindex3 = new integer[NuniqueD];
	dihedralindex4 = new integer[NuniqueD];
	dihedralV = new Vdihedral[NuniqueD];
	dihedralassigned = new bool[NuniqueD];
	for (int i = 0; i < NuniqueD; i++) {
		dihedral[i] = m.dihedral[i];
		dihedralindex1[i] = m.dihedralindex1[i];
		dihedralindex2[i] = m.dihedralindex2[i];
		dihedralindex3[i] = m.dihedralindex3[i];
		dihedralindex4[i] = m.dihedralindex4[i];
		dihedralassigned[i] = m.dihedralassigned[i];
		dihedralV[i].N = m.dihedralV[i].N;
		for (int j = 0; j < m.dihedralV[i].N; j++) {
			dihedralV[i].Kchi[j] = m.dihedralV[i].Kchi[j];
			dihedralV[i].nD[j] = m.dihedralV[i].nD[j];
			dihedralV[i].delta[j] = m.dihedralV[i].delta[j];
			dihedralV[i].psi0[j] = m.dihedralV[i].psi0[j];
			dihedralV[i].Kpsi[j] = m.dihedralV[i].Kpsi[j];
		}
		
	}
	if (hascmap) {
		cmapnumber = new integer[Ncmap];
		cmapassigned = new bool[Ncmap];
		cindex = new cmapindex[Ncmap];
		for (int i = 0; i < Ncmap; i++) {
			cmapnumber[i] = m.cmapnumber[i];
			cmapassigned[i] = m.cmapassigned[i];
			for (int j = 0; j < 8; j++) cindex[i][j] = m.cindex[i][j];
		}
		cmaplength = new integer[ncmaps];
		cmap = new real**[ncmaps];
		for (int i = 0; i < ncmaps; i++) {
			cmaplength[i] = m.cmaplength[i];
			cmap[i] = new real*[cmaplength[i]];
			for (int I = 0; I < cmaplength[i]; I++) {
				cmap[i][I] = new real[cmaplength[i]];
				for (int J = 0; J < cmaplength[i]; J++) {
					cmap[i][I][J] = m.cmap[i][I][J];
				}
			}
		}
	}
	
	resname = new groupname[Natomsbase];
	atomname = new groupname[Natomsbase];
	atomtype = new groupname[Natomsbase];
	resid = new integer[Natomsbase];
	for (int i = 0; i < Natomsbase; i++) {
		strcpy(resname[i], m.resname[i]);
		strcpy(atomname[i], m.atomname[i]);
		strcpy(atomtype[i], m.atomtype[i]);
		resid[i] = m.resid[i];
	}
	
}


realmolecule::~realmolecule() {
	delete [] r;
	if (Natoms > 0) {
		delete [] x;
		delete [] y;
		delete [] z;
		delete [] coordinateassigned;
		delete [] charge;
		delete [] sigma;
		delete [] epsilon;
		delete [] vdwassigned;
	}
	if (Nbonds > 0) {
		delete [] bondindex1;
		delete [] bondindex2;
		delete [] Kb;
		delete [] b0;
		delete [] bond;
		delete [] bondassigned;
	}
	if (Nangles > 0) {
		delete [] angleindex1;
		delete [] angleindex2;
		delete [] angleindex3;
		delete [] theta0;
		delete [] Ktheta;
		delete [] Kub;
		delete [] S0;
		delete [] angle;
		delete [] angleassigned;
	}
	if (Ndihedrals > 0) {
		delete [] dihedral;
		delete [] dihedralindex1;
		delete [] dihedralindex2;
		delete [] dihedralindex3;
		delete [] dihedralindex4;
		delete [] dihedralV;
	}
}



// assumes it is the same CMAP if both have it!!!
// TO BE FINISHED
realmolecule& realmolecule::operator+=(const realmolecule& m) {

	totalmass += m.totalmass;
	Natoms += m.Natoms;
	Nbonds += m.Nbonds;
	Nangles += m.Nangles;
	NuniqueD += m.Nangles;
	Ndihedrals += m.Ndihedrals;
	Ncmap += m.Ncmap;
	
	if (!hascmap && m.hascmap) {
		// fix up CMAP
	}



	x = new real[Natoms];
	y = new real[Natoms];
	z = new real[Natoms];

	charge = new real[Natoms];
	sigma = new real[Natoms];
	epsilon = new real[Natoms];
	mass = new real[Natoms];
	coordinateassigned = new bool[Natoms];
	vdwassigned = new bool[Natoms];

	
	for (int i = 0; i < Natoms; i++) {
		x[i] = m.x[i];
		y[i] = m.y[i];
		z[i] = m.z[i];
		charge[i] = m.charge[i];
		sigma[i] = m.sigma[i];
		epsilon[i] = m.epsilon[i];
		mass[i] = m.mass[i];
		coordinateassigned[i] = m.coordinateassigned[i];
		vdwassigned[i] = m.vdwassigned[i];
	}
	
	bond = new real[Nbonds];
	bondindex1 = new integer[Nbonds];
	bondindex2 = new integer[Nbonds];
	b0 = new real[Nbonds];
	Kb = new real[Nbonds];
	bondassigned = new bool[Nbonds];
	
	for (int i = 0; i < Nbonds; i++) {
		bond[i] = m.bond[i];
		bondindex1[i] = m.bondindex1[i];
		bondindex2[i] = m.bondindex2[i];
		b0[i] = m.b0[i];
		Kb[i] = m.Kb[i];
		bondassigned[i] = m.bondassigned[i];
	}
	
	angle = new real[Nangles];
	angleindex1 = new integer[Nangles];
	angleindex2 = new integer[Nangles];
	angleindex3 = new integer[Nangles];
	Ktheta = new real[Nangles];
	theta0 = new real[Nangles];
	angleassigned = new bool[Nangles];

	for (int i = 0; i < Nangles; i++) {
		angle[i] = m.angle[i];
		angleindex1[i] = m.angleindex1[i];
		angleindex2[i] = m.angleindex2[i];
		angleindex3[i] = m.angleindex3[i];
		Ktheta[i] = m.Ktheta[i];
		theta0[i] = m.theta0[i];
		Kub[i] = m.Kub[i];
		S0[i] = m.S0[i];
		angleassigned[i] = m.angleassigned[i];
	}

	dihedral = new real[NuniqueD];
	dihedralV = new Vdihedral[NuniqueD];
	
	for (int i = 0; i < NuniqueD; i++) {
		dihedral[i] = m.dihedral[i];
		dihedralindex1[i] = m.dihedralindex1[i];
		dihedralindex2[i] = m.dihedralindex2[i];
		dihedralindex3[i] = m.dihedralindex3[i];
		dihedralindex4[i] = m.dihedralindex4[i];
		dihedralV[i].N = m.dihedralV[i].N;
		for (int j = 0; j < dihedralV[i].N; j++) {
			dihedralV[i].Kchi[j] = m.dihedralV[i].Kchi[j];
			dihedralV[i].nD[j] = m.dihedralV[i].nD[j];
			dihedralV[i].delta[j] = m.dihedralV[i].delta[j];
			dihedralV[i].psi0[j] = m.dihedralV[i].psi0[j];
			dihedralV[i].Kpsi[j] = m.dihedralV[i].Kpsi[j];
		}
		
	}
	
	if (hascmap) {
		cmapnumber = new integer[Ncmap];
		cmapassigned = new bool[Ncmap];
		cindex = new cmapindex[Ncmap];
		for (int i = 0; i < Ncmap; i++) {
			cmapnumber[i] = m.cmapnumber[i];
			cmapassigned[i] = m.cmapassigned[i];
			for (int j = 0; j < 8; j++) cindex[i][j] = m.cindex[i][j];
		}
		cmaplength = new integer[ncmaps];
		cmap = new real**[ncmaps];
		for (int i = 0; i < ncmaps; i++) {
			cmaplength[i] = m.cmaplength[i];
			cmap[i] = new real*[cmaplength[i]];
			for (int I = 0; I < cmaplength[i]; I++) {
				cmap[i][I] = new real[cmaplength[i]];
				for (int J = 0; J < cmaplength[i]; I++) {
					cmap[i][I][J] = m.cmap[i][I][J];
				}
			}
		}
	}

	return *this;
	
}


realmolecule realmolecule::operator+(realmolecule& m1) {
	realmolecule mcopy(*this);
	return mcopy += m1;
}


void realmolecule::expandcoordinates(int size) {
	if (size < Natoms) THROWERROR("Not enough space to allocate molecule in expandcoordinates()");
	real *xnew = new real[size];
	real *ynew = new real[size];
	real *znew = new real[size];
	real *mnew = new real[size];
	real *snew = new real[size];
	for (int i = 0; i < Natoms; i++) {
		xnew[i] = x[i];
		ynew[i] = y[i];
		znew[i] = z[i];
		mnew[i] = mass[i];
		snew[i] = sigma[i];
	}
	delete [] x;
	delete [] y;
	delete [] z;
	delete [] mass;
	delete [] sigma;
	x = new real[size]; r[0] = x;
	y = new real[size]; r[1] = y;
	z = new real[size]; r[2] = z;
	mass = new real[size];
	sigma = new real[size];
	for (int i = 0; i < Natoms; i++) {
		x[i] = xnew[i];
		y[i] = ynew[i];
		z[i] = znew[i];
		mass[i] = mnew[i];
		sigma[i] = snew[i];
	}
	for (int i = Natoms; i < size; i++) {
		x[i] = x[i % Natoms];
		y[i] = y[i % Natoms];
		z[i] = z[i % Natoms];
		mass[i] = mass[i % Natoms];
		sigma[i] = sigma[i % Natoms];

	}
	delete [] xnew;
	delete [] ynew;
	delete [] znew;
	delete [] mnew;
	delete [] snew;
}


void realmolecule::adddihedral(integer index, real Kchi, integer nD, real delta, real psi0, real Kpsi) {
	if (index >= NuniqueD) THROWERROR("Out of range dihedral value!");
	dihedralassigned[index] = true;

	dihedralV[index].N++;
	if (dihedralV[index].N > 10) THROWERROR("Too many dihedral parameters for single dihedral, greater than 10");
	dihedralV[index].Kchi[dihedralV[index].N-1] = Kchi;
	dihedralV[index].nD[dihedralV[index].N-1] = nD;
	dihedralV[index].delta[dihedralV[index].N-1] = delta;
	dihedralV[index].psi0[dihedralV[index].N-1] = psi0;
	dihedralV[index].Kpsi[dihedralV[index].N-1] = Kpsi;
	return;
}


#define MISSINGPARAMETER(msg,num) { printf("Missing %s parameter for number %d\n", msg, num); exit(-1); }

integer realmolecule::checkparameters() {
	for (int i = 0; i < Natoms; i++) {
		if (!coordinateassigned[i]) MISSINGPARAMETER("coordinates", i);
		if (!vdwassigned[i]) MISSINGPARAMETER("VDW", i); 
	}
	for (int i = 0; i < Nbonds; i++) if (!bondassigned[i]) MISSINGPARAMETER("bond", i);
	for (int i = 0; i < Nangles; i++) if (!angleassigned[i]) MISSINGPARAMETER("angle", i);
	for (int i = 0; i < NuniqueD; i++) if (!dihedralassigned[i]) MISSINGPARAMETER("dihedral", i);
	return 1;
}


void realmolecule::drawmolecule(const char *hostname, integer port, real lx, real ly, real lz, bool minimize) {
	if (visualize == false) return;
	int size;
	double *cell = new double[3];
	double *radii = new double[Natoms];
	int *indices = new int[Natoms];
	double *coords = new double[3*Natoms];
	cell[0] = lx;
	cell[1] = ly;
	cell[2] = lz;
	size = Natoms;
	for (int i = 0; i < Natoms; i++) {
		if (minimize) {
			coords[3*i] = x[i];
			coords[3*i+1] = y[i];
			coords[3*i+2] = z[i];
		} else {
			coords[3*i] = x[i]-lx/2.0;
			coords[3*i+1] = y[i]-ly/2.0;
			coords[3*i+2] = z[i]-lz/2.0;		
		}
		if (mass[i] < 2.0) indices[i] = 0;
		else if (fabs(mass[i]-12.0) < 1.0) indices[i] = 1;
		else if (fabs(mass[i]-14.5) < 1.0) indices[i] = 2;
		else if (fabs(mass[i]-16.0) < 0.2) indices[i] = 3;
		else indices[i] = 4;
		radii[i] = sigma[i]/2.0;
		if (radii[i] < 0.3) radii[i] += 0.8;
	}
	jsocket *js = GetJSocket(hostname, port, 0);
	bool end = false;
	if (js != NULL) {
		if (WriteJSocket(js->socket, cell, 3*sizeof(double)) == -1) end = true;
		if (WriteJSocket(js->socket, &size, sizeof(int)) == -1) end = true;
		if (WriteJSocket(js->socket, radii, size*sizeof(double)) == -1) end = true;
		if (WriteJSocket(js->socket, indices, size*sizeof(int)) == -1) end = true;
		if (WriteJSocket(js->socket, coords, 3*size*sizeof(double)) == -1) end = true;
		FreeJSocket(js);
	}
	delete [] cell;
	delete [] radii;
	delete [] indices;
	delete [] coords;
}





void realmolecule::drawmolecule(const char *hostname, integer port) {
	if (visualize == false) return;
	int size;
	double *cell = new double[3];
	double *radii = new double[Natoms];
	int *indices = new int[Natoms];
	double *coords = new double[3*Natoms];
	cell[0] = 0.0;
	cell[1] = 0.0;
	cell[2] = 0.0;
	size = Natoms;
	for (int i = 0; i < Natoms; i++) {
		coords[3*i] = x[i];
		coords[3*i+1] = y[i];
		coords[3*i+2] = z[i];
		if (mass[i] < 2.0) indices[i] = 0;
		else if (fabs(mass[i]-12.0) < 1.0) indices[i] = 1;
		else if (fabs(mass[i]-14.0) < 1.0) indices[i] = 2;
		else if (fabs(mass[i]-16.0) < 1.0) indices[i] = 3;
		else indices[i] = 4;
		radii[i] = sigma[i]/2.0;
		if (radii[i] < 0.3) radii[i] += 0.8;
	}
	jsocket *js = GetJSocket(hostname, port, 0);
	bool end = false;
	if (js != NULL) {
		if (WriteJSocket(js->socket, cell, 3*sizeof(double)) == -1) end = true;
		if (!end) if (WriteJSocket(js->socket, &size, sizeof(int)) == -1) end = true;
		if (!end) if (WriteJSocket(js->socket, radii, size*sizeof(double)) == -1) end = true;
		if (!end) if (WriteJSocket(js->socket, indices, size*sizeof(int)) == -1) end = true;
		if (!end) if (WriteJSocket(js->socket, coords, 3*size*sizeof(double)) == -1) end = true;
		FreeJSocket(js);
	}
	delete [] cell;
	delete [] radii;
	delete [] indices;
	delete [] coords;
}


void realmolecule::loadPDBcoordinates(const char *pdbname, bool verbose) {

	char *p;
	char s[256];
	char word[32];
	real xx, yy, zz, occ, beta;
	char atom[8];
	char atomnum[8];
	char atomtype[8];
	char chain[8];
	char residue[8];
	char resn[8];
	char none[8];
	if (Natoms == 0) THROWERROR2("Attempt to read PDB file for empty molecule", pdbname);
	FILE *f = fopen(pdbname, "r");
	if (f == NULL) THROWERROR2("Empty PDB file", pdbname);
	int i = 0;
	do {
		
		p = fgets(s, 1024, f);

      	if (p == NULL) THROWERROR2("Invalid data in PDB file", pdbname);
		sscanf(s, "%s", word);
		if (strcmp(word, "CRYST1") == 0) {
			if (sscanf(s, "%s %lf %lf %lf", word, &xx, &yy, &zz) == 4) {
				if (verbose) printf("Read box dimensions %7.3lf %7.3lf %7.3lf from PDB file\n", xx, yy, zz);
			} else {
				if (verbose) printf("Error in CRYST1 line format in PDB file %s, continuing...\n", pdbname);
			}
		}
		else if (strcmp(word, "ATOM") == 0) {
			sscanf(s, "%6c%5c%2c%4c%4c%1c%4c%4c%8lf%8lf%8lf%6lf%6lf", atom, atomnum, none, atomtype, residue, chain, resn, none, &xx, &yy, &zz, &occ, &beta);
			x[i] = xx; y[i] = yy; z[i] = zz; 
			coordinateassigned[i] = true;
			if (verbose) printf("Assigned coordinates to molecule %d\n", i);
			i++;
		}
	} while (i < Natoms);
	
	return;
}


void realmolecule::loadPDBmolecule(const char *pdbname, int N, int index, bool verbose) {

	char *p;
	char s[256];
	char word[32];
	real xx, yy, zz, occ, beta;
	char atom[8];
	char atomnum[8];
	char atomtype[8];
	char chain[8];
	char residue[8];
	char resn[8];
	char none[8];
	if (Natoms == 0) THROWERROR2("Attempt to read PDB file for empty molecule", pdbname);
	FILE *f = fopen(pdbname, "r");
	if (f == NULL) THROWERROR2("Empty PDB file", pdbname);
	int i = 0;
	do {
		
		p = fgets(s, 1024, f);

      	if (p == NULL) THROWERROR2("Invalid data in PDB file", pdbname);
		sscanf(s, "%s", word);
		if (strcmp(word, "CRYST1") == 0) {
			if (sscanf(s, "%s %lf %lf %lf", word, &xx, &yy, &zz) == 4) {
				if (verbose) printf("Read box dimensions %7.3lf %7.3lf %7.3lf from PDB file\n", xx, yy, zz);
			} else {
				if (verbose) printf("Error in CRYST1 line format in PDB file %s, continuing...\n", pdbname);
			}
		}
		else if (strcmp(word, "ATOM") == 0) {
			sscanf(s, "%6c%5c%2c%4c%4c%1c%4c%4c%8lf%8lf%8lf%6lf%6lf", atom, atomnum, none, atomtype, residue, chain, resn, none, &xx, &yy, &zz, &occ, &beta);
			x[index+i] = xx; y[index+i] = yy; z[index+i] = zz; 
			//coordinateassigned[i] = true;
			if (verbose) printf("Assigned coordinates to molecule %d\n", i);
			i++;
		}
	} while (i < N);
	fclose(f);
	return;
}


void realmolecule::loadgamma(FILE *f, bool verbose) {
	gamma = new real[Natoms];
	for (int i = 0; i < Natoms; i++) {
		fscanf(f, "%lf", &gamma[i]);
		if (verbose) printf("GBVI gamma for atom %d is %lf kcal/mol\n", i, gamma[i]);
	}
}


void realmolecule::loadgamma2(FILE *f, bool verbose) {
        gamma2 = new real[Natoms];
        for (int i = 0; i < Natoms; i++) {
                fscanf(f, "%lf", &gamma2[i]);
                if (verbose) printf("GBVI gamma2 for atom %d is %lf kcal/mol\n", i, gamma2[i]);
        }

}


void realmolecule::loadGBradii(FILE *f, bool verbose) {
	radii = new real[Natoms];
	for (int i = 0; i < Natoms; i++) {
		fscanf(f, "%lf", &radii[i]);
		if (verbose) printf("GB radius for atom %d is %lf Å\n", i, radii[i]);
	}
}


void realmolecule::writePDBmolecule(const char *pdbname, int N, bool verbose) {
	FILE *f = fopen(pdbname, "w");
	printf("%d\n", Natomsbase);
	for (int i = 0; i < N; i++) {
		fprintf(f, "ATOM  %5d %4s %3s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", i+1, atomname[i % Natomsbase], resname[i % Natomsbase], resid[i % Natomsbase], x[i], y[i], z[i], 1.0, 1.0);
	}
	fclose(f);
}
