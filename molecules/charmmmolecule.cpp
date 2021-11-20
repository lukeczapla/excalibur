#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "charmmmolecule.h"
#include "../util/misc.h"

#define TWOROOTONESIXTH 1.12246205


using namespace molecules;
using namespace std;


charmmmolecule::charmmmolecule() : realmolecule() { 
	scaling14es = 1.0;
	scaling14lj = 1.0;
}


charmmmolecule::charmmmolecule(const char *psfname) : realmolecule() {
	if (psfname != NULL) loadmoleculePSF(psfname, false);
}


charmmmolecule::charmmmolecule(const char *psfname, const char *parname) : realmolecule() {
	if (psfname != NULL) loadmoleculePSF(psfname, true);
	if (parname != NULL) loadmoleculePRM(parname, true);
}


charmmmolecule::charmmmolecule(const charmmmolecule& m) : realmolecule(m), Nproper(m.Nproper), Nimproper(m.Nimproper) {
		resname = new groupname[Natoms];
		atomtype = new groupname[Natoms];
		atomname = new groupname[Natoms];
		segname = new groupname[Natoms];
		resid = new integer[Natoms];
		for (int i = 0; i < Natoms; i++) {
			strcpy(resname[i], m.resname[i]);
			strcpy(atomtype[i], m.atomtype[i]);
			strcpy(atomname[i], m.atomname[i]);
			strcpy(segname[i], m.segname[i]);
			resid[i] = m.resid[i];
		}
		isspecificD = new bool[NuniqueD];
		for (int i = 0; i < NuniqueD; i++) isspecificD[i] = true;
}


charmmmolecule::~charmmmolecule() {
	if (Natoms != 0) {
		delete [] resname;
		delete [] atomname;
		delete [] atomtype;
		delete [] segname;
		delete [] resid;
		if (hascmap && (cmap != NULL)) delete [] cmap;
	}
}


void charmmmolecule::loadmoleculePSF(const char *psfname, bool verbose) {
	FILE* f = fopen(psfname, "r");
	if (f == NULL) THROWERROR2("Could not open PSF file", psfname);
	char* p;
	char s[1024];

	int count;

	char* space = new char[3];
	strcpy(space, " \t");
	char** strings = new char*[32];
	for (int i = 0; i < 32; i++) strings[i] = new char[32];

	p = fgets(s, 1024, f);
	count = TokenizeString(s, space, strings, 32);
	if ((count == 0) || (strcasecmp(strings[0], "PSF") != 0)) THROWERROR2("Invalid PSF file", psfname);

	for (int i = 1; (i < count) && !hascmap; i++) if (strcasecmp(strings[i], "CMAP") == 0) hascmap = true;
	
	bool done = false;
	do {
		p = fgets(s, 1024, f);
		if (p == NULL) THROWERROR2("Could not find atoms in PSF file", psfname);
		count = TokenizeString(s, space, strings, 32);
		if ((count >= 2) && (strncmp(strings[1], "!NATOM", 6) == 0)) done = true;
	} while (!done);

	Natoms = StringToInteger(strings[0], 10);
	if (verbose) printf("%d atoms to read from PSF file\n", Natoms);
	resname = new groupname[Natoms];
	atomname = new groupname[Natoms];
	atomtype = new groupname[Natoms];
	segname = new groupname[Natoms];
	resid = new integer[Natoms];
	x = new real[Natoms];
	y = new real[Natoms];
	z = new real[Natoms];
	coordinateassigned = new bool[Natoms];
	vdwassigned = new bool[Natoms];
	scale = 1.0;
	scaling14es = 1.0;
	scaling14lj = 1.0;

	for (int i = 0; i < Natoms; i++) {
		x[i] = 0.0;
		y[i] = 0.0;
		z[i] = 0.0;
		coordinateassigned[i] = false;
		vdwassigned[i] = false;
	}

	mass = new real[Natoms];
	charge = new real[Natoms];
	
	sigma = new real[Natoms];
	epsilon = new real[Natoms];
	totalmass = 0.0;

	for (int i = 0; i < Natoms; i++) {
		p = fgets(s, 1024, f);
		if (p == NULL) THROWERROR2("Could not read all atoms in PSF file", psfname);
		count = TokenizeString(s, space, strings, 32);

		if (count < 8) THROWERROR2("Not enough parameters for atom in PSF file", psfname);

		strncpy(segname[i], strings[1], 7);
		resid[i] = StringToInteger(strings[2], 10);

		strncpy(resname[i], strings[3], 7);
		strncpy(atomname[i], strings[4], 7);
		strncpy(atomtype[i], strings[5], 7);
		charge[i] = StringToReal(strings[6]);
		mass[i] = StringToReal(strings[7]);
		totalmass += mass[i];
		if (verbose) printf("Atom %d : %s %s %s, %lf %lf\n", i+1, atomname[i], atomtype[i], resname[i], charge[i], mass[i]);
	}

	done = false;
	do {
		p = fgets(s, 1024, f);
		if (p == NULL) THROWERROR2("Could not find bonds in PSF file", psfname);
		count = TokenizeString(s, space, strings, 32);
		if ((count >= 2) && (strncmp(strings[1], "!NBOND", 6) == 0)) done = true;
	} while (!done);
	
	Nbonds = StringToInteger(strings[0], 10);
	
	bondindex1 = new integer[Nbonds];
	bondindex2 = new integer[Nbonds];
	Kb = new real[Nbonds];
	b0 = new real[Nbonds];
	bond = new real[Nbonds];
	bondassigned = new bool[Nbonds];

	for (int i = 0; i < Nbonds; i++) {
		Kb[i] = 0.0;
		b0[i] = 1.0;
		bond[i] = 1.5;
		bondassigned[i] = false;
	}
	
	for (int i = 0; i < Nbonds;) {
		p = fgets(s, 1024, f);
		if (p == NULL) THROWERROR2("Could not find all bonds in PSF file", psfname);
		count = TokenizeString(s, space, strings, 32);
		if ((count % 2) != 0) THROWERROR2("Uneven bond spacing in PSF file", psfname);
		for (int j = 0; j < count; j += 2) {
			bondindex1[i] = StringToReal(strings[j]);
			bondindex2[i] = StringToReal(strings[j+1]);
			if (verbose) printf("Bond %d: %d %d\n", i+1, bondindex1[i], bondindex2[i]);
			i++;
		}

	}
	
	done = false;
	do {
		p = fgets(s, 1024, f);
		if (p == NULL) THROWERROR2("Could not find angles in PSF file", psfname);
		count = TokenizeString(s, space, strings, 32);
		if ((count >= 2) && (strncmp(strings[1], "!NTHET", 6) == 0)) done = true;
	} while (!done);

	Nangles = StringToInteger(strings[0], 10);
	angleindex1 = new integer[Nangles];
	angleindex2 = new integer[Nangles];
	angleindex3 = new integer[Nangles];
	angle = new real[Nangles];
	theta0 = new real[Nangles];
	Ktheta = new real[Nangles];
	Kub = new real[Nangles];
	S0 = new real[Nangles];
	angleassigned = new bool[Nangles];
	for (int i = 0; i < Nangles; i++) {
		theta0[i] = 120.0;
		Ktheta[i] = 1.0;
		Kub[i] = 0.0;
		S0[i] = 0.0;
		angleassigned[i] = false;
	}

	for (int i = 0; i < Nangles;) {
		p = fgets(s, 1024, f);
		if (p == NULL) THROWERROR2("Could not find all angles in PSF file", psfname);
		count = TokenizeString(s, space, strings, 32);
		if ((count % 3) != 0) THROWERROR2("Uneven angle spacing in PSF file", psfname);
		for (int j = 0; j < count; j += 3) {
			angleindex1[i] = StringToReal(strings[j]);
			angleindex2[i] = StringToReal(strings[j+1]);
			angleindex3[i] = StringToReal(strings[j+2]);
			if (verbose) printf("Angle %d: %d %d %d\n", i+1, angleindex1[i], angleindex2[i], angleindex3[i]);
			i++;
		}
	}

	done = false;
	do {
		p = fgets(s, 1024, f);
		if (p == NULL) THROWERROR2("Could not find dihedrals in PSF file", psfname);
		count = TokenizeString(s, space, strings, 32);
		if ((count >= 2) && (strncmp(strings[1], "!NPHI", 5) == 0)) done = true;
	} while (!done);

	NuniqueD = StringToInteger(strings[0], 10);
	Ndihedrals = NuniqueD;
	Nproper = NuniqueD;
	integer *tempindex1 = new integer[NuniqueD];
	integer *tempindex2 = new integer[NuniqueD];
	integer *tempindex3 = new integer[NuniqueD];
	integer *tempindex4 = new integer[NuniqueD];

	for (int i = 0; i < Ndihedrals;) {
		p = fgets(s, 1024, f);
		if (p == NULL) THROWERROR2("Could not find all dihedrals in PSF file", psfname);
		count = TokenizeString(s, space, strings, 32);
		if ((count % 4) != 0) THROWERROR2("Uneven dihedral spacing in PSF file", psfname);
		for (int j = 0; j < count; j += 4) {
			tempindex1[i] = StringToReal(strings[j]);
			tempindex2[i] = StringToReal(strings[j+1]);
			tempindex3[i] = StringToReal(strings[j+2]);
			tempindex4[i] = StringToReal(strings[j+3]);
			i++;
		}
	}

	done = false;
	do {
		p = fgets(s, 1024, f);
		if (p == NULL) THROWERROR2("Could not find impropers in PSF file", psfname);
		count = TokenizeString(s, space, strings, 32);
		if ((count >= 2) && (strncmp(strings[1], "!NIMPH", 6) == 0)) done = true;
	} while (!done);

	Nimproper = StringToInteger(strings[0], 10);
	NuniqueD += Nimproper;
	Ndihedrals += Nimproper;
	dihedral = new real[NuniqueD];
	dihedralindex1 = new integer[NuniqueD];
	dihedralindex2 = new integer[NuniqueD];
	dihedralindex3 = new integer[NuniqueD];
	dihedralindex4 = new integer[NuniqueD];
	dihedralV = new Vdihedral[NuniqueD];
	isspecificD = new bool[NuniqueD];
	dihedralassigned = new bool[NuniqueD];

	for (int i = 0; i < NuniqueD; i++) {
		dihedralV[i].N = 0;
		dihedralassigned[i] = false;
		isspecificD[i] = false;
	}

	for (int i = 0; i < Nproper; i++) {
		dihedralindex1[i] = tempindex1[i];
		dihedralindex2[i] = tempindex2[i];
		dihedralindex3[i] = tempindex3[i];
		dihedralindex4[i] = tempindex4[i];
		if (verbose) printf("dihedral %d: %d %d %d %d\n", i+1, dihedralindex1[i], dihedralindex2[i], dihedralindex3[i], dihedralindex4[i]);
	}
	delete [] tempindex1;  delete [] tempindex2;  delete [] tempindex3;  delete [] tempindex4;

	for (int i = 0; i < Nimproper;) {
		p = fgets(s, 1024, f);
		if (p == NULL) THROWERROR2("Could not find all impropers in PSF file", psfname);
		count = TokenizeString(s, space, strings, 32);
		if ((count % 4) != 0) THROWERROR2("Uneven improper spacing in PSF file", psfname);
		for (int j = 0; j < count; j += 4) {
			dihedralindex1[Nproper+i] = StringToReal(strings[j]);
			dihedralindex2[Nproper+i] = StringToReal(strings[j+1]);
			dihedralindex3[Nproper+i] = StringToReal(strings[j+2]);
			dihedralindex4[Nproper+i] = StringToReal(strings[j+3]);
			if (verbose) printf("dihedral %d (improper), %d %d %d %d\n", Nproper+i+1, dihedralindex1[Nproper+i], dihedralindex2[Nproper+i], dihedralindex3[Nproper+i], dihedralindex4[Nproper+i]);
			i++;
		}
	}

	if (hascmap) {

		done = false;
		do {
			p = fgets(s, 1024, f);
			if (p == NULL) THROWERROR2("Could not find CMAP cross terms in PSF file", psfname);
			count = TokenizeString(s, space, strings, 32);
			if ((count >= 2) && (strncmp(strings[1], "!NCRTE", 6) == 0)) done = true;
		} while (!done);	

		Ncmap = StringToInteger(strings[0], 10);
		if (Ncmap == 0) {
			hascmap = false;
		} else {
			cindex = new cmapindex[Ncmap];
			cmapnumber = new integer[Ncmap];
			cmapassigned = new bool[Ncmap];
			for (int i = 0; i < Ncmap; i++) {
				cmapnumber[i] = 0;
				cmapassigned[i] = false;
			}



			for (int i = 0; i < Ncmap;) {
				p = fgets(s, 1024, f);
				if (p == NULL) THROWERROR2("Could not find all CMAP cross terms in PSF file", psfname);
				count = TokenizeString(s, space, strings, 32);
				if ((count % 8) != 0) THROWERROR2("Uneven CMAP term spacing in PSF file", psfname);
				for (int j = 0; j < count; j += 8) {
					cindex[i][0] = StringToReal(strings[j]);
					cindex[i][1] = StringToReal(strings[j+1]);
					cindex[i][2] = StringToReal(strings[j+2]);
					cindex[i][3] = StringToReal(strings[j+3]);
					cindex[i][4] = StringToReal(strings[j+4]);
					cindex[i][5] = StringToReal(strings[j+5]);	
					cindex[i][6] = StringToReal(strings[j+6]);
					cindex[i][7] = StringToReal(strings[j+7]);
					if (verbose) printf("CMAP term %d\n", i+1);
					i++;
				}
			}
		}
	}

	fclose(f);
// finished reading PSF file

	if (verbose) printf("Finished reading CHARMM PSF file\n");
	Natomsbase = Natoms;

	delete [] strings;

}


#define BONDLINE 0
#define ANGLELINE 1
#define DIHEDRALLINE 2
#define IMPROPERLINE 3
#define CMAPLINE 4
#define NBLINE 5
#define HBONDLINE 6
#define ENDLINE 7

void charmmmolecule::loadmoleculePRM(const char *parname, bool verbose) {

	integer Nkeywords = 8;
	const char *keywords[] = { "BOND", "ANGL", "DIHE", "IMPR", "CMAP", "NONB", "HBON", "END" };

	if (verbose) printf("Opening parameter file %s\n", parname);

	FILE* f = fopen(parname, "r");
	if (f == NULL) THROWERROR2("Could not open CHARMM parameter file", parname);
	char *p;
	char s[1024];

	int count;

	char *space = new char[8];
	strcpy(space, " ");
	char **strings = new char*[32];
	for (int i = 0; i < 32; i++) strings[i] = new char[32];

	// new
	int keyword = -1;
	do {
		p = fgets(s, 1024, f);
		if (p == NULL) return;
		count = TokenizeString(s, space, strings, 32);
	} while ((count == 0) || ((keyword = IsKeyword(strings[0], Nkeywords, keywords, 4)) < 0));
	// end new

	integer indices[32];

	char **atomtypenames = new char*[Natoms];
	for (int i = 0; i < Natoms; i++) {
		atomtypenames[i] = new char[8];
		strncpy(atomtypenames[i], atomtype[i], 6);	
	}

	bool done;

#define RETURN() { 	delete [] atomtypenames; delete [] strings; delete [] space; fclose(f); return;  }

	do {

	switch (keyword) {

		case BONDLINE : {

			done = false;
			do {
				p = fgets(s, 1024, f);
				if (p == NULL) RETURN();
				count = TokenizeString(s, space, strings, 32);
				if ((count >= 4) && (strings[0][0] != '!')) {
					for (int i = 0; i < Nbonds; i++) {
						indices[0] = bondindex1[i]-1;
						indices[1] = bondindex2[i]-1;
						if (SymmetricMatch(strings, atomtypenames, indices, 2)) {
							Kb[i] = StringToReal(strings[2]);
							b0[i] = StringToReal(strings[3]);
							bondassigned[i] = true;
							if (verbose) printf("bond %d: %d %d, %lf %lf\n", i+1, indices[0], indices[1], Kb[i], b0[i]);
						}
					}
				}
				if ((count > 0) && ((keyword = IsKeyword(strings[0], Nkeywords, keywords, 4)) >= 0)) done = true; 


			} while (!done);
			
			break;
		}
		
		case ANGLELINE : {

			done = false;
			do {
				p = fgets(s, 1024, f);
				if (p == NULL) RETURN();
				count = TokenizeString(s, space, strings, 32);
				if ((count >= 5) && (strings[0][0] != '!')) {
					for (int i = 0; i < Nangles; i++) {
						indices[0] = angleindex1[i]-1;
						indices[1] = angleindex2[i]-1;
						indices[2] = angleindex3[i]-1;

						if (SymmetricMatch(strings, atomtypenames, indices, 3)) {
							if (verbose) printf("angle %d: ", i+1);
							Ktheta[i] = StringToReal(strings[3]);
							theta0[i] = StringToReal(strings[4]);
							angleassigned[i] = true;
							if (verbose) printf("%lf %lf ", Ktheta[i], theta0[i]);
							if ((count > 5) && (strings[5][0] != '!')) {
								Kub[i] = StringToReal(strings[5]);
								S0[i] = StringToReal(strings[6]);
								if (verbose) printf("%lf %lf ", Kub[i], S0[i]);
							}
							if (verbose) printf("\n");
						}
					}
				} 
				if ((count > 0) && ((keyword = IsKeyword(strings[0], Nkeywords, keywords, 4)) >= 0)) done = true; 
			} while (!done);
			break;
		}
		
		case DIHEDRALLINE : {

			done = false;
			do {
				p = fgets(s, 1024, f);
				if (p == NULL) RETURN();
				count = TokenizeString(s, space, strings, 32);
				if ((count >= 7) && (strings[0][0] != '!')) {
					for (int i = 0; i < Nproper; i++) {
						indices[0] = dihedralindex1[i]-1;
						indices[1] = dihedralindex2[i]-1;
						indices[2] = dihedralindex3[i]-1;
						indices[3] = dihedralindex4[i]-1;
						if (SpecificSymmetricMatch(strings, atomtypenames, indices, 4)) {
							if (verbose) printf("dihedral %d\n", i+1);
							if ((isspecificD[i] == false) && (dihedralassigned[i] == true)) {
								dihedralV[i].N = 0;
							}
							adddihedral(i, StringToReal(strings[4]), StringToInteger(strings[5], 10), StringToReal(strings[6]), 0.0, 0.0);
							dihedralassigned[i] = true;
							isspecificD[i] = true;
						}
						else if ((SymmetricMatch(strings, atomtypenames, indices, 4)) && (!isspecificD[i])) {
							if (verbose) printf("dihedral %d\n", i+1);
							dihedralassigned[i] = true;
							adddihedral(i, StringToReal(strings[4]), StringToInteger(strings[5], 10), StringToReal(strings[6]), 0.0, 0.0);
						}
					}
				} 
				if ((count > 0) && ((keyword = IsKeyword(strings[0], Nkeywords, keywords, 4)) >= 0)) done = true; 		
			} while (!done);
			break;
			
		}
		
		case IMPROPERLINE : {
			done = false;
			do {
				p = fgets(s, 1024, f);
				if (p == NULL) RETURN();
				count = TokenizeString(s, space, strings, 32);
				if ((count >= 7) && (strings[0][0] != '!')) {
					for (int i = Nproper; i < Nproper+Nimproper; i++) {
						indices[0] = dihedralindex1[i]-1;
						indices[1] = dihedralindex2[i]-1;
						indices[2] = dihedralindex3[i]-1;
						indices[3] = dihedralindex4[i]-1;
						if (SpecificSymmetricMatch(strings, atomtypenames, indices, 4)) {
							if (verbose) printf("improper dihedral %d\n", i+1);
							if ((isspecificD[i] == false) && (dihedralassigned[i] == true)) {
								dihedralV[i].N = 0;
							}
							adddihedral(i, 0.0, 0, 0.0, StringToReal(strings[6]), StringToReal(strings[4]));
							dihedralassigned[i] = true;
							isspecificD[i] = true;
						}
						else if ((SymmetricMatch(strings, atomtypenames, indices, 4)) && (!isspecificD[i])) {
							if (verbose) printf("improper dihedral %d\n", i+1);
							dihedralassigned[i] = true;
							adddihedral(i, 0.0, 0, 0.0, StringToReal(strings[6]), StringToReal(strings[4]));
						}
					}
				} 
				if ((count > 0) && ((keyword = IsKeyword(strings[0], Nkeywords, keywords, 4)) >= 0)) done = true;
			} while (!done);			
			break;
		}

		case CMAPLINE : {
			integer mlength;
			cmap = new real**[CMAPMAXN];
			cmaplength = new integer[CMAPMAXN];
			ncmaps = 0;
			char **oldstrings = new char*[8];
			for (int i = 0; i < 8; i++) oldstrings[i] = new char[32];
			if (verbose) printf("Looking for CMAP parameters\n");
			done = false;
			do {
				p = fgets(s, 1024, f);
				if (p == NULL) RETURN();
				count = TokenizeString(s, space, strings, 32);
				if ((count >= 9) && (strings[0][0] != '!') && (IsKeyword(strings[0], Nkeywords, keywords, 4) < 0)) {
					for (int i = 0; i < 8; i++) strcpy(oldstrings[i], strings[i]);
					mlength = StringToInteger(strings[8], 10);
					cmap[ncmaps] = new real*[mlength];
					cmaplength[ncmaps] = mlength;
					for (integer X = 0; X < mlength; X++) {
						cmap[ncmaps][X] = new real[mlength];
						for (integer Y = 0; Y < mlength; Y++) cmap[ncmaps][X][Y] = 0.0;
					}
					integer total = 0;
					do {
						p = fgets(s, 1024, f);
						if (p == NULL) THROWERROR2("End of line while looking for CMAP parameters in file", parname);
						count = TokenizeString(s, space, strings, 32);
						real value;
						if (count > 3) for (integer j = 0; (j < count) && (strings[j][0] != '!') && (strlen(strings[j]) != 0) && (total < mlength*mlength); j++) {
							if (strlen(strings[j]) > 3) {
								value = StringToReal(strings[j]);
								cmap[ncmaps][total/mlength][total%mlength] = value;
								total++;
							}
						}
					} while (total < mlength*mlength);					
					for (int i = 0; i < Ncmap; i++) {
						indices[0] = cindex[i][0]-1;
						indices[1] = cindex[i][1]-1;
						indices[2] = cindex[i][2]-1;
						indices[3] = cindex[i][3]-1;
						indices[4] = cindex[i][4]-1;
						indices[5] = cindex[i][5]-1;
						indices[6] = cindex[i][6]-1;
						indices[7] = cindex[i][7]-1;
						if (SpecificMatch(oldstrings, atomtypenames, indices, 8)) {
							if (verbose) printf("CMAP term %d\n", i+1);
							cmapassigned[i] = true;
							cmapnumber[i] = ncmaps;
						}
					}
					ncmaps++;
				}
				if ((count > 0) && ((keyword = IsKeyword(strings[0], Nkeywords, keywords, 4)) >= 0)) done = true;
			} while (!done);
			delete [] oldstrings;
			for (int i = 0; i < ncmaps; i++) {
				if (verbose) printf("CMAP %d\n", i);
				for (int j = 0; j < cmaplength[i]; j++) {
					for (int k = 0; k < cmaplength[i]; k++) if (verbose) printf("%lf ", cmap[i][j][k]);
					if (verbose) printf("\n");
				}
				if (verbose) printf("\n");
			}
			break;
		}

		case NBLINE : {
			done = false;
			do {
				p = fgets(s, 1024, f);
				if (p == NULL) RETURN();
				count = TokenizeString(s, space, strings, 32);
				if ((count >= 4) && (strings[0][0] != '!')) {
					for (int i = 0; i < Natoms; i++) {
						indices[0] = i;
						if (SymmetricMatch(strings, atomtypenames, indices, 1)) {
							if (verbose) printf("VDW parameters %d\n", i+1);
							vdwassigned[i] = true;
							sigma[i] = 2.0*StringToReal(strings[3])/TWOROOTONESIXTH;
							epsilon[i] = -StringToReal(strings[2]);
							if (verbose) printf("VDW atom %d: %lf %lf\n", i+1, sigma[i], epsilon[i]);
						}
					}
				} 
				if ((count > 0) && ((keyword = IsKeyword(strings[0], Nkeywords, keywords, 4)) >= 0)) done = true;
			} while (!done);
			break;
		}
		
		case ENDLINE : {
			RETURN(); break;
		}
		
		default : {
			if (verbose) printf("Valid but unimplemented line in PSF file, %s\n", strings[0]);
			done = false;
			do {
				p = fgets(s, 1024, f);
				if (p == NULL) RETURN();
				if ((count > 0) && ((keyword = IsKeyword(strings[0], Nkeywords, keywords, 4)) >= 0)) done = true;		
			} while (!done);
			break;
		}
	}
	} while (p != NULL);

	if (verbose) printf("Finished reading CHARMM parameter data\n");
	fclose(f);
	RETURN();
	
}


vector <int> charmmmolecule::atomtypeatoms(const char *s) {
	vector <int> atlist;
	for (int i = 0; i < Natoms; i++) {
		if (strcasecmp(atomtype[i], s) == 0) atlist.push_back(i);
	}
	return atlist;	
}


vector <int> charmmmolecule::segmentresidueatoms(const char *s, integer res) {
	vector <int> segreslist;
	for (int i = 0; i < Natoms; i++) {
		if ((strcasecmp(segname[i], s) == 0) && (resid[i] == res)) segreslist.push_back(i);
	}
	return segreslist;
}


vector <int> charmmmolecule::segmentatoms(const char *s) {
	vector <int> seglist;
	for (int i = 0; i < Natoms; i++) {
		if (strcasecmp(segname[i], s) == 0) seglist.push_back(i);
	}
	return seglist;
}


vector <int> charmmmolecule::getNsegments() {
	vector <int> segments;
	if (Natoms == 0) return segments;
	segments.push_back(0);
	for (int i = 1; i < Natoms; i++) {
		if (strcasecmp(segname[i]-1, segname[i]) != 0) segments.push_back(i);
	}
	return segments;
}

