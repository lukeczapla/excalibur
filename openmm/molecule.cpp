#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "molecule.h"

using namespace OpenMM;
using namespace molecules;
using namespace std;


molecule::molecule(const realmolecule& m) : systemBuilt(false), Ngroups(1), usegbsa(false), usegbvi(false), usecm(false), realmolecule(m) {
	atomgroups = new int[Natoms];
	for (int i = 0; i < Natoms; i++) atomgroups[i] = 0;
}


void molecule::randomizePosition() {

	double xt = 0.0, yt = 0.0, zt = 0.0;
	for (int i = 0; i < Natoms; i++) {
		xt += x[i];
		yt += y[i];
		zt += z[i];
	}
	xt /= (double)Natoms;
	yt /= (double)Natoms;
	zt /= (double)Natoms;
	
	double xr = boxL[0] * drand48();
	double yr = boxL[1] * drand48();
	double zr = boxL[2] * drand48();

	double gam = M_PI*drand48();
	double phi = M_PI*drand48();
	double omega = 2.0*M_PI*drand48();

	double sp, cp, sm, cm, sg, cg;

	sp = sin(omega/2.0+phi); cp = cos(omega/2.0+phi); sm = sin(omega/2.0-phi);
	cm = cos(omega/2.0-phi); sg = sin(gam); cg = cos(gam);

	double xn, yn, zn;

	for (int i = 0; i < Natoms; i++) {
		xn = (x[i] - xt) * (cm*cg*cp-sm*sp) + (y[i] - yt) * (sm*cg*cp+cm*sp) + (z[i] - zt) * -sg*cp + xr;
		yn = (x[i] - xt) * (-cm*cg*sp-sm*cp) + (y[i] - yt) * (-sm*cg*sp+cm*cp) + (z[i] - zt) * sg*sp + yr;
		zn = (x[i] - xt) * cm*sg + (y[i] - yt) * sm*sg + (z[i] - zt) * cg + zr;
		x[i] = xn; y[i] = yn; z[i] = zn;
	}

}


void molecule::addConstraintGroup(int index, double scaling) {
	constraintGroups.push_back(index);
	constraintScaling.push_back(scaling);
}


double molecule::getPotentialEnergy() {
	if (!hasContext) {
		printf("Error, no context\n"); exit(0);
	}
	State state = context->getState(State::Energy);
	return state.getPotentialEnergy();
}



void molecule::removeConstraints(int index) {
	int element = -1;
	for (int i = 0; (i < constraintGroups.size()) && (element == -1); i++) {
		if (constraintGroups[i] == index) element = i;
	}
	if (element == -1) {
		cout << "Constraint to remove not found" << endl;
		return;
	}
	if (systemBuilt) {
		atomconstraints[element]->setGlobalParameterDefaultValue(0, 0.0);
		if (hasContext) atomconstraints[element]->updateParametersInContext(*context);
	} else {
		constraintGroups.erase(constraintGroups.begin() + element);
		constraintScaling.erase(constraintScaling.begin() + element);
	}
}


int molecule::addGroup(vector <int> a) {
	for (int i = 0; i < Natoms; i++) {
		for (int j = 0; j < a.size(); j++) if (a[j] == i) atomgroups[i] = Ngroups;
	}
	return Ngroups++;
}


void molecule::removeGroup(int index) {
	if (index >= Ngroups) {
		cout << "Out of range group attempted to be removed" << endl;
		return;
	}
	for (int i = 0; i < Natoms; i++) {
		if (atomgroups[i] == index) atomgroups[i] = 0;
		if (atomgroups[i] > index) --atomgroups[i];
	}
}


void molecule::setDihedralGroup(int group, int n) {
	if (group >= Ngroups) {
		cout << "Illegal assignment of non-existant group " << group << endl;
		return;
	}
	torsion[group]->setForceGroup(n);
}


System* molecule::buildSystem(bool verbose) {

	system = new System();
	nonbonded = new NonbondedForce();
	if (usegbsa && usegbvi) {
		THROWERROR("GBSA OBC2 and GBVI are mutually exclusive");
	}
	if (usegbsa) gbsa = new GBSAOBCForce();
	if (usegbvi) obc = new OBC2Force();
	bonded = new HarmonicBondForce();
	angle = new HarmonicAngleForce();
	torsion = new PeriodicTorsionForce*[Ngroups];
	impropertorsion = new CustomTorsionForce("Kpsi*(theta-psi0)^2");
	if (hascmap) cmaptorsion = new CMAPTorsionForce();
	for (int i = 0; i < Ngroups; i++) torsion[i] = new PeriodicTorsionForce();
	
	cout << "Total mass " << totalmass << " and " << Ngroups << " groups are defined" << endl;

	if (verbose) printf("Adding atoms\n");
	for (int n = 0; n < Natoms; n++) {
		if (verbose) printf("Atom %d, mass %lf charge %lf sigma %lf epsilon %lf\n", n, mass[n], charge[n], sigma[n], epsilon[n]);
		system->addParticle(mass[n]);
		if (epsilon[n] <= 0.0) THROWERROR("Invalid epsilon value <= 0");
		if (usegbsa) {
			gbsa->addParticle(charge[n], radii[n] * NmPerAngstrom, 1.0);
		}
		else if (usegbvi) {
			vector <double> parameters(5);
			parameters[0] = charge[n];
			parameters[1] = radii[n] * NmPerAngstrom;
			parameters[2] = 1.0;
			parameters[3] = gamma[n] * KJPerKcal;
			parameters[4] = gamma2[n] * KJPerKcal;
			obc->addParticle(parameters);
		}
		nonbonded->addParticle(charge[n], sigma[n] * NmPerAngstrom, epsilon[n] * KJPerKcal);
		const Vec3 pos(x[n] * NmPerAngstrom, y[n] * NmPerAngstrom, z[n] * NmPerAngstrom);
		rvec.push_back(pos);
	}

	// add non-hydrogen atom constraints
	if (constraintGroups.size() > 0) {
		atomconstraints = new CustomExternalForce*[constraintGroups.size()];
		vector<double> constraintparameters(3);
		for (int n = 0; n < constraintGroups.size(); n++) {
			atomconstraints[n] = new CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)");
			atomconstraints[n]->addGlobalParameter("k", constraintScaling[n]);
			atomconstraints[n]->addPerParticleParameter("x0");
			atomconstraints[n]->addPerParticleParameter("y0");
			atomconstraints[n]->addPerParticleParameter("z0");
			for (int i = 0; i < Natoms; i++) {
				if ((atomgroups[i] == n) && (mass[i] > 2.0)) {
					constraintparameters[0] = x[i] * NmPerAngstrom;
					constraintparameters[1] = y[i] * NmPerAngstrom;
					constraintparameters[2] = z[i] * NmPerAngstrom;
					atomconstraints[n]->addParticle(i, constraintparameters);
				}
			}
		}
	}

	vector < pair<int,int> > bondPairs;
	for (int n = 0; n < Nbonds; n++) {
		if (verbose) printf("Bond %d, b0 %lf Kb %lf\n", n, b0[n], Kb[n]);
		if ((mass[bondindex1[n]-1] < 2.0) || (mass[bondindex2[n]-1] < 2.0)) system->addConstraint(bondindex1[n]-1, bondindex2[n]-1, b0[n] * NmPerAngstrom);
		if (Kb[n] > 0.0) {
			bonded->addBond(bondindex1[n]-1, bondindex2[n]-1, b0[n] * NmPerAngstrom, 2.0*Kb[n] * KJPerKcal * AngstromsPerNm * AngstromsPerNm);
//			if (usegbvi) gbvi->addBond(bondindex1[n]-1, bondindex2[n]-1, b0[n] * NmPerAngstrom);
		}
		bondPairs.push_back(make_pair(bondindex1[n]-1, bondindex2[n]-1));
	}

	nonbonded->createExceptionsFromBonds(bondPairs, scaling14es, scaling14lj);
	
    for (int n = 0; n < Nangles; n++) {
		if (verbose) printf("Angle %d, theta0 %lf Ktheta %lf Kub %lf S0 %lf\n", n, theta0[n], Ktheta[n], Kub[n], S0[n]);
		if (Ktheta[n] > 0.0) angle->addAngle(angleindex1[n]-1, angleindex2[n]-1, angleindex3[n]-1, theta0[n] * RadiansPerDegree, 2.0*Ktheta[n] * KJPerKcal);
		if ((Kub[n] > 0.0) && (S0[n] > 0.0)) {
			if ((mass[angleindex1[n]-1] < 2.0) && (mass[angleindex3[n]-1] < 2.0)) system->addConstraint(angleindex1[n]-1, angleindex3[n]-1, S0[n] * NmPerAngstrom);
			bonded->addBond(angleindex1[n]-1, angleindex3[n]-1, S0[n] * NmPerAngstrom, 2.0*Kub[n] * KJPerKcal * AngstromsPerNm * AngstromsPerNm);
		}
	}

	
	if (verbose) printf("Adding dihedrals %d\n", NuniqueD);
	impropertorsion->addPerTorsionParameter("Kpsi");
	impropertorsion->addPerTorsionParameter("psi0");
	vector<double> improperparameters(2);
	for (int j = 0; j < Ngroups; j++) for (int n = 0; n < NuniqueD; n++) {
		if (atomgroups[dihedralindex1[n]-1] == j) for (int i = 0; i < dihedralV[n].N; i++) {
			if (dihedralV[n].Kchi[i] != 0.0) {
				if (verbose) printf("Proper Dihedral %d, n %d delta %lf Kchi %lf\n", n, dihedralV[n].nD[i], dihedralV[n].delta[i], dihedralV[n].Kchi[i]);
				torsion[j]->addTorsion(dihedralindex1[n]-1, dihedralindex2[n]-1, dihedralindex3[n]-1, dihedralindex4[n]-1,
					(int)dihedralV[n].nD[i], (dihedralV[n].delta[i]) * RadiansPerDegree, dihedralV[n].Kchi[i] * KJPerKcal);
			}
			if (dihedralV[n].Kpsi[i] != 0.0) {
				if (verbose) printf("Improper Dihedral %d, Kpsi %lf psi0 %lf\n", n, dihedralV[n].Kpsi[i], dihedralV[n].psi0[i]);
				improperparameters[0] = dihedralV[n].Kpsi[i] * KJPerKcal;
				improperparameters[1] = dihedralV[n].psi0[i] * RadiansPerDegree;
				impropertorsion->addTorsion(dihedralindex1[n]-1, dihedralindex2[n]-1, dihedralindex3[n]-1, dihedralindex4[n]-1, improperparameters);		
			}
		}
	}

	if (hascmap) {
		if (verbose) printf("Adding CMAP cross terms\n");
		for (int n = 0; n < ncmaps; n++) {
			vector<double> emap(cmaplength[n] * cmaplength[n]);
			for (int X = cmaplength[n]/2; X < cmaplength[n]; X++) {
				for (int Y = cmaplength[n]/2; Y < cmaplength[n]; Y++) {
					emap[(X-cmaplength[n]/2)+(Y-cmaplength[n]/2)*cmaplength[n]] = KJPerKcal*cmap[n][X][Y];
				}
				for (int Y = 0; Y < cmaplength[n]/2; Y++) {
					emap[(X-cmaplength[n]/2)+(Y+cmaplength[n]/2)*cmaplength[n]] = KJPerKcal*cmap[n][X][Y];
				}
			}
			for (int X = 0; X < cmaplength[n]/2; X++) {
				for (int Y = cmaplength[n]/2; Y < cmaplength[n]; Y++) {
					emap[(X+cmaplength[n]/2)+(Y-cmaplength[n]/2)*cmaplength[n]] = KJPerKcal*cmap[n][X][Y];
				}
				for (int Y = 0; Y < cmaplength[n]/2; Y++) {
					emap[(X+cmaplength[n]/2)+(Y+cmaplength[n]/2)*cmaplength[n]] = KJPerKcal*cmap[n][X][Y];
				}
			}
	
			cmaptorsion->addMap(cmaplength[n], emap);
		}
	}
	for (int n = 0; n < Ncmap; n++) {
		if (cmapassigned[n]) cmaptorsion->addTorsion(cmapnumber[n], cindex[n][0]-1, cindex[n][1]-1, cindex[n][2]-1, cindex[n][3]-1, 
			cindex[n][4]-1, cindex[n][5]-1, cindex[n][6]-1, cindex[n][7]-1);
		if (verbose) printf("Assigned CMAP term %d to crossterm %d\n", cmapnumber[n]+1, n+1);
	}

	nonbonded->setForceGroup(0);
	bonded->setForceGroup(1);
	angle->setForceGroup(1);
	system->addForce(nonbonded);
	system->addForce(bonded);
	system->addForce(angle);
	for (int i = 0; i < Ngroups; i++) system->addForce(torsion[i]);
	system->addForce(impropertorsion);
	if (hascmap) {
		system->addForce(cmaptorsion);
	}
	
	if (usegbsa) {
		gbsa->setForceGroup(4);
		system->addForce(gbsa);
	} else if (usegbvi) {
		obc->setForceGroup(4);
		system->addForce(obc);
	}
	
	if (usecm) {
		cmmotion = new CMMotionRemover(1);
		system->addForce(cmmotion);
	}

	if (verbose) printf("Finishing up\n");
	systemBuilt = true;
	return system;

}



System* molecule::buildSystemClone(int N, bool verbose) {

	expandcoordinates(N*Natoms);

	system = new System();
	nonbonded = new NonbondedForce();
	if (usegbsa && usegbvi) {
		THROWERROR("GBSA OBC2 and GBVI are mutually exclusive");
	}
	if (usegbsa) gbsa = new GBSAOBCForce();
	if (usegbvi) obc = new OBC2Force();
	bonded = new HarmonicBondForce();
	angle = new HarmonicAngleForce();
	torsion = new PeriodicTorsionForce*[Ngroups];
	impropertorsion = new CustomTorsionForce("Kpsi*(theta-psi0)^2");
	if (hascmap) cmaptorsion = new CMAPTorsionForce();
	for (int i = 0; i < Ngroups; i++) torsion[i] = new PeriodicTorsionForce();
	
	cout << "Total mass " << totalmass << " and " << Ngroups << " groups are defined" << endl;

	for (int X = 0; X < N; X++) {
	if (verbose) printf("Adding atoms\n");
	randomizePosition();
	for (int n = 0; n < Natoms; n++) {
		if (verbose) printf("Atom %d, mass %lf charge %lf sigma %lf epsilon %lf\n", n, mass[n], charge[n], sigma[n], epsilon[n]);
		system->addParticle(mass[n]);
		if (epsilon[n] <= 0.0) THROWERROR("Invalid epsilon value <= 0");
		if (usegbsa) {
			gbsa->addParticle(charge[n], radii[n] * NmPerAngstrom, 1.0);
		}
		else if (usegbvi) {
			vector <double> parameters(5);
			parameters[0] = charge[n];
			parameters[1] = radii[n] * NmPerAngstrom;
			parameters[2] = 1.0;
			parameters[3] = gamma[n] * KJPerKcal;
			parameters[4] = gamma2[n] * KJPerKcal;
			// debugging statement below
			printf("GBVI Particle %d : charge %lf, radius %lf, gamma %lf, gamma2 %lf, sigma %lf, epsilon %lf\n", X*Natoms+n, charge[n], radii[n], gamma[n], gamma2[n], sigma[n], epsilon[n]);
			obc->addParticle(parameters);
		}
		nonbonded->addParticle(charge[n], sigma[n] * NmPerAngstrom, epsilon[n] * KJPerKcal);
		Vec3 pos(x[n] * NmPerAngstrom, y[n] * NmPerAngstrom, z[n] * NmPerAngstrom);
		rvec.push_back(pos);
	}

	// add non-hydrogen atom constraints
	if (constraintGroups.size() > 0) {
		atomconstraints = new CustomExternalForce*[constraintGroups.size()];
		vector<double> constraintparameters(3);
		for (int n = 0; n < constraintGroups.size(); n++) {
			atomconstraints[n] = new CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)");
			atomconstraints[n]->addGlobalParameter("k", constraintScaling[n]);
			atomconstraints[n]->addPerParticleParameter("x0");
			atomconstraints[n]->addPerParticleParameter("y0");
			atomconstraints[n]->addPerParticleParameter("z0");
			for (int i = 0; i < Natoms; i++) {
				if ((atomgroups[i] == n) && (mass[i] > 2.0)) {
					constraintparameters[0] = x[i] * NmPerAngstrom;
					constraintparameters[1] = y[i] * NmPerAngstrom;
					constraintparameters[2] = z[i] * NmPerAngstrom;
					atomconstraints[n]->addParticle(X*Natoms+i, constraintparameters);
				}
			}
		}
	}

	vector < pair<int,int> > bondPairs;
	for (int n = 0; n < Nbonds; n++) {
		if (verbose) printf("Bond %d, b0 %lf Kb %lf\n", n, b0[n], Kb[n]);
		if ((mass[bondindex1[n]-1] < 2.0) || (mass[bondindex2[n]-1] < 2.0)) system->addConstraint(X*Natoms+bondindex1[n]-1, X*Natoms+bondindex2[n]-1, b0[n] * NmPerAngstrom);
		if (Kb[n] > 0.0) {
			bonded->addBond(X*Natoms+bondindex1[n]-1, X*Natoms+bondindex2[n]-1, b0[n] * NmPerAngstrom, 2.0*Kb[n] * KJPerKcal * AngstromsPerNm * AngstromsPerNm);
//			if (usegbvi) gbvi->addBond(bondindex1[n]-1, bondindex2[n]-1, b0[n] * NmPerAngstrom);
		}
		bondPairs.push_back(make_pair(X*Natoms+bondindex1[n]-1, X*Natoms+bondindex2[n]-1));
	}

	nonbonded->createExceptionsFromBonds(bondPairs, scaling14es, scaling14lj);
	
    for (int n = 0; n < Nangles; n++) {
		if (verbose) printf("Angle %d, theta0 %lf Ktheta %lf Kub %lf S0 %lf\n", n, theta0[n], Ktheta[n], Kub[n], S0[n]);
		if (Ktheta[n] > 0.0) angle->addAngle(X*Natoms+angleindex1[n]-1, X*Natoms+angleindex2[n]-1, X*Natoms+angleindex3[n]-1, theta0[n] * RadiansPerDegree, 2.0*Ktheta[n] * KJPerKcal);
		if ((Kub[n] > 0.0) && (S0[n] > 0.0)) {
			if ((mass[angleindex1[n]-1] < 2.0) && (mass[angleindex3[n]-1] < 2.0)) system->addConstraint(X*Natoms+angleindex1[n]-1, X*Natoms+angleindex3[n]-1, S0[n] * NmPerAngstrom);
			bonded->addBond(X*Natoms+angleindex1[n]-1, X*Natoms+angleindex3[n]-1, S0[n] * NmPerAngstrom, 2.0*Kub[n] * KJPerKcal * AngstromsPerNm * AngstromsPerNm);
		}
	}

	
	if (verbose) printf("Adding dihedrals %d\n", NuniqueD);
	impropertorsion->addPerTorsionParameter("Kpsi");
	impropertorsion->addPerTorsionParameter("psi0");
	vector<double> improperparameters(2);
	for (int j = 0; j < Ngroups; j++) for (int n = 0; n < NuniqueD; n++) {
		if (atomgroups[dihedralindex1[n]-1] == j) for (int i = 0; i < dihedralV[n].N; i++) {
			if (dihedralV[n].Kchi[i] != 0.0) {
				if (verbose) printf("Proper Dihedral %d, n %d delta %lf Kchi %lf\n", n, dihedralV[n].nD[i], dihedralV[n].delta[i], dihedralV[n].Kchi[i]);
				torsion[j]->addTorsion(X*Natoms+dihedralindex1[n]-1, X*Natoms+dihedralindex2[n]-1, X*Natoms+dihedralindex3[n]-1, X*Natoms+dihedralindex4[n]-1,
					(int)dihedralV[n].nD[i], (dihedralV[n].delta[i]) * RadiansPerDegree, dihedralV[n].Kchi[i] * KJPerKcal);
			}
			if (dihedralV[n].Kpsi[i] != 0.0) {
				if (verbose) printf("Improper Dihedral %d, Kpsi %lf psi0 %lf\n", n, dihedralV[n].Kpsi[i], dihedralV[n].psi0[i]);
				improperparameters[0] = dihedralV[n].Kpsi[i] * KJPerKcal;
				improperparameters[1] = dihedralV[n].psi0[i] * RadiansPerDegree;
				impropertorsion->addTorsion(X*Natoms+dihedralindex1[n]-1, X*Natoms+dihedralindex2[n]-1, X*Natoms+dihedralindex3[n]-1, X*Natoms+dihedralindex4[n]-1, improperparameters);		
			}
		}
	}

	if (hascmap) {
		if (verbose) printf("Adding CMAP cross terms\n");
		for (int n = 0; n < ncmaps; n++) {
			vector<double> emap(cmaplength[n] * cmaplength[n]);
			for (int X = cmaplength[n]/2; X < cmaplength[n]; X++) {
				for (int Y = cmaplength[n]/2; Y < cmaplength[n]; Y++) {
					emap[(X-cmaplength[n]/2)+(Y-cmaplength[n]/2)*cmaplength[n]] = KJPerKcal*cmap[n][X][Y];
				}
				for (int Y = 0; Y < cmaplength[n]/2; Y++) {
					emap[(X-cmaplength[n]/2)+(Y+cmaplength[n]/2)*cmaplength[n]] = KJPerKcal*cmap[n][X][Y];
				}
			}
			for (int X = 0; X < cmaplength[n]/2; X++) {
				for (int Y = cmaplength[n]/2; Y < cmaplength[n]; Y++) {
					emap[(X+cmaplength[n]/2)+(Y-cmaplength[n]/2)*cmaplength[n]] = KJPerKcal*cmap[n][X][Y];
				}
				for (int Y = 0; Y < cmaplength[n]/2; Y++) {
					emap[(X+cmaplength[n]/2)+(Y+cmaplength[n]/2)*cmaplength[n]] = KJPerKcal*cmap[n][X][Y];
				}
			}
	
			cmaptorsion->addMap(cmaplength[n], emap);
		}
	}
	for (int n = 0; n < Ncmap; n++) {
		if (cmapassigned[n]) cmaptorsion->addTorsion(cmapnumber[n], X*Natoms+cindex[n][0]-1, X*Natoms+cindex[n][1]-1, X*Natoms+cindex[n][2]-1, X*Natoms+cindex[n][3]-1, 
			X*Natoms+cindex[n][4]-1, X*Natoms+cindex[n][5]-1, X*Natoms+cindex[n][6]-1, X*Natoms+cindex[n][7]-1);
		if (verbose) printf("Assigned CMAP term %d to crossterm %d\n", cmapnumber[n]+1, n+1);
	}

	}

	nonbonded->setForceGroup(0);
	bonded->setForceGroup(1);
	angle->setForceGroup(1);
	system->addForce(nonbonded);
	system->addForce(bonded);
	system->addForce(angle);
	for (int i = 0; i < Ngroups; i++) system->addForce(torsion[i]);
	system->addForce(impropertorsion);
	if (hascmap) {
		system->addForce(cmaptorsion);
	}
	
	if (usegbsa) {
		gbsa->setForceGroup(4);
		system->addForce(gbsa);
	} else if (usegbvi) {
		obc->setForceGroup(0);
		system->addForce(obc);
	}
	
	if (usecm) {
		cmmotion = new CMMotionRemover(1);
		system->addForce(cmmotion);
	}

	if (verbose) printf("Finishing up\n");
	systemBuilt = true;
	Natoms *= N;
	
	return system;

}

