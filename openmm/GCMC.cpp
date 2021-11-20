#include "GCMC.h"
#include <iostream>
#include <math.h>

using namespace OpenMM;


gcmc::gcmc(const molecule& m) : simulation(m), species(NULL), Nadded(0) {
	baseAtoms = Natoms;
}

gcmc::~gcmc() {
	
}


void gcmc::setup() {
	
//	CMMotionRemover &cmmotion = *new CMMotionRemover(1);
	cout << "Adding thermostat and barostat" << endl;
	thermostat = new AndersenThermostat(temperature, temperatureFrequency);
	system->addForce(thermostat);
	if (usePBCbox) {

		nonbonded->setNonbondedMethod(NonbondedForce::PME);
		nonbonded->setCutoffDistance(cutoff);
		nonbonded->setEwaldErrorTolerance(1e-3);
		nonbonded->setUseDispersionCorrection(true);
	
		system->setDefaultPeriodicBoxVectors(Vec3(Lx,0,0), Vec3(0,Ly,0), Vec3(0,0,Lz));
		cout << "Box vectors (nm) " << Lx << ' ' << Ly << ' ' << Lz << endl;
		if (usePressure) {
			barostat = new MonteCarloBarostat(pressure, temperature, pressureFrequency);
			system->addForce(barostat);
		}
/*		gbsa = new GBSAOBCForce();
		gbsa->setSolventDielectric(80);
		gbsa->setSoluteDielectric(1);
		for (int i = 0; i < Natoms; i++) {
			gbsa->addParticle(charge[i], (sigma[i]/2.0)*NmPerAngstrom, 1.0);
		}
		gbsa->setNonbondedMethod(GBSAOBCForce::CutoffPeriodic);
		gbsa->setCutoffDistance(1.2);
		system->addForce(gbsa);*/
	} else {
		nonbonded->setNonbondedMethod(NonbondedForce::NoCutoff);
	}
	
	integrator = new VerletIntegrator(timestep_fs * PsPerFs);

	cout << "Setting up context" << endl;
    if (platformSpecified) context = new Context(*system, *integrator, Platform::getPlatformByName(platformName));
	else context = new Context(*system, *integrator);
	hasContext = true;

	cout << "Running on platform " << context->getPlatform().getName() << endl;

	cout << "Setting up positions" << endl;
    context->setPositions(rvec);

	State state = context->getState(State::Energy+State::Positions);
	cout << "Potential Energy = " << state.getPotentialEnergy() * KcalPerKJ <<  " kcal/mol" << endl;
	if (visualize) drawmolecule("localhost", 7001);
	
}


void gcmc::setSpecies(simulation *m) { 
	species = m;
	oldx = new double[species->Natoms];
	oldy = new double[species->Natoms];
	oldz = new double[species->Natoms];
}


void gcmc::addMolecule(bool restore) {

	cout << "Inserting" << endl;

	double randX = Lx * (drand48()-0.5);
	double randY = Ly * (drand48()-0.5);
	double randZ = Lz * (drand48()-0.5);

	double avgX = 0.0, avgY = 0.0, avgZ = 0.0;
	for (int n = 0; n < species->Natoms; n++) {
		avgX += species->x[n];
		avgY += species->y[n];
		avgZ += species->z[n];
	}
	avgX /= species->Natoms; avgY /= species->Natoms; avgZ /= species->Natoms;

	list<int> particleList;
	list<int> nonbondedList;
	for (int n = 0; n < species->Natoms; n++) {
		particleList.push_back(system->addParticle(species->mass[n]));
		if (epsilon[n] <= 0.0) THROWERROR("Invalid epsilon value <= 0");
		nonbondedList.push_back(nonbonded->addParticle(species->charge[n], species->sigma[n] * NmPerAngstrom, species->epsilon[n] * KJPerKcal));
		if (restore) {
			Vec3 pos(oldx[n], oldy[n], oldz[n]);
			rvec.push_back(pos);
		}
		else {
			Vec3 pos((species->x[n] - avgX) * NmPerAngstrom + randX, (species->y[n] - avgY) * NmPerAngstrom + randY, (species->z[n] - avgZ) * NmPerAngstrom + randZ);
			rvec.push_back(pos);
		}
	}

	
	vector < pair<int,int> > bondPairs;
	list <int> bondList;
	list <int> constraintList;
	for (int n = 0; n < species->Nbonds; n++) {
		if ((species->mass[species->bondindex1[n]-1] < 2.0) || (species->mass[species->bondindex2[n]-1] < 2.0)) constraintList.push_back(system->addConstraint(Natoms+species->bondindex1[n]-1, Natoms+species->bondindex2[n]-1, species->b0[n] * NmPerAngstrom));
		if (species->Kb[n] > 0.0) bondList.push_back(bonded->addBond(Natoms+species->bondindex1[n]-1, Natoms+species->bondindex2[n]-1, 
			species->b0[n] * NmPerAngstrom, 2.0*species->Kb[n] * KJPerKcal * AngstromsPerNm * AngstromsPerNm));
		bondPairs.push_back(make_pair(Natoms+species->bondindex1[n]-1, Natoms+species->bondindex2[n]-1));
	}
	
	int Nprior = nonbonded->getNumExceptions();
	nonbonded->createExceptionsFromBonds(bondPairs, scaling14es, scaling14lj);
	int Nnew = nonbonded->getNumExceptions();
	
	list <int> exceptionList;
	for (int i = Nprior; i < Nnew; i++) exceptionList.push_back(i);
	
	list <int> angleList;
    for (int n = 0; n < species->Nangles; n++) {
		if (species->Ktheta[n] > 0.0) angleList.push_back(angle->addAngle(Natoms+species->angleindex1[n]-1, Natoms+species->angleindex2[n]-1, Natoms+species->angleindex3[n]-1, species->theta0[n] * RadiansPerDegree, 2.0*species->Ktheta[n] * KJPerKcal));
		if ((species->Kub[n] > 0.0) && (species->S0[n] > 0.0)) {
			if ((species->mass[species->angleindex1[n]-1] < 2.0) && (species->mass[species->angleindex3[n]-1] < 2.0)) constraintList.push_back(system->addConstraint(Natoms+species->angleindex1[n]-1, Natoms+species->angleindex3[n]-1, species->S0[n] * NmPerAngstrom));
			bondList.push_back(bonded->addBond(Natoms+species->angleindex1[n]-1, Natoms+species->angleindex3[n]-1, species->S0[n] * NmPerAngstrom, 2.0*species->Kub[n] * KJPerKcal * AngstromsPerNm * AngstromsPerNm));
		}
	}
	
	vector<double> improperparameters(2);
	list <int> torsionList;
	list <int> improperList;
	for (int n = 0; n < species->NuniqueD; n++) {
		for (int i = 0; i < species->dihedralV[n].N; i++) {
			if (species->dihedralV[n].Kchi[i] != 0.0) {
				torsionList.push_back(torsion[0]->addTorsion(Natoms+species->dihedralindex1[n]-1, Natoms+species->dihedralindex2[n]-1, Natoms+species->dihedralindex3[n]-1, Natoms+species->dihedralindex4[n]-1,
					(int)species->dihedralV[n].nD[i], species->dihedralV[n].delta[i] * RadiansPerDegree, species->dihedralV[n].Kchi[i] * KJPerKcal));
			}
			if (species->dihedralV[n].Kpsi[i] != 0.0) {
				improperparameters[0] = species->dihedralV[n].Kpsi[i] * KJPerKcal;
				improperparameters[1] = species->dihedralV[n].psi0[i] * RadiansPerDegree;
				improperList.push_back(impropertorsion->addTorsion(Natoms+species->dihedralindex1[n]-1, Natoms+species->dihedralindex2[n]-1, Natoms+species->dihedralindex3[n]-1, Natoms+species->dihedralindex4[n]-1, improperparameters));
			}
		}
	}


	list <int> cmapList;
	for (int n = 0; n < species->Ncmap; n++) {
		if (species->cmapassigned[n]) cmapList.push_back(cmaptorsion->addTorsion(species->cmapnumber[n], Natoms+species->cindex[n][0]-1, Natoms+species->cindex[n][1]-1, Natoms+species->cindex[n][2]-1, Natoms+species->cindex[n][3]-1, 
			Natoms+species->cindex[n][4]-1, Natoms+species->cindex[n][5]-1, Natoms+species->cindex[n][6]-1, Natoms+species->cindex[n][7]-1));
	}

	particleIndex.push_back(particleList);
	constraintIndex.push_back(constraintList);
	exceptionIndex.push_back(exceptionList);
	nonbondedIndex.push_back(nonbondedList);
	bondIndex.push_back(bondList);
	angleIndex.push_back(angleList);
	torsionIndex.push_back(torsionList);
	improperIndex.push_back(improperList);
	cmapIndex.push_back(cmapList);
	
	if (rvec.size() != Natoms + species->Natoms) THROWERROR("Sizes do not match");
	for (int i = Natoms; i < Natoms + species->Natoms; i++) {
		x[i] = rvec[i][0] * AngstromsPerNm;
		y[i] = rvec[i][1] * AngstromsPerNm;
		z[i] = rvec[i][2] * AngstromsPerNm;
		mass[i] = species->mass[i - Natoms];
		sigma[i] = species->sigma[i - Natoms];
	}
	Nadded++;
	Natoms += species->Natoms;
	if (rvec.size() != Natoms) THROWERROR("Rvec does not match Natoms");
	
	return;
	
}


void gcmc::deleteMolecule(bool undo) {

	if (Nadded == 0) return;
	cout << "Deleting" << endl;
	State state;

	int index = Nadded-1;
	state = context->getState(State::Positions);
	vector<Vec3> pos = state.getPositions();

	// swap a random species molecule to the end
	if (!undo) {
		index = (int)(Nadded*drand48());

		if (index != Nadded-1) {
			for (int i = 0; i < species->Natoms; i++) {
				Vec3 tmp = pos[baseAtoms + index*species->Natoms + i];
				pos[baseAtoms + index*species->Natoms + i] = pos[baseAtoms + (Nadded-1)*species->Natoms + i];
				pos[baseAtoms + (Nadded-1)*species->Natoms + i] = tmp;
			}
		}
	}
	index = Nadded-1;
	for (int i = 0; i < species->Natoms; i++) {
		oldx[i] = pos[baseAtoms + index*species->Natoms + i][0];
		oldy[i] = pos[baseAtoms + index*species->Natoms + i][1];
		oldz[i] = pos[baseAtoms + index*species->Natoms + i][2];
	}	
	rvec = pos;

	nonbondedIndex[index].reverse();
	for (list<int>::iterator i = nonbondedIndex[index].begin(); i != nonbondedIndex[index].end(); i++) {
		int p = *i;
		nonbonded->deleteParticle(p);
	}
	constraintIndex[index].reverse();
	for (list<int>::iterator i = constraintIndex[index].begin(); i != constraintIndex[index].end(); i++) {
		int p = *i;
		system->deleteConstraint(p);
	}
	exceptionIndex[index].reverse();
	for (list<int>::iterator i = exceptionIndex[index].begin(); i != exceptionIndex[index].end(); i++) {
		int p = *i;
		nonbonded->deleteException(p);
	}
	bondIndex[index].reverse();
	for (list<int>::iterator i = bondIndex[index].begin(); i != bondIndex[index].end(); i++) {
		int p = *i;
		bonded->deleteBond(p);
	}
	angleIndex[index].reverse();
	for (list<int>::iterator i = angleIndex[index].begin(); i != angleIndex[index].end(); i++) {
		int p = *i;
		angle->deleteAngle(p);
	}
	torsionIndex[index].reverse();
	for (list<int>::iterator i = torsionIndex[index].begin(); i != torsionIndex[index].end(); i++) {
		int p = *i;
		torsion[0]->deleteTorsion(p);
	}
	improperIndex[index].reverse();
	for (list<int>::iterator i = improperIndex[index].begin(); i != improperIndex[index].end(); i++) {
		int p = *i;
		impropertorsion->deleteTorsion(p);
	}
	cmapIndex[index].reverse();
	for (list<int>::iterator i = cmapIndex[index].begin(); i != cmapIndex[index].end(); i++) {
		int p = *i;
		cmaptorsion->deleteTorsion(p);
	}
	particleIndex[index].reverse();
	cout << "Deleting system particles" << endl;
	for (list<int>::iterator i = particleIndex[index].begin(); i != particleIndex[index].end(); i++) {
		int p = *i;
		system->deleteParticle(p);
	}
	for (int i = 0; i < species->Natoms; i++) rvec.pop_back();
	particleIndex.pop_back();
	nonbondedIndex.pop_back();
	constraintIndex.pop_back();
	exceptionIndex.pop_back();
	bondIndex.pop_back();
	angleIndex.pop_back();
	torsionIndex.pop_back();
	improperIndex.pop_back();
	cmapIndex.pop_back();
	Nadded--;
	Natoms -= species->Natoms;
	if (rvec.size() != Natoms) THROWERROR("Rvec does not match Natoms");
}


#define MASS_U 1.66054e-27
#define KM_J 6.9477e-21
#define PLANCK 6.62607e-34
#define CUBE(m) (m)*(m)*(m)

void gcmc::runGCMC(int Niterations, int Nsteps) {

	double factor = 0.0;
	double volume = Lx*Ly*Lz;
	double mu = species->mu;
	double beta = 1.6774 * 300.0 / temperature;
	// Lambda in nanometers
	double lambda = 1e9 * PLANCK / sqrt(3.0*species->totalmass*MASS_U*(KM_J/beta));  
	bool inserted;
	bool mademove;
	
	for (int i = 0; i < Niterations; i++) {
		run(Nsteps,1000);
		species->run(Nsteps,1000);
		const State& state = context->getState(State::Energy);
		double Eold = state.getPotentialEnergy() * KcalPerKJ;
		mademove = true;
		inserted = false;
		if (drand48() < 0.5) {
			addMolecule();
			factor = volume*exp(beta*mu)/(CUBE(lambda)*(Nadded));
			inserted = true;
		}
		else {
			cout << "Attempt delete" << endl;
			if (Nadded > 0) {
				deleteMolecule();
				factor = exp(-beta*mu)*CUBE(lambda)*(Nadded+1)/volume;
			} else mademove = false;
			inserted = false;
		}
		if (mademove) {
			renewContext();
			const State& state = context->getState(State::Energy);
			double Enew = state.getPotentialEnergy() * KcalPerKJ;
			double P = factor*exp(-beta*(Enew - Eold));
			if ((P >= 1.0) || (drand48() < P)) {
				cout << "Accepted move with deltaE equal to " << (Enew - Eold) << " kcal/mol" << endl;
			} else {
				cout << "Rejected move with deltaE equal to " << (Enew - Eold) << " kcal/mol" << endl;
				if (inserted) deleteMolecule(true);
				else addMolecule(true);
				renewContext();
			}
		}
		cout << Nadded << " species in the system and " << Natoms << " atoms" << endl;
	}

}
