#include <stdio.h>
#include <string.h>
#include <iostream>

#include "simulation.h"
#include "AMDIntegrator.h"
#include "OBC2Force.h"

using namespace OpenMM;
using namespace std;


simulation::simulation(const molecule& m) : temperature(300.0), cutoff(1.6), timestep_fs(1.0), usePBCbox(false), useBoost(false), usePressure(false),
	pressureFrequency(5), temperatureFrequency(50), writeDCD(false), platformSpecified(false), platformName(""), molecule(m) {

}


void simulation::dihedralBoostSystem(int group, double alphavalue, double Evalue) {
	useBoost = true;
	alpha = alphavalue;
	E = Evalue;
	boostGroup = group;
}


void simulation::setupDCD(const char *DCDfilename, int freq) {
	writeDCD = true;
	DCDfreq = freq;
	DCD = new dcdtool(Natoms);
	DCD->openfile(DCDfilename);
	DCD->writeDCDheader("Written by MC-OpenMM", DCDfreq, timestep_fs);
}


void simulation::reportEnergy() {
	State state = context->getState(State::Energy, true, 0x00000001);
	cout << state.getPotentialEnergy() * KcalPerKJ << " kcal/mol" << endl;
}



void simulation::setup() {
	
//	CMMotionRemover &cmmotion = *new CMMotionRemover(1);
	cout << "Adding thermostat and barostat" << endl;
	thermostat = new AndersenThermostat(temperature, temperatureFrequency);
	system->addForce(thermostat);
	if (usePBCbox) {
		if (usegbsa) {
			gbsa->setNonbondedMethod(GBSAOBCForce::CutoffPeriodic);
			gbsa->setCutoffDistance(cutoff);
			gbsa->setSolventDielectric(78.5);
			gbsa->setSoluteDielectric(1.0);
			nonbonded->setReactionFieldDielectric(1.0);
		} 
		if (usegbvi) {
			obc->setNonbondedMethod(CustomGBForce::CutoffPeriodic);
			obc->setCutoffDistance(cutoff);
			obc->setGlobalParameterDefaultValue(3, cutoff);
			nonbonded->setReactionFieldDielectric(1.0);
		}
		nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
		nonbonded->setCutoffDistance(cutoff);
//		nonbonded->setEwaldErrorTolerance(1e-3);
//		nonbonded->setUseDispersionCorrection(true);
		system->setDefaultPeriodicBoxVectors(Vec3(Lx,0,0), Vec3(0,Ly,0), Vec3(0,0,Lz));
		cout << "Box vectors (nm) " << Lx << ' ' << Ly << ' ' << Lz << endl;
		if (usePressure) {
			barostat = new MonteCarloBarostat(pressure, temperature, pressureFrequency);
			system->addForce(barostat);
		}
	} else {
		if (usegbsa) {
			gbsa->setNonbondedMethod(GBSAOBCForce::CutoffNonPeriodic);
			gbsa->setCutoffDistance(cutoff);
			gbsa->setSolventDielectric(78.5);
			gbsa->setSoluteDielectric(1.0);
			nonbonded->setReactionFieldDielectric(1.0);
		}
		if (usegbvi) {
			obc->setNonbondedMethod(CustomGBForce::CutoffNonPeriodic);
			obc->setCutoffDistance(cutoff);
			obc->setGlobalParameterDefaultValue(3, cutoff);
			nonbonded->setReactionFieldDielectric(1.0);
		}
		nonbonded->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
		nonbonded->setCutoffDistance(cutoff);
	}
	
	if (useBoost) {
		cout << "Enabling AMD dihedral boost potential" << endl;
		integrator = new AMDDihedralIntegrator(boostGroup, timestep_fs * PsPerFs, alpha, E);
	} 
	else integrator = new VerletIntegrator(timestep_fs * PsPerFs);

	cout << "Setting up context" << endl;
    if (platformSpecified) context = new Context(*system, *integrator, Platform::getPlatformByName(platformName));
	else context = new Context(*system, *integrator);
	hasContext = true;

	cout << "Running on platform " << context->getPlatform().getName() << endl;

	cout << "Setting up positions" << endl;
    context->setPositions(rvec);
	for (int i = 0; i < Natoms; i++) {
		x[i] = rvec[i][0] * AngstromsPerNm;
		y[i] = rvec[i][1] * AngstromsPerNm;
		z[i] = rvec[i][2] * AngstromsPerNm;
	}

	State state = context->getState(State::Energy);
	cout << "Potential Energy = " << state.getPotentialEnergy() * KcalPerKJ <<  " kcal/mol" << endl;
	if (visualize) drawmolecule("localhost", 7001);
}

void simulation::setupXML() {
	integrator = new VerletIntegrator(timestep_fs * PsPerFs);
	if (platformSpecified) context = new Context(*system, *integrator, Platform::getPlatformByName(platformName));
        else context = new Context(*system, *integrator);
        hasContext = true;

        cout << "Running on platform " << context->getPlatform().getName() << endl;

        cout << "Setting up positions" << endl;
        for (int i = 0; i < Natoms; i++) {
                const Vec3 pos(x[i] * NmPerAngstrom, y[i] * NmPerAngstrom, z[i] * NmPerAngstrom);
                rvec.push_back(pos);
        }
	context->setPositions(rvec);
	State state = context->getState(State::Energy);
        cout << "Potential Energy = " << state.getPotentialEnergy() * KcalPerKJ <<  " kcal/mol" << endl;
        if (visualize) drawmolecule("localhost", 7001);
}
	
	
void simulation::renewContext() {
	context->reinitialize();
	if (rvec.size() != Natoms) THROWERROR("Not enough coordinates defined");
	context->setPositions(rvec);
}


void simulation::steepestDescent(int Nsteps, int ofreq) {
	
	double stepsize = 0.0005; // nanometers
	double Fmaxtol = 50000.0;
	double F;

	reportEnergy();

	for (int i = 0; i < Nsteps; i++) {
		if (visualize) drawmolecule("localhost", 7001);
		const State& state = context->getState(State::Positions+State::Forces+State::Energy);
		vector<Vec3> forces = state.getForces();
		vector<Vec3> positions = state.getPositions();
		double Fmax = 0.1;
		for (int j = 0; j < positions.size(); j++) {
			F = sqrt(forces[j].dot(forces[j]));
			if (F > Fmax) Fmax = F;
		}
		for (int j = 0; j < positions.size(); j++) {
			forces[j] /= Fmax;
			positions[j] += forces[j]*stepsize;
		}

		context->setPositions(positions);
		if (i % ofreq == 0) { 
			cout << "Step " << i << " potential energy " << state.getPotentialEnergy() * KcalPerKJ << " kcal/mol, Fmax = " << Fmax * KcalPerKJ << " kcal/mol*nm" << endl;
		}
	}		
	context->applyConstraints(0.0001);
	
}



void simulation::minimize(int Nsteps, int ofreq) {
	
	int step = 0;
	for (int iter = 0; iter < Nsteps/ofreq; iter++) {
		LocalEnergyMinimizer::minimize(*context, 0.1, ofreq);
		step += ofreq;
		const State& state = context->getState(State::Positions+State::Energy);
		vector<Vec3> pos = state.getPositions();
		for (int i = 0; i < pos.size(); i++) {
			if (i >= Natoms) THROWERROR("Out of range atom assigned during minimization");
			x[i] = pos[i][0] * AngstromsPerNm;
			y[i] = pos[i][1] * AngstromsPerNm;
			z[i] = pos[i][2] * AngstromsPerNm;
		}
		if (visualize && usePBCbox) drawmolecule("localhost", 7001, 10*Lx, 10*Ly, 10*Lz, true);
		if (visualize && !usePBCbox) drawmolecule("localhost", 7001);
		cout << "Minimization Step " << step << ", Potential Energy = " << state.getPotentialEnergy() * KcalPerKJ << " kcal/mol" << endl;

	}
	if (Nsteps % ofreq != 0) LocalEnergyMinimizer::minimize(*context, 1.0, Nsteps % ofreq);
	const State& state = context->getState(State::Positions, true);
	rvec = state.getPositions();
	for (int i = 0; i < Natoms; i++) {
		x[i] = rvec[i][0] * AngstromsPerNm;
		y[i] = rvec[i][1] * AngstromsPerNm;
		z[i] = rvec[i][2] * AngstromsPerNm;
	}
	if (visualize && usePBCbox) drawmolecule("localhost", 7001, 10*Lx, 10*Ly, 10*Lz, true);
	if (visualize && !usePBCbox) drawmolecule("localhost", 7001);

	return;
}


void simulation::run(int Nsteps, int ofreq) {

	context->setVelocitiesToTemperature(temperature);

	if (writeDCD) DCD->writeframe(x, y, z);
	int step = 0;
	for (int iter = 0; iter < Nsteps/ofreq; iter++) {

		integrator->step(ofreq);
		step += ofreq;
		const State& state = context->getState(State::Positions+State::Energy+State::Forces+State::Parameters, true);
		
		if (useBoost) {
			
		}

		cout << "Step " << step << " (" << step*timestep_fs << " fs) " << "Potential Energy " << state.getPotentialEnergy() * KcalPerKJ << " kcal/mol" << endl;
		cout << "Kinetic Energy " << state.getKineticEnergy() * KcalPerKJ << " kcal/mol at temperature " << temperature << " K" << endl;

		const vector<Vec3>& pos = state.getPositions();
		const vector<Vec3>& forces = state.getForces();
		const map<string, double>& parameters = state.getParameters();

		for (int i = 0; i < pos.size(); i++) {
			if (i >= Natoms) THROWERROR("Out of range atom assigned during MD");
			x[i] = pos[i][0] * AngstromsPerNm;
			y[i] = pos[i][1] * AngstromsPerNm;
			z[i] = pos[i][2] * AngstromsPerNm;
		}

		if (usePBCbox) {
			Vec3 a, b, c;
			state.getPeriodicBoxVectors(a, b, c);
			cout << "Box vectors (nm): " << a[0] << ' ' << b[1] << ' ' << c[2] << endl;
			cout << "Density: " << totalmass/(602.2*a[0]*b[1]*c[2]) << " g/mL" << endl;
//			double dotp = 0.0;
//			int nheavy = 0;
//			for (int i = 0; i < pos.size(); i++) {
//				if (mass[i] > 2.0) {
//					++nheavy;
//					dotp += 1e3 * pos[i].dot(forces[i]) / 6.022e23;
//				}
//			}
//			if (usePressure) cout << "Instantaneous pressure = " << 1e-5 * ((nheavy * 1.3806488e-23 * temperature) + (dotp / 3.0))/(a[0]*b[1]*c[2]*1e-27) << " bar" << endl;
		}

		reportEnergy();
		
		if (writeDCD && ((step % DCDfreq) == 0)) DCD->writeframe(x, y, z);
		if (visualize && usePBCbox) drawmolecule("localhost", 7001, 10*Lx, 10*Ly, 10*Lz);
		if (visualize && !usePBCbox) drawmolecule("localhost", 7001);

	}
	
	if (Nsteps % ofreq != 0) integrator->step(Nsteps % ofreq);

	const State& state = context->getState(State::Positions);
	rvec = state.getPositions();
	for (int i = 0; i < Natoms; i++) {
		x[i] = rvec[i][0] * AngstromsPerNm;
		y[i] = rvec[i][1] * AngstromsPerNm;
		z[i] = rvec[i][2] * AngstromsPerNm;
	}
	if (visualize && usePBCbox) drawmolecule("localhost", 7001, 10*Lx, 10*Ly, 10*Lz);
	if (visualize && !usePBCbox) drawmolecule("localhost", 7001);

	return;

}

