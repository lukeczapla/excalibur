
#include "hexane.h"
#include "OBC2Force.h"
#include "OpenMM.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

/*
CustomNonbondedForce* gaussian(int N, int *p1, int *p2, double *h, double *r_h, double *sigma) {
	CustomNonbondedForce *nbforce = new CustomNonbondedForce();
}
*/


void hexane::buildCyclo(int M, const char *pdbname, int N, double Lx, double Ly, double Lz, bool fullPDB) {
	Natoms = M * N;
	totalmass = 0.0;
	x = new double[Natoms];
	y = new double[Natoms];
	z = new double[Natoms];
	sigma = new double[Natoms];
	mass = new double[Natoms];
	obc = new OBC2Force();
	nonbonded = new NonbondedForce();
	bonded = new HarmonicBondForce();
	angle = new HarmonicAngleForce();
	torsion = new PeriodicTorsionForce*[1];
	torsion[0] = new PeriodicTorsionForce();

	system = new System();
	usegbvi = false;
	usegbsa = false;
	vector < pair<int,int> > bondPairs;
	for (int i = 0; i < N; i++) {
		cout << "Hexane #" << i << " out of " << "N" << endl;
		double xavg = 0.0; double yavg = 0.0; double zavg = 0.0;
		if (fullPDB) loadPDBmolecule(pdbname, M*N, 0); 
		else {
			loadPDBmolecule(pdbname, M, M*i);
			for (int j = 0; j < M; j++) {
				xavg += x[M*i+j];
				yavg += y[M*i+j];
				zavg += z[M*i+j];
			}
			xavg /= M; yavg /= M; zavg /= M;
		}
		double xr = Lx*drand48();
		double yr = Ly*drand48();
		double zr = Lz*drand48();
		vector < double > parameters(4);
		for (int j = 0; j < M; j++) {
			if (!fullPDB) {
				x[M*i+j] -= xavg;
				y[M*i+j] -= yavg;
				z[M*i+j] -= zavg;
				x[M*i+j] += xr;
				y[M*i+j] += yr;
				z[M*i+j] += zr;
			}
			const Vec3 pos(x[M*i+j] * NmPerAngstrom, y[M*i+j] * NmPerAngstrom, z[M*i+j] * NmPerAngstrom);
			rvec.push_back(pos);
			system->addParticle(mass[M*i+j] = 16.03);
			totalmass += mass[M*i+j];
			nonbonded->addParticle(0.0, 0.400, 0.12 * KJPerKcal);
			sigma[M*i+j] = 4.00;
			parameters[0] = 0.0;
			parameters[1] = 0.24;
			parameters[2] = 1.0;
			parameters[3] = 1.0*KJPerKcal;
			obc->addParticle(parameters);
		}
// bonds
		for (int j = 0; j < M; j++) {
			bonded->addBond(M*i+j, M*i+(j+1)%M, 0.153, 20.0 * KJPerKcal * AngstromsPerNm * AngstromsPerNm);
			bondPairs.push_back(make_pair(M*i+j, M*i+(j+1)%M));
		}
// angles
		for (int j = 0; j < M; j++) {
			angle->addAngle(M*i+j, M*i+(j+1)%M, M*i+(j+2)%M, 112.0 * RadiansPerDegree, 20.0 * KJPerKcal);
		}
// torsions (CHARMM36 lipid)
		for (int j = 0; j < M; j++) {
			torsion[0]->addTorsion(M*i+j, M*i+(j+1)%M, M*i+(j+2)%M, M*i+(j+3)%M, 2, 0.0, 0.162 * KJPerKcal);
			torsion[0]->addTorsion(M*i+j, M*i+(j+1)%M, M*i+(j+2)%M, M*i+(j+3)%M, 3, 180.0 * RadiansPerDegree, 0.047 * KJPerKcal);
			torsion[0]->addTorsion(M*i+j, M*i+(j+1)%M, M*i+(j+2)%M, M*i+(j+3)%M, 4, 0.0, 0.105 * KJPerKcal);
			torsion[0]->addTorsion(M*i+j, M*i+(j+1)%M, M*i+(j+2)%M, M*i+(j+3)%M, 5, 0.0, 0.177 * KJPerKcal);
		}
	}
	nonbonded->createExceptionsFromBonds(bondPairs, 1.0, 1.0);
	system->addForce(nonbonded);
	system->addForce(bonded);
	system->addForce(angle);
	system->addForce(torsion[0]);
	system->addForce(obc);
	obc->setForceGroup(4);
	nonbonded->setForceGroup(4);
	systemBuilt = true;
}



void hexane::buildHexane(int M, const char *pdbname, int N, double Lx, double Ly, double Lz, bool fullPDB) {
	Natoms = M * N;
	totalmass = 0.0;
	x = new double[Natoms];
	y = new double[Natoms];
	z = new double[Natoms];
	sigma = new double[Natoms];
	mass = new double[Natoms];
	obc = new OBC2Force();
	nonbonded = new NonbondedForce();
	bonded = new HarmonicBondForce();
	angle = new HarmonicAngleForce();
	torsion = new PeriodicTorsionForce*[1];
	torsion[0] = new PeriodicTorsionForce();

	system = new System();
	usegbvi = false;
	usegbsa = false;
	vector < pair<int,int> > bondPairs;
	for (int i = 0; i < N; i++) {
		cout << "Hexane #" << i << " out of " << "N" << endl;
		double xavg = 0.0; double yavg = 0.0; double zavg = 0.0;
		if (fullPDB) loadPDBmolecule(pdbname, M*N, 0); 
		else {
			loadPDBmolecule(pdbname, M, M*i);
			for (int j = 0; j < M; j++) {
				xavg += x[M*i+j];
				yavg += y[M*i+j];
				zavg += z[M*i+j];
			}
			xavg /= M; yavg /= M; zavg /= M;
		}
		double xr = Lx*drand48();
		double yr = Ly*drand48();
		double zr = Lz*drand48();
		vector < double > parameters(4);
		for (int j = 0; j < M; j++) {
			if (!fullPDB) {
				x[M*i+j] -= xavg;
				y[M*i+j] -= yavg;
				z[M*i+j] -= zavg;
				x[M*i+j] += xr;
				y[M*i+j] += yr;
				z[M*i+j] += zr;
			}
			const Vec3 pos(x[M*i+j] * NmPerAngstrom, y[M*i+j] * NmPerAngstrom, z[M*i+j] * NmPerAngstrom);
			rvec.push_back(pos);
			if ((j == 0) || (j == M-1)) system->addParticle(mass[M*i+j] = 17.04);
			else system->addParticle(mass[M*i+j] = 16.03);
			totalmass += mass[M*i+j];
			nonbonded->addParticle(0.0, 0.400, 0.12 * KJPerKcal);
			sigma[M*i+j] = 4.00;
			parameters[0] = 0.0;
			parameters[1] = 0.22;
			parameters[2] = 1.0;
			if ((j == 0) || (j == M-1)) parameters[3] = 1.3*KJPerKcal;
			else parameters[3] = 1.0*KJPerKcal;
			obc->addParticle(parameters);
		}
// bonds
		for (int j = 0; j < M-1; j++) {
			bonded->addBond(M*i+j, M*i+j+1, 0.153, 20.0 * KJPerKcal * AngstromsPerNm * AngstromsPerNm);
			bondPairs.push_back(make_pair(M*i+j, M*i+j+1));
		}
// angles
		for (int j = 0; j < M-2; j++) {
			angle->addAngle(M*i+j, M*i+j+1, M*i+j+2, 112.0 * RadiansPerDegree, 20.0 * KJPerKcal);
		}
// torsions (CHARMM36 lipid)
		for (int j = 0; j < M-3; j++) {
			if (j == 1) {
				torsion[0]->addTorsion(M*i+j, M*i+j+1, M*i+j+2, M*i+j+3, 2, 0.0, 0.101 * KJPerKcal);
				torsion[0]->addTorsion(M*i+j, M*i+j+1, M*i+j+2, M*i+j+3, 3, 180.0 * RadiansPerDegree, 0.142 * KJPerKcal);
				torsion[0]->addTorsion(M*i+j, M*i+j+1, M*i+j+2, M*i+j+3, 4, 0.0, 0.074 * KJPerKcal);
				torsion[0]->addTorsion(M*i+j, M*i+j+1, M*i+j+2, M*i+j+3, 5, 0.0, 0.097 * KJPerKcal);
			} else {
				torsion[0]->addTorsion(M*i+j, M*i+j+1, M*i+j+2, M*i+j+3, 2, 0.0, 0.162 * KJPerKcal);
				torsion[0]->addTorsion(M*i+j, M*i+j+1, M*i+j+2, M*i+j+3, 3, 180.0 * RadiansPerDegree, 0.047 * KJPerKcal);
				torsion[0]->addTorsion(M*i+j, M*i+j+1, M*i+j+2, M*i+j+3, 4, 0.0, 0.105 * KJPerKcal);
				torsion[0]->addTorsion(M*i+j, M*i+j+1, M*i+j+2, M*i+j+3, 5, 0.0, 0.177 * KJPerKcal);
			}
		}
	}
	nonbonded->createExceptionsFromBonds(bondPairs, 1.0, 1.0);
	system->addForce(nonbonded);
	system->addForce(bonded);
	system->addForce(angle);
	system->addForce(torsion[0]);
	system->addForce(obc);
	obc->setForceGroup(4);
	nonbonded->setForceGroup(4);
	systemBuilt = true;
}



void hexane::setupDCD(const char *DCDfilename, int freq) {
	writeDCD = true;
	DCDfreq = freq;
	DCD = new dcdtool(Natoms);
	DCD->openfile(DCDfilename);
	DCD->writeDCDheader("OpenMM", DCDfreq, timestep_fs);
}


void hexane::reportEnergy() {
	State state = context->getState(State::Energy, true, 0xFFFFFFFC);
	cout << state.getPotentialEnergy() * KcalPerKJ << " kcal/mol nonbonded" << endl;
}



void hexane::setup() {
	
//	CMMotionRemover &cmmotion = *new CMMotionRemover(1);
	cout << "Adding thermostat and barostat" << endl;
	thermostat = new AndersenThermostat(temperature, temperatureFrequency);
	system->addForce(thermostat);
	if (usePBCbox) {

		nonbonded->setReactionFieldDielectric(1.0);	
		nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
		nonbonded->setCutoffDistance(cutoff);
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
			gbvi->setNonbondedMethod(GBN2Force::CutoffNonPeriodic);
		//	gbvi->setCutoffDistance(cutoff);
		//	gbvi->setSolventDielectric(78.5);
		//	gbvi->setSoluteDielectric(1.0);
		//	gbvi->setBornRadiusScalingMethod(OBC2Force::QuinticSpline);
		//	gbvi->setQuinticLowerLimitFactor(0.4);
		//	gbvi->setQuinticUpperBornRadiusLimit(cutoff);
			nonbonded->setReactionFieldDielectric(1.0);
		}
		nonbonded->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
		nonbonded->setCutoffDistance(cutoff);
	}
	
	integrator = new VerletIntegrator(timestep_fs * PsPerFs);

	cout << "Setting up context" << endl;
    if (platformSpecified) context = new Context(*system, *integrator, Platform::getPlatformByName(platformName));
	else context = new Context(*system, *integrator);
	hasContext = true;

	cout << "Running on platform " << context->getPlatform().getName() << endl;

	cout << "Setting up positions" << endl;
	context->setPositions(rvec);

	State state = context->getState(State::Energy);
	cout << "Potential Energy = " << state.getPotentialEnergy() * KcalPerKJ <<  " kcal/mol" << endl;
	if (visualize) drawmolecule("localhost", 7001);
}
	
	
void hexane::renewContext() {
	context->reinitialize();
	if (rvec.size() != Natoms) THROWERROR("Not enough coordinates defined");
	context->setPositions(rvec);
}


void hexane::steepestDescent(int Nsteps, int ofreq) {
	double stepsize = 0.00005; // nanometers
	double Fmaxtol = 50000.0;
	double F;

	reportEnergy();

	for (int i = 0; i < Nsteps; i++) {
		if (visualize) drawmolecule("localhost", 7001);
		const State& state = context->getState(State::Positions+State::Forces+State::Energy);
		vector<Vec3> forces = state.getForces();
		vector<Vec3> positions = state.getPositions();
		double Fmax = 0.0;
		for (int j = 0; j < positions.size(); j++) {
			F = sqrt(forces[j].dot(forces[j]));
			if (F > Fmax) Fmax = F;
		}
		for (int j = 0; j < positions.size(); j++) {
			forces[j] *= (10.0*stepsize/Fmax);
			positions[j] += forces[j];
		}

		context->setPositions(positions);
		if (i % ofreq == 0) { 
			cout << "Step " << i << " potential energy " << state.getPotentialEnergy() * KcalPerKJ << " kcal/mol, Fmax = " << Fmax * KcalPerKJ << " kcal/mol*nm" << endl;
		}
	}
	
}



void hexane::minimize(int Nsteps, int ofreq) {
	
	int step = 0;
	for (int iter = 0; iter < Nsteps/ofreq; iter++) {
		LocalEnergyMinimizer::minimize(*context, 2, ofreq);
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
	const State& state = context->getState(State::Positions);
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


double hexane::score() {

	vector <double> parameters(4);

        for (int i = 0; i < Natoms; i++) {
		obc->getParticleParameters(i, parameters);
		if ((i == 0) || (i == Natoms-1)) parameters[1] = 0.22;
		else parameters[1] = 0.22;
		parameters[3] = 0.0;
		obc->setParticleParameters(i, parameters);
	}
	obc->updateParametersInContext(*context);

	integrator->step(10000);
	State state = context->getState(State::Energy, true);
	double PE0 = state.getPotentialEnergy();

        for (int i = 0; i < Natoms; i++) {
                obc->getParticleParameters(i, parameters);
                if ((i == 0) || (i == Natoms-1)) parameters[3] = 1.3 * KJPerKcal;
		else parameters[3] = 1.0 * KJPerKcal;
                obc->setParticleParameters(i, parameters);
        }
	obc->updateParametersInContext(*context);

	state = context->getState(State::Positions+State::Energy, true);
	double PE = state.getPotentialEnergy();

	vector<Vec3> pos = state.getPositions();
	for (int i = 0; i < pos.size(); i++) {
		x[i] = pos[i][0] * AngstromsPerNm;                
		y[i] = pos[i][1] * AngstromsPerNm;
		z[i] = pos[i][2] * AngstromsPerNm;
        }
	if (visualize) drawmolecule("localhost", 7001, 10*Lx, 10*Ly, 10*Lz);

	cout << "Solvation free energy = " << (PE-PE0)*KcalPerKJ << " kcal/mol" << endl;

	return (PE-PE0)*KcalPerKJ;

}


void hexane::runTI(int Nsteps, int Navg, int Nwindows) {

	context->setVelocitiesToTemperature(temperature);

	double *PEwindow = new double[Nwindows+1];
	for (int i = 0; i < Nwindows+1; i++) PEwindow[i] = 0.0;

	vector <double> parameters(4);

	for (int i = 0; i < Natoms; i++) {
		obc->getParticleParameters(i, parameters);
		parameters[3] = 0.0;
		obc->setParticleParameters(i, parameters);
	}
	obc->updateParametersInContext(*context);

	for (int iter = 0; iter < Nwindows+1; iter++) {
		cout << "Window " << iter << endl;
		for (int i = 0; i < Natoms; i++) {
			obc->getParticleParameters(i, parameters);
			parameters[3] = 2.0 * KJPerKcal * (double)iter/(double)Nwindows;
			obc->setParticleParameters(i, parameters);
		}
		obc->updateParametersInContext(*context);
		double lambda = (double)iter/(double)Nwindows;
		for (int iter2 = 0; iter2 < Navg; iter2++) {
			integrator->step(Nsteps/Navg);
			const State& state = context->getState(State::Positions+State::Energy, true);
			const vector<Vec3>& pos = state.getPositions();
			double PE = KcalPerKJ * state.getPotentialEnergy();
			cout << "Potential energy = " << PE << " kcal/mol" << endl;
			PEwindow[iter] += PE;
			for (int i = 0; i < pos.size(); i++) {
				x[i] = pos[i][0] * AngstromsPerNm;
				y[i] = pos[i][1] * AngstromsPerNm;
				z[i] = pos[i][2] * AngstromsPerNm;
			}
			if (visualize) drawmolecule("localhost", 7001, 10*Lx, 10*Ly, 10*Lz);
		}
		PEwindow[iter] /= Navg;
		cout << "Window " << iter << ", lambda = " << lambda << " <Potential Energy> = " << PEwindow[iter] << " kcal/mol" << endl;
	}

	for (int i = 0; i < Nwindows+1; i++) cout << "Window " << i << ", lambda = " << (double)i/(double)Nwindows << " <Potential Energy> = " << PEwindow[i] << " kcal/mol" << endl;
	for (int i = 0; i < Nwindows+1; i++) cout << (double)i/(double)Nwindows << " " << PEwindow[i];
}


void hexane::run(int Nsteps, int ofreq) {

	context->setVelocitiesToTemperature(temperature);

	if (writeDCD) DCD->writeframe(x, y, z);
	int step = 0;
	for (int iter = 0; iter < Nsteps/ofreq; iter++) {

		integrator->step(ofreq);
		step += ofreq;
		const State& state = context->getState(State::Positions+State::Energy+State::Forces+State::Parameters, true);
		

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
			double dotp = 0.0;
			int nheavy = 0;
			for (int i = 0; i < pos.size(); i++) {
				if (mass[i] > 2.0) {
					++nheavy;
					dotp += 1e3 * pos[i].dot(forces[i]) / 6.022e23;
				}
			}
			if (visualize) drawmolecule("localhost", 7001, 10*a[0], 10*b[1], 10*c[2]);
//			if (usePressure) cout << "Instantaneous pressure = " << 1e-5 * ((nheavy * 1.3806488e-23 * temperature) + (dotp / 3.0))/(a[0]*b[1]*c[2]*1e-27) << " bar" << endl;
		}

		reportEnergy();
		
		if (writeDCD && (step % DCDfreq == 0)) DCD->writeframe(x, y, z);
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

