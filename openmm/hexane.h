#ifndef _HEXANE_H
#define _HEXANE_H

#include "molecule.h"
#include "OBC2Force.h"
#include "../util/dcdtool.h"

namespace OpenMM {

class hexane : public molecule {

	Integrator* integrator;
	OBC2Force* obc;
	AndersenThermostat* thermostat;
	MonteCarloBarostat* barostat;
	bool platformSpecified;
	string platformName;
		
	double temperature; double temperatureFrequency;
	bool usePressure; double pressure; double pressureFrequency;
	bool usePBCbox;
	double Lx, Ly, Lz, cutoff;
	double timestep_fs;
	
	bool writeDCD;
	int DCDfreq;
	dcdtool *DCD;

public:
	hexane() { platformSpecified = false; }
	void setTemperature(double T, double tfreq = 50) { temperature = T; temperatureFrequency = tfreq; }
	void setPressure(double P, double pfreq = 10) { usePressure = true;  pressure = P;  pressureFrequency = pfreq;  }
	void setTimestep_fs(double ts_fs) { timestep_fs = ts_fs; }
	void setPBCbox(double a, double b, double c, double cut) { usePBCbox = true; Lx = a; Ly = b; Lz = c; cutoff = cut; }
	void setcutoff(double cut) { cutoff = cut; }
	void setPlatform(string pname) { platformSpecified = true; platformName = pname; }
	
	void buildHexane(int M, const char *pdbfile, int N = 1, double Lx = 40.0, double Ly = 40.0, double Lz = 40.0, bool fullPDB = false);
	void buildCyclo(int M, const char *pdbfile, int N = 1, double Lx = 40.0, double Ly = 40.0, double Lz = 40.0, bool fullPDB = false);
	void setupDCD(const char *DCDfilename, int freq = 500);
	
	void setup();
	void reportEnergy();
	
	void renewContext();
	
	void minimize(int Nsteps, int ofreq = 10);
	void steepestDescent(int Nsteps, int ofreq = 10);

	double score();
	void runTI(int Nsteps = 100000, int Navg = 1000, int Nwindows = 10);
	void run(int Nsteps, int ofreq = 100);
	
};

}

#endif
