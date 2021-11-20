#ifndef _NPTSIMULATION_H
#define _NPTSIMULATION_H

#include "OpenMM.h"
#include "molecule.h"
#include "hexane.h"
#include "../util/dcdtool.h"

using namespace molecules;

namespace OpenMM {
	
	class simulation : public molecule {

	protected: 
		
		Integrator* integrator;
		bool platformSpecified;
		string platformName;
		
		AndersenThermostat* thermostat;
		MonteCarloBarostat* barostat;
		
		// AMD parameters
		bool useBoost; double alpha; double E; int boostGroup;
		
		double temperature; double temperatureFrequency;
		bool usePressure; double pressure; double pressureFrequency;
		bool usePBCbox;
		double Lx, Ly, Lz, cutoff;
		double timestep_fs;
		
		bool writeDCD;
		int DCDfreq;
		dcdtool *DCD;

	public:

		simulation(const molecule& m);

		void setTemperature(double T, double tfreq = 50) { temperature = T; temperatureFrequency = tfreq; }
		void setPressure(double P, double pfreq = 10) { usePressure = true;  pressure = P;  pressureFrequency = pfreq;  }
		void setTimestep_fs(double ts_fs) { timestep_fs = ts_fs; }
		void setPBCbox(double a, double b, double c, double cut) { usePBCbox = true; Lx = a; Ly = b; Lz = c; cutoff = cut; }
		void setcutoff(double cut) { cutoff = cut; }
		Context *getContext() { return context; }
		void setPlatformProperty(std::string property, std::string value) { context->getPlatform().setPropertyValue(*context, property, value); }
		Platform &getPlatform() { return context->getPlatform(); }
		void setPlatform(string pname) { platformSpecified = true; platformName = pname; }
		
		void dihedralBoostSystem(int group, double alpha, double E);
		
		void setupDCD(const char *DCDfilename, int freq = 500);
		
		void setup();
		void setupXML();
		void reportEnergy();
		
		void renewContext();
		
		void minimize(int Nsteps, int ofreq = 10);
		void steepestDescent(int Nsteps, int ofreq = 10);
		void run(int Nsteps, int ofreq = 100);

	};
	
}


#endif

