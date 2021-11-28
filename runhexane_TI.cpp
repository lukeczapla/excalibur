#include <stdio.h>
#include <iostream>

#include "molecules/charmmmolecule.h"
#include "openmm/simulation.h"
#include "openmm/hexane.h"

// MD simulation code using OpenMM written by Luke Czapla (czaplaluke@gmail.com)

using namespace OpenMM;

int main(int argc, char **argv) {

	if (argc < 2) {
		printf("This program needs PDB coordinates\n\n");
		printf("Usage: runhexane [pdb_filename]\n\n");
		return -1;
	}

	cout << "Constructing molecule from PDB file " << argv[1] << endl;
	hexane *mysim = new hexane();
	mysim->buildHexane(6, argv[1], 1, 80.0, 80.0, 80.0);
	
	cout << "constructing system" << endl;

// here's where the fun starts	
	// here is the settings for the system, temperature, pressure, periodic box
	mysim->setTemperature(303.0, 25);
//	mysim->setPressure(1.0135);
	mysim->setPBCbox(8.0, 6.4, 6.4, 3.0);
	mysim->setTimestep_fs(5.0);

	mysim->visualize = true;
	cout << "Drawing" << endl;
	if (mysim->visualize) mysim->drawmolecule("localhost", 7001, 80.0, 64.0, 64.0);

	cout << "Setting up" << endl;
		// for maintaining a zero center of mass motion 
//		mysim->maintainCM();

	Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
#ifdef USE_PLATFORM
	mysim->setPlatform(USE_PLATFORM);
#endif
	mysim->setup();
		
	cout << "Minimizing the system" << endl;
//	mysim->steepestDescent(1000, 100);
	mysim->minimize(3000, 100);
		
	cout << "Simulating..." << endl;

		// this is for outputting the data
	mysim->setupDCD("output-run.dcd");
	double Gsolv = 0.0;
	int Ntrials = 20;
	for (int i = 0; i < Ntrials; i++) Gsolv += mysim->score();
	cout << "Estimated solvation free energy = " << Gsolv / (double)Ntrials << endl;

//	mysim->runTI(100000, 1000, 10);
	delete mysim;
		
	
	return 1;

}

