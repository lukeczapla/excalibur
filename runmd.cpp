#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "molecules/charmmmolecule.h"
#include "openmm/simulation.h"

// MD simulation code using OpenMM written by Luke Czapla (luke.czapla@yahoo.com)

using namespace OpenMM;

int main(int argc, char **argv) {

	if (argc < 4) {
		printf("This program needs a PSF file, PDB coordinates, and CHARMM-formatted parameter files with the parameters\n\n");
		printf("Usage: runcharmm [psf_filename] [pdb_filename] [parfile1 parfile2 ...]\n\n");
		return -1;
	}

	cout << "Constructing molecule from PSF file " << argv[1] << endl;
	charmmmolecule *charmmmol = new charmmmolecule(argv[1]);
	for (int i = 3; i < argc; i++) {
		cout << "Reading parameter file " << argv[i] << endl;
		charmmmol->loadmoleculePRM(argv[i]);
	}
	cout << "Reading PDB file " << argv[2] << endl;
	charmmmol->loadPDBcoordinates(argv[2]);

	if (charmmmol->checkparameters()) {

		simulation *mysim = new simulation(*charmmmol);

// here's where the fun starts	
		mysim->visualize = true;
		if (mysim->visualize) mysim->drawmolecule("localhost", 7001);

		cout << "Building OpenMM system from molecule data" << endl;

		// These are Generalized Born parameters that are system-specific
		// for the implicit solvent model (water represented as bulk dielectric)
//		mysim->loadGBradii(fopen("gbradii", "r"));
//		mysim->useGBVI();

		// This is required to construct the system
		mysim->buildSystem();

		// TO SETUP THE aMD BOOST POTENTIAL for enhanced sampling
//		mysim->setDihedralGroup(0, 1);
//		mysim->dihedralBoostSystem(1, 90.0, 360.0);
		
		// here is the settings for the system, temperature, pressure, periodic box
		mysim->setTemperature(310.0, 25);
//		mysim->setPressure(1.0135, 10);
		mysim->setPBCbox(6.4, 6.4, 6.4, 2.0);
		mysim->setTimestep_fs(2.0);

		// for maintaining a zero center of mass motion 
//		mysim->maintainCM();

		Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
#ifdef USE_PLATFORM
		mysim->setPlatform(USE_PLATFORM);
		mysim->setPlatformProperty("CudaDeviceIndex", "0,1");
#endif
		mysim->setup();
		mysim->serializeSystem();
		
		cout << "Minimizing the system" << endl;
		mysim->steepestDescent(200, 100);
		mysim->minimize(1000, 100);
		
		cout << "Simulating..." << endl;

		// this is for outputting the data
		mysim->setupDCD("output-run.dcd");
		mysim->run(2500000, 100);
		delete mysim;
		
	}
	
	delete charmmmol;
	return 1;

}


