#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "molecules/charmmmolecule.h"
#include "openmm/simulation.h"

// MD simulation code using OpenMM written by Luke Czapla (luke.czapla@yahoo.com)
// Serialization will work with CHARMM System object (Nhp6A test included)


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

		cout << "Building OpenMM system from XML data" << endl;

// DESERIALIZE THE XML SYSTEM OBJECT
		mysim->deserializeSystem();

		// here is the settings for the system, temperature, pressure, periodic box

		Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
#ifdef USE_PLATFORM
		mysim->setPlatform(USE_PLATFORM);
		mysim->setPlatformProperty("CudaDeviceIndex", "0,1");
#endif
		mysim->setTimestep_fs(2.0);
		mysim->setupXML();
		
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


