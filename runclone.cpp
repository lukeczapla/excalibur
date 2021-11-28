#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "molecules/charmmmolecule.h"
#include "openmm/simulation.h"

// MD simulation code using OpenMM written by Luke Czapla (luke.czapla@yahoo.com)

using namespace OpenMM;

int main(int argc, char **argv) {

	if (argc < 5) {
		printf("Runs a simulation with N clones of the molecule data placed randomly\n");
		printf("This program needs a PSF file, PDB coordinates, and CHARMM-formatted parameter files with the parameters\n\n");
		printf("Usage: runclone [N] [psf_filename] [pdb_filename] [parfile1 parfile2 ...]\n\n");
		return -1;
	}

	srand48(time(0));

	cout << "Constructing molecule from PSF file " << argv[2] << endl;
	charmmmolecule *charmmmol = new charmmmolecule(argv[2]);
	for (int i = 4; i < argc; i++) {
		cout << "Reading parameter file " << argv[i] << endl;
		charmmmol->loadmoleculePRM(argv[i]);
	}
	cout << "Reading PDB file " << argv[3] << endl;
	charmmmol->loadPDBcoordinates(argv[3]);

	if (charmmmol->checkparameters()) {

		simulation *mysim = new simulation(*charmmmol);

// here's where the fun starts	
		mysim->visualize = true;
		if (mysim->visualize) mysim->drawmolecule("localhost", 7001);

		cout << "Building OpenMM system from molecule data" << endl;

		// These are Generalized Born parameters that are system-specific
		// for the implicit solvent model (water represented as bulk dielectric)
		mysim->loadGBradii(fopen("gbradii", "r"));
		// volume coefficients
		mysim->loadgamma(fopen("gbvigamma", "r"));
		// SASA coefficients
		mysim->loadgamma2(fopen("gbvigamma2", "r"));
		mysim->useGBSA();

		// This is required to construct the system clones
		mysim->setbox(96.0, 96.0, 96.0);
		mysim->buildSystemClone(atoi(argv[1]));

		// TO SETUP THE aMD BOOST POTENTIAL for enhanced sampling
//		mysim->setDihedralGroup(0, 1);
//		mysim->dihedralBoostSystem(1, 90.0, 360.0);
		
		// here is the settings for the system, temperature, pressure, periodic box
		mysim->setTemperature(298.0, 10);
//		mysim->setPressure(1.0135, 10);
		mysim->setPBCbox(9.6, 9.6, 9.6, 3.0);
		mysim->setTimestep_fs(3.0);

		// for maintaining a zero center of mass motion 
//		mysim->maintainCM();

		Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
#ifdef USE_PLATFORM
		mysim->setPlatform(USE_PLATFORM);
#endif
		mysim->setup();
		std::vector <std::string> properties = mysim->getPlatform().getPropertyNames();
		mysim->setPlatformProperty("OpenCLPrecision", "double");
		for (int i = 0; i < properties.size(); i++) {
			cout << properties[i] << ": " << mysim->getPlatform().getPropertyValue(*mysim->getContext(), properties[i]) << endl;
		}

		cout << "Serializing" << endl;
		mysim->serializeSystem();
		//exit(0);
		
		cout << "Minimizing the system" << endl;
		mysim->steepestDescent(1000, 100);
		mysim->minimize(2500, 100);
		mysim->steepestDescent(1000, 100);
		mysim->writePDBmolecule("systemrun.pdb", mysim->Natoms);
		cout << "Simulating..." << endl;

		// this is for outputting the data
		mysim->setupDCD("output-run.dcd");
		mysim->run(100000000, 2500); // 300 ns simulation
		mysim->writePDBmolecule("endrun.pdb", mysim->Natoms);
		delete mysim;	
	}
	
	delete charmmmol;
	return 1;

}


