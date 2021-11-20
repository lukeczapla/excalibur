#include <stdio.h>
#include <time.h>

#include <iostream>

#include "molecules/charmmmolecule.h"
#include "openmm/GCMC.h"

// GCMC simulation code using OpenMM written by Luke Czapla (luke.czapla@yahoo.com)

using namespace OpenMM;


int main(int argc, char **argv) {

	if (argc < 6) {
		printf("This program needs a PSF file, PDB coordinates, and CHARMM-formatted parameter files with the parameters\n\n");
		printf("Usage: runcharmm [psf_filename] [pdb_filename] [species_psf_filename] [species_psf_filename] [parfile1 parfile2 ...]\n\n");
		return -1;
	}

	srand48(time(NULL));

	cout << "Constructing molecule from PSF file " << argv[1] << endl;
	charmmmolecule *charmmmol = new charmmmolecule(argv[1]);
	cout << "Constructing GCMC species from PSF file " << argv[3] << endl;
	charmmmolecule *charmmspecies = new charmmmolecule(argv[3]);
	for (int i = 5; i < argc; i++) {
		cout << "Reading parameter file " << argv[i] << endl;
		charmmmol->loadmoleculePRM(argv[i]);
		charmmspecies->loadmoleculePRM(argv[i]);
	}
	cout << "Reading PDB file " << argv[2] << endl;
	charmmmol->loadPDBcoordinates(argv[2]);
	cout << "Reading species PDB file " << argv[4] << endl;
	charmmspecies->loadPDBcoordinates(argv[4]);
	
	charmmspecies->mu = 20.0;
	
	if ((charmmmol->checkparameters()) && (charmmspecies->checkparameters())) {

		gcmc *mysim = new gcmc(*charmmmol);
		simulation *buffer = new simulation(*charmmspecies);

// here's where the fun starts	
		mysim->visualize = true;
		buffer->visualize = false;
		if (mysim->visualize) mysim->drawmolecule("localhost", 7001);

		cout << "Building OpenMM system from molecule data" << endl;
		mysim->buildSystem();
		cout << "Building OpenMM system from species data" << endl;
		buffer->buildSystem();

		mysim->setTemperature(303.0, 25);
		buffer->setTemperature(303.0, 25);

		mysim->setPBCbox(6.4, 6.4, 6.4, 1.0);
		mysim->setTimestep_fs(2.0);
		buffer->setTimestep_fs(2.0);

		Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
#ifdef USE_PLATFORM
		mysim->setPlatform(USE_PLATFORM);
		buffer->setPlatform(USE_PLATFORM);
#endif
		mysim->setup();
		buffer->setup();
		
		cout << "Minimizing the system" << endl;
		mysim->steepestDescent(1000, 100);
		mysim->minimize(1000, 100);

		cout << "Minimizing the buffer" << endl;
		buffer->steepestDescent(1000, 100);
		buffer->minimize(1000, 100);
		
		mysim->expandcoordinates(100000);
		mysim->setSpecies(buffer);

		mysim->runGCMC(1000, 1000);

		delete mysim;
		delete buffer;
		
	}
	
	delete charmmmol;
	delete charmmspecies;
	return 1;

}

