#ifndef _GCMC_H
#define _GCMC_H

#include <list>
#include "OpenMM.h"

#include "simulation.h"
#include "realmolecule.h"

using namespace std;
using namespace molecules;

namespace OpenMM {

	
class gcmc : public simulation {

	simulation *species;
	int Nadded;
	int baseAtoms;
	
	vector < list<int> > particleIndex;
	vector < list<int> > constraintIndex;
	vector < list<int> > exceptionIndex;
	vector < list<int> > nonbondedIndex;
	vector < list<int> > bondIndex;
	vector < list<int> > angleIndex;
	vector < list<int> > torsionIndex;
	vector < list<int> > improperIndex;
	vector < list<int> > cmapIndex;
	
	double *oldx;
	double *oldy;
	double *oldz;

public:
	
	gcmc();
	~gcmc();
	
	gcmc(const molecule& m);

	void setup();

	void addMolecule(bool restore = false);
	void deleteMolecule(bool undo = false);

	void setSpecies(simulation *m);
	
	void runGCMC(int Niterations, int Nsteps);
	
};
	

}

#endif



