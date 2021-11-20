#include <string>
#include <sstream>
#include "AMDIntegrator.h"

using namespace OpenMM;
using namespace std;

// alpha and E are in kJ/mol now!

AMDDihedralIntegrator::AMDDihedralIntegrator(int group, double dt, double alphavalue, double Evalue) : Ngroup(group), alpha(alphavalue), E(Evalue), CustomIntegrator(dt) {
	addGlobalVariable("alphaGroup", alpha);
	addGlobalVariable("EGroup", E);
	addGlobalVariable("groupEnergy", 0);
	addPerDofVariable("oldx", 0);
	addPerDofVariable("fg", 0);
	addUpdateContextState();
	addComputeGlobal("groupEnergy", "energy"+to_string(group));
	addComputePerDof("fg", "f"+to_string(group));
	addComputePerDof("v", "v+dt*fprime/m; fprime=fother + fg*((1-modify) + modify*(alphaGroup/(alphaGroup+EGroup-groupEnergy))^2); fother=f-fg; modify=step(EGroup-groupEnergy)");
	addComputePerDof("oldx", "x");
	addComputePerDof("x", "x+dt*v");
	addConstrainPositions();
	addComputePerDof("v", "(x-oldx)/dt");
}


double AMDDihedralIntegrator::getBoostedEnergy(double groupEnergy) {
	
	if (groupEnergy < E) {
		return groupEnergy + (E-groupEnergy)*(E-groupEnergy)/(alpha+E-groupEnergy);
	}

	return groupEnergy;
	
}


