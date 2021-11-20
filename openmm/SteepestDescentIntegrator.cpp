#include "SteepestDescentIntegrator.h"

using namespace OpenMM;

SteepestDescentIntegrator::SteepestDescentIntegrator(double initial) : CustomIntegrator(1.0) {

	initialStepSize = initial;
	addGlobalVariable("step_size", initialStepSize);
	addGlobalVariable("energy_old", 0);
	addGlobalVariable("energy_new", 0);
	addGlobalVariable("delta_energy", 0);
	addGlobalVariable("accept", 0);
	addGlobalVariable("fnorm2", 0);
	addPerDofVariable("x_old", 0);

	addComputePerDof("x_old", "x");

	addComputeSum("fnorm2", "f^2");

	addComputePerDof("x", "x+step_size*f/sqrt(fnorm2 + delta(fnorm2))");
	addConstrainPositions();

	addComputeGlobal("energy_new", "energy");
	addComputeGlobal("delta_energy", "energy_new-energy_old");
	addComputeGlobal("accept", "step(-delta_energy) * delta(energy - energy_new)");

	addComputePerDof("x", "accept*x + (1-accept)*x_old");

	addComputeGlobal("step_size", "step_size * (2.0*accept + 0.5*(1-accept))");

}

