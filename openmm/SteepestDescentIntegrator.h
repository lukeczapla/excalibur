#ifndef _STEEPESTDESCENTINTEGRATOR
#define _STEEPESTDESCENTINTEGRATOR

#include "OpenMM.h"

namespace OpenMM {


class SteepestDescentIntegrator : public CustomIntegrator {

	double initialStepSize;

public:
	
	SteepestDescentIntegrator(double initial = 0.001);
	
};


}

#endif