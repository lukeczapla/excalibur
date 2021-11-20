#ifndef _OBC2FORCE_H
#define _OBC2FORCE_H

#include "OpenMM.h"

namespace OpenMM {

class OBC2Force : public CustomGBForce {
public:
	OBC2Force(double solventDielectric = 78.5, double soluteDielectric = 1.0, double cutoff = 3.0, double kappa = 0.0);
};

class GBN2Force : public CustomGBForce {
public:
	GBN2Force(double solventDielectric = 78.5, double soluteDielectric = 1.0, double cutoff = 3.0, double kappa = 0.0);
};


}

#endif

