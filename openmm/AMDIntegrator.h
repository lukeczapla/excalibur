#ifndef _AMDINTEGRATOR_H
#define _AMDINTEGRATOR_H

#include "OpenMM.h"

namespace OpenMM {

// alpha and E are in kJ/mol

class AMDDihedralIntegrator : public CustomIntegrator {

	double alpha, E;
	int Ngroup;

public:
	
	AMDDihedralIntegrator(int group, double dt, double alpha, double E);
	
	inline double getAlpha() {  return alpha;  }
	inline void setAlpha(double alphavalue) {  alpha = alphavalue;  }
	
	double getBoostedEnergy(double groupEnergy);
	double getDeltaV();
	inline double getE() {  return E;  }
	inline void setE(double Evalue) {  E = Evalue;  }

	
};


}


#endif
