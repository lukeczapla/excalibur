#ifndef _CHARMMMOLECULE_H
#define _CHARMMMOLECULE_H

#include "OpenMM.h"
#include "realmolecule.h"

#define CMAPMAXN 20  // maximum numbers of CMAP maps, CHARMM36 has 6 as of 2014


namespace molecules {

class charmmmolecule : public realmolecule {

protected: 
	
	groupname *segname;
	
	bool *isspecificD;
	integer Nproper, Nimproper;

public:
	
	charmmmolecule();
	charmmmolecule(const char *psfname);
	charmmmolecule(const char *psfname, const char *parname);
	charmmmolecule(const charmmmolecule& m);
	
	~charmmmolecule();


	void loadmoleculePSF(const char *psfname, bool verbose = false);
	void loadmoleculePRM(const char *parname, bool verbose = false);

	std::vector <int> atomtypeatoms(const char *atname);
	std::vector <int> segmentresidueatoms(const char *segname, integer res);
	std::vector <int> segmentatoms(const char *segname);
	// getNsegments gets list of starting index of segments, with list size equal to number of segments
	std::vector <int> getNsegments();
	
};


}

#endif

