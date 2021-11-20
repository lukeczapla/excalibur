#ifndef _CHARMMMOLECULE_OPENMM_H
#define _CHARMMMOLECULE_OPENMM_H

#define MAXGROUPS 32

#include <iostream>
#include <fstream>

#include "openmm/serialization/XmlSerializer.h"
#include "OpenMM.h"
#include "realmolecule.h"
#include "OBC2Force.h"

using namespace molecules;
using namespace std;

namespace OpenMM {

	class molecule : public realmolecule {

	protected:

		vector<Vec3> rvec;
		System* system;
		bool systemBuilt;
		Context* context;
		bool hasContext;
		bool usecm;

		vector<int> constraintGroups;
		vector<double> constraintScaling;

		CustomExternalForce** atomconstraints;
		NonbondedForce* nonbonded;
		HarmonicBondForce* bonded;
		HarmonicAngleForce* angle;
		PeriodicTorsionForce** torsion;
		CustomTorsionForce* impropertorsion;
		CMAPTorsionForce* cmaptorsion;
		CMMotionRemover* cmmotion;
		GBSAOBCForce* gbsa;
		GBSAOBCForce* gbvi;
		OBC2Force* obc;
		bool usegbsa, usegbvi;
		int* atomgroups;
		int Ngroups;
		Vec3 boxL;
		
	public:

		molecule() {}
		molecule(const realmolecule& m);

		void serializeSystem() {
			ofstream fout; fout.open("system.xml");
			XmlSerializer::serialize(system, "charmm", fout);
			fout.close();
		}

		void deserializeSystem() {
			ifstream fin; fin.open("system.xml");
			system = XmlSerializer::deserialize<System>(fin);
			fin.close();
		}

		void addConstraintGroup(int index, double scaling = 10.0);
		void removeConstraints(int index);

		void setDihedralGroup(int group, int n);
		int addGroup(vector <int> a);
		void removeGroup(int index);
		
		void maintainCM() { usecm = true; }
		void useGBSA() { usegbsa = true; usegbvi = false; }
		void useGBVI() { usegbsa = false; usegbvi = true; }
		
		double getPotentialEnergy();

		int getNgroups() { return Ngroups; }
		int* getatomgroups() { return atomgroups; }
		
		System* createEmptySystem() { return system = new System(); }
		System* buildSystem(bool verbose = false);
		System* buildSystemClone(int N = 1, bool verbose = false);

		void randomizePosition();
		void setbox(double X, double Y, double Z) { boxL[0] = X; boxL[1] = Y; boxL[2] = Z; }
		System* getSystem() { return system; }
		bool IsSystemBuilt() const { return systemBuilt; }
		void setSystem(System* s) { system = s; }

	};
	
}

#endif
