# ----------------------------------------------------------------------
# Makefile for OpenMM Preview Release 4 workshop "hello world" examples.
# August 18, 2009
# See https://simtk.org/home/openmm.
# ----------------------------------------------------------------------
# This assumes you have gcc compilers for whatever language you are
# using: g++ for C++ and C, gfortran for Fortran 95.
# 
# For the C and Fortran examples, we're depending on your version of
# OpenMM to have been built with the automatically-generated API
# wrappers.
#
# This has had only minimal testing, although it has been known to
# work. It is likely to work fine for C and C++. For Fortran, you
# may need to add some of the C/C++ libraries: 
#    -lc -lm -lstdc++ (or -lstdc++.6) -lgcc -lgcc_s
# but this wasn't required for these examples on Centos 5.2 using
# gcc 4.1.2.

# Check whether this is the right capitalization for your install directory.
OpenMM_INSTALL_DIR=/usr/local/openmm
CFLAGS = -std=c++0x -Iopenmm -Imolecules -g 
FFLAGS = -g -ffree-line-length-none

# This one is to specify the platform, uncomment and specify or it will detect
#DEFINES = -DUSE_PLATFORM="\"OpenCL\""


LIB_DIR=$(OpenMM_INSTALL_DIR)/lib
INCLUDE_DIR=$(OpenMM_INSTALL_DIR)/include -I$(OpenMM_INSTALL_DIR)/include/openmm/internal
LIBS= -lOpenMM

ALL_PROGS = runserialize runclone runclone-resume runmd runhexane runhexane_TI # rungcmc

MYSOURCE = molecules/realmolecule.cpp molecules/charmmmolecule.cpp openmm/molecule.cpp openmm/OBC2Force.cpp openmm/SteepestDescentIntegrator.cpp openmm/AMDIntegrator.cpp openmm/simulation.cpp util/dcdtool.cpp util/misc.cpp util/socket.cpp

GCMCSOURCE = openmm/gcmc.cpp
HEXANESOURCE = openmm/hexane.cpp

default: all

all : $(ALL_PROGS)

# Treat all .cpp source files the same way.
.cpp : 
	g++ $(CFLAGS) -I$(INCLUDE_DIR) -c

runclone: runclone.cpp $(MYSOURCE)
	g++ $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) runclone.cpp $(MYSOURCE) \
	-L$(LIB_DIR) $(LIBS) -o runclone

runclone-resume: runclone-resume.cpp $(MYSOURCE)
	g++ $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) runclone-resume.cpp $(MYSOURCE) \
	-L$(LIB_DIR) $(LIBS) -o runclone-resume


runmd: runmd.cpp $(MYSOURCE)
	g++ $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) runmd.cpp $(MYSOURCE) \
	-L$(LIB_DIR) $(LIBS) -o runmd

runserialize: runserialize.cpp $(MYSOURCE)
	g++ $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) runserialize.cpp $(MYSOURCE) \
	-L$(LIB_DIR) $(LIBS) -o runserialize

rungcmc: rungcmc.cpp $(MYSOURCE) $(GCMCSOURCE)
	g++ $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) rungcmc.cpp $(MYSOURCE) $(GCMCSOURCE) \
	-L$(LIB_DIR) $(LIBS) -o rungcmc

runhexane: runhexane.cpp $(MYSOURCE) $(HEXANESOURCE)
	g++ $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) runhexane.cpp $(MYSOURCE) $(HEXANESOURCE) \
	-L$(LIB_DIR) $(LIBS) -o runhexane

runhexane_TI: runhexane_TI.cpp $(MYSOURCE) $(HEXANESOURCE)
	g++ $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) runhexane_TI.cpp $(MYSOURCE) $(HEXANESOURCE) \
	-L$(LIB_DIR) $(LIBS) -o runhexane_TI

clean : 
	rm $(ALL_PROGS) *.o *.mod *.obj *.exe ; rm -rf *.dSYM; rm -rf ._*

