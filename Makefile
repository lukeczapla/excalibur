
# Check whether this is the right capitalization for your install directory.
OPENMM_INSTALL_DIR=/home/luke/miniconda3
LIB_DIR=/home/luke/miniconda3/lib
INCLUDE_DIR=/home/luke/miniconda3/include

CFLAGS = -std=c++0x -Iopenmm -Imolecules -g -fopenmp -DCL_TARGET_OPENCL_VERSION=220

# DEFUNCT! This one WAS to specify the platform, uncomment and specify or it will use OpenCL
DEFINES = -DUSE_PLATFORM="\"OpenCL\""

LIBS=-lOpenMM

ALL_PROGS = runserialize runclone runclone-resume runmd runhexane runhexane_TI # rungcmc

MYSOURCE = molecules/realmolecule.cpp molecules/charmmmolecule.cpp openmm/molecule.cpp openmm/OBC2Force.cpp openmm/SteepestDescentIntegrator.cpp openmm/AMDIntegrator.cpp openmm/simulation.cpp util/dcdtool.cpp util/misc.cpp util/socket.cpp

GCMCSOURCE = openmm/gcmc.cpp
HEXANESOURCE = openmm/hexane.cpp

COMPILER=g++

default: all

all : $(ALL_PROGS)

# Treat all .cpp source files the same way.
.cpp :
	$(COMPILER) $(CFLAGS) -I$(INCLUDE_DIR) -c

runclone: runclone.cpp $(MYSOURCE)
	$(COMPILER) $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) runclone.cpp $(MYSOURCE) \
	-L$(LIB_DIR) $(LIBS) -o runclone

runclone-resume: runclone-resume.cpp $(MYSOURCE)
	$(COMPILER) $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) runclone-resume.cpp $(MYSOURCE) \
	-L$(LIB_DIR) $(LIBS) -o runclone-resume


runmd: runmd.cpp $(MYSOURCE)
	$(COMPILER) $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) runmd.cpp $(MYSOURCE) \
	-L$(LIB_DIR) $(LIBS) -o runmd

runserialize: runserialize.cpp $(MYSOURCE)
	$(COMPILER) $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) runserialize.cpp $(MYSOURCE) \
	-L$(LIB_DIR) $(LIBS) -o runserialize

rungcmc: rungcmc.cpp $(MYSOURCE) $(GCMCSOURCE)
	$(COMPILER) $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) rungcmc.cpp $(MYSOURCE) $(GCMCSOURCE) \
	-L$(LIB_DIR) $(LIBS) -o rungcmc

runhexane: runhexane.cpp $(MYSOURCE) $(HEXANESOURCE)
	$(COMPILER) $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) runhexane.cpp $(MYSOURCE) $(HEXANESOURCE) \
	-L$(LIB_DIR) $(LIBS) -o runhexane

runhexane_TI: runhexane_TI.cpp $(MYSOURCE) $(HEXANESOURCE)
	$(COMPILER) $(CFLAGS) $(DEFINES) -I$(INCLUDE_DIR) runhexane_TI.cpp $(MYSOURCE) $(HEXANESOURCE) \
	-L$(LIB_DIR) $(LIBS) -o runhexane_TI

clean :
	rm $(ALL_PROGS) *.o *.mod *.obj *.exe ; rm -rf *.dSYM; rm -rf ._*
