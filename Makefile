# Note you may also need something like this ``` export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/lczapla/miniconda3/lib ```
# Just edit me here, put your plugins in the plugins subfolder here, and good to go!
HOME=/Users/lczapla

# For Mac, you may want to brew install g++ for a real g++ compiler and replace this with g++-11 or whatever version you get
# You could also just drop the -fopenmp flag but if you want to use it, it does work on Mac!
COMPILER=g++

# Check whether this is the right capitalization for your install directory.
OPENMM_INSTALL_DIR=$(HOME)/miniconda3
LIB_DIR=$(OPENMM_INSTALL_DIR)/lib
INCLUDE_DIR=$(OPENMM_INSTALL_DIR)/include

CFLAGS = -std=c++0x -Iopenmm -Imolecules -g -DCL_TARGET_OPENCL_VERSION=220 # -fopenmp

# This one WAS to specify the platform, uncomment and specify or it will use OpenCL
DEFINES = -DUSE_PLATFORM="\"OpenCL\""

LIBS=-lOpenMM

ALL_PROGS = runserialize runclone runclone-resume runmd runhexane runhexane_TI # rungcmc

MYSOURCE = molecules/realmolecule.cpp molecules/charmmmolecule.cpp openmm/molecule.cpp openmm/OBC2Force.cpp openmm/SteepestDescentIntegrator.cpp openmm/AMDIntegrator.cpp openmm/simulation.cpp util/dcdtool.cpp util/misc.cpp util/socket.cpp

GCMCSOURCE = openmm/gcmc.cpp
HEXANESOURCE = openmm/hexane.cpp

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

# runmd was the only one needed to run the proteins (e.g. Abl1 kinase in H2O with KCl)
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
