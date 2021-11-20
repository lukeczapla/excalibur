
The Sword in the Stone was locked in place for 8 years until Tabitha the Cat came around to weild..... EXCALIBUR!


Standalone MD simulation run SYNTAX:

Load your PSF and PDB file and all the parameter files needed to grab the
parameters, XPLOR PSF format is used on PSF (atom type names, not ID numbers)

"./runmd my_system.psf my_system.pdb par_file1.prm par_file2 par_file3"

runmd.cpp has all the setting and steps defined before the "run()" method

Set up your temperature, pressure, box vectors (if PBC), cutoffs, etc.


Have fun and enjoy!  This is meant for development purposes and running systems constructed in other codes.  The methodology with Volumetric integration and SASA will be adapted for GBn (to do).


*** Compilation of Visualizer and Visualization on TCP/IP port 7001:

Type either "make linux" or "make osx" in the visualizer folder.
it would run in a separate terminal left open, just type "./visualizer 7001"


*** XML Serialization of CHARMM output

"./runserialize system.psf system.pdb par_file1 par_file2 ..."

The System object will be read from system.xml at the moment.



