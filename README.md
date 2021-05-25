# Charge-Driven-Import-of-Bacteriocins

## CPPTraj
prerequisite:
To run the scripts cpptraj from AmberTools must be installed on the local machine

### RMSD
This script uses cpptraj to calculate the RMSD of a protein from a trajectory. 
The topology file is loaded using the command 'parm'. 
Each of the trajectory output files are loaded using the command 'trajin'. 
The command 'rms' is used to specify that the RMSD is to be calculated. 
A mask is created which specifies the residues and atoms to include within the calculation.
The mask can be changed to accomodate labile or non labile residues
An output file is create using the command 'out' followed by a file name. 
The 'norotate' command prevents co-ordinates from rotating which may otherwise give false results.

The script can be run using the command:
cpptraj<cpptraj_rmsd 

### RMSF
This script uses cpptraj to calculate the RMSF of a protein from a trajectory. 
The command 'atomicfluct' is used to specify that the RMSF is to be calculated. 
A mask is created which specifies the residues and atoms to include within the calculation. 
An output file is create using the command 'out' followed by a file name. 
The 'byres' command specifies that the RMSF is calculaded by residue and not by atom.

The script can be run using the command:
cpptraj<cpptraj_rmsf

### Distance 
This script uses cpptraj to calculate the distance of a selection of residues from a different selection of residues. 
The command 'distance' is used to specify that distance measurements are to be calculated. 
The name of the mask is specified e.g. distance 185
A mask is created which specifies the residues and atoms to include within the calculation. 
An output file is create using the command 'out' followed by a file name. 
When the same file name is used for each distance measurment
the results are appended to the same file in separate columns 

The script can be run using the command:
cpptraj<cpptraj_distance

## Electrostatic Potential

prerequiste: 
Pandas, matplotlib, numpy, bio, mdtraj

This script calculates the physiological isoelectric point (pI) and the APBS electrostatics 
### pI calculation:
The difference in solvent accessible surface area is calculated with MDAnalysis using PDB files with and without the labile subdomain. 
The protein is sliced along the Z-axis into 12Å-height regions (run_avg), with a shift of 2.5Å (parameter shift). 
A biopython library function is used to calulate the pI which considers residues with a greater than 20% change in SASA (parameter threshold). 
The script returns the average Z and the pI for each sliced region. 

### Electrostatics calculation:
Adaptive Poisson-Boltzmann Solver(APBS) is used to calculate the electrostaic potential.
PDB2PQR converts the PDB files into PQR files.
APBS returns the electrostatic potential calculated at the nodes of a 3D grid,  each grid cube has a spacing of ~0.5 Å. 
The Scipy library was used to interpolate the electrostatic potential along the pore. 
The pore co-ordinates are defined (prior to script usage) using ChexViS (an webserve that determines transmembrane pores availabile at http://vgl.serc.iisc.ernet.in/chexvis/.
The β-factors in the PDB are replaced with the electrostatic potential values at each co-ordinate.  

To run this script:
python3 pore_mde.py (protein name)_(PDB identifyer)
e.g. python3 pore_mde.py FpvAI_5ODW
 
## Graphs 

Pandas, Numpy and Matplotlib must be installed in local environment. 
The file name is specified to load the data.
The title, x and y axis titles can be adjucted accordingly. 
The line colour, width and lable can be spcified following 'plt.plot'.
The y limitation can be removed or adjusted to fit the data range. 
The figure file name can be adjusted to suit the user.

The distance graph code skip the first row as this row contained headers
and not numerical data

Graph files; RSMD, RMSF and distance are complementary to the data
outputted from the CPPTraj scripts.

Graph files; electrostatic potential +(notscaled), and pI are
complementary to the data outputted from pore_mde.py

