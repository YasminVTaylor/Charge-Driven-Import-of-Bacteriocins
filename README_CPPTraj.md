## Charge-Driven-Import-of-Bacteriocins: CPPTaj scripts 

# RMSD
This script uses cpptraj to calculate the RMSD of a protein from a trajectory. 
The topology file is loaded using the command 'parm'. 
Each of the trajectory output files are loaded using the command 'trajin'. 
The command 'rms' is used to specify that the RMSD is to be calculated. 
A mask is created which specifies the residues and atoms to include within the calculation.
An output file is create using the command 'out' followed by a file name. 
The 'norotate' command prevents co-ordinates from rotating which may otherwise give false results.

To run the script cpptraj from AmberTools must be installed on the local machine
The script can be run using the command:
cpptraj<cpptraj_rmsd_script.txt

# RMSF
This script uses cpptraj to calculate the RMSF of a protein from a trajectory. 
The topology file is loaded using the command 'parm'. 
Each of the trajectory output files are loaded using the command 'trajin'. 
The command 'atomicfluct' is used to specify that the RMSF is to be calculated. 
A mask is created which specifies the residues and atoms to include within the calculation. 
An output file is create using the command 'out' followed by a file name. 
The 'byres' command specifies that the RMSF is calculaded by residue and not by atom.

To run the script cpptraj from AmberTools must be installed on the local machine
The script can be run using the command:
cpptraj<cpptraj_rmsf_script.txt

# Distance 
