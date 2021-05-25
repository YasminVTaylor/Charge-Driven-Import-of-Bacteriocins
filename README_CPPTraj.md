# Charge-Driven-Import-of-Bacteriocins: CPPTaj scripts 
prerequisite:
To run the scripts cpptraj from AmberTools must be installed on the local machine

## RMSD
This script uses cpptraj to calculate the RMSD of a protein from a trajectory. 
The topology file is loaded using the command 'parm'. 
Each of the trajectory output files are loaded using the command 'trajin'. 
The command 'rms' is used to specify that the RMSD is to be calculated. 
A mask is created which specifies the residues and atoms to include within the calculation.
An output file is create using the command 'out' followed by a file name. 
The 'norotate' command prevents co-ordinates from rotating which may otherwise give false results.

The script can be run using the command:
cpptraj<cpptraj_rmsd 

## RMSF
This script uses cpptraj to calculate the RMSF of a protein from a trajectory. 
The command 'atomicfluct' is used to specify that the RMSF is to be calculated. 
A mask is created which specifies the residues and atoms to include within the calculation. 
An output file is create using the command 'out' followed by a file name. 
The 'byres' command specifies that the RMSF is calculaded by residue and not by atom.

The script can be run using the command:
cpptraj<cpptraj_rmsf

## Distance 
This script uses cpptraj to calculate the distance of a selection of residues from a different selection of residues. 
The command 'distance' is used to specify that distance measurements are to be calculated. 
The name of the mask is specified e.g. distance 185
A mask is created which specifies the residues and atoms to include within the calculation. 
An output file is create using the command 'out' followed by a file name. 
When the same file name is used for each distance measurment
the results are appended to the same file in separate columns 

The script can be run using the command:
cpptraj<cpptraj_distance
