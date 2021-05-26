import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
import Bio
from subprocess import Popen
import os
from scipy.spatial import cKDTree
from itertools import product
import mdtraj as md

# PARAMETERS FOR RUNNING AVERAGE
run_avg = 12    # takes region with delta z = 10
shift = 2.5      # previous region is shifetd by z = 1 
                # e.g.: first region z=10-20, then 11-21, then 12-22, ...
threshold = 0.2 # The difference in ABS should be greater than 10%

def read_z_from_pdb(pdb,res):
    ''' extract average z for residue given as input from pdb file '''
    atoms = []
    with open(pdb) as f:
        for line in f.readlines():
            if line.split()[0] == 'ATOM':
                atoms.append(line.split())

    residues = []; z = []; y = []; x =[] ; amino_acid = []
    for atom in atoms:
        residues.append(int(atom[4]))
        z.append(float(atom[7]))
        y.append(float(atom[6]))
        x.append(float(atom[5]))
        amino_acid.append(atom[3])
    
    zs = []; xs = []; ys = []; aa = []
    for i in range(len(res)):
        z_tmp = []; aa_tmp = []; x_tmp= [] ; y_tmp = []
        for j in range(len(residues)):
            if residues[j] == res[i]:
                z_tmp.append(z[j])
                y_tmp.append(y[j])
                x_tmp.append(x[j])
                if amino_acid[j] not in aa_tmp:
                    aa_tmp.append(amino_acid[j])
        zs.append(np.average(z_tmp))
        xs.append(np.average(x_tmp))
        ys.append(np.average(y_tmp))
        a_acid = aa_tmp[0]
        if a_acid == "HSD":
            a_acid = "HIS"
        aa.append(Bio.SeqUtils.IUPACData.protein_letters_3to1[a_acid.lower().capitalize()])
    return res, aa, zs, ys, xs

def running_average(res,run_avg,shift):
    print("\tCalculating running average")
    start = min(res['z'])
    end = start + run_avg

    mid_values = []; IPs = []; sequences = []; center_of_mass = []; k = 1
    xr = []; yr = []
    while end <= max(res['z']):
        mid_values.append((start+end)/2)
        z_in_region = []; y_in_region = []; x_in_region = []; ab_in_region = []; sequence_in_region = ""
        for i in range(len(res)):
            if res['z'][i] <= end and res['z'][i] >= start:
#                if res['ABS'][i] >= max(res['ABS'])*threshold:
                    x_in_region.append(res['x'][i])
                    y_in_region.append(res['y'][i])
                    z_in_region.append(res['z'][i])
                    ab_in_region.append(res['ABS'][i])
                    sequence_in_region += res['aa'][i]            
        print('Region '+str(k)+': '+sequence_in_region)
        try:
            IPs.append(IP(sequence_in_region).pi())
        except:
            IPs.append(0)
        sequences.append(sequence_in_region)
    
        start += shift
        end += shift
        k += 1
        
        com = [np.average(x_in_region), np.average(y_in_region), (start+end)/2]
        center_of_mass.append(com) 
        
        xr.append(x_in_region)
        yr.append(y_in_region)
    return mid_values, IPs, sequences, center_of_mass, xr, yr

def calculate_electrostatics(protein):
    pdb_file = "pdb/PDB_"+protein+"_NP.pdb"
    print("\t Electrostatics calculation - pdb2pqr")
    os.system("pdb2pqr30 --ff AMBER --apbs-input "+protein+"_APBSINPUT.in "+pdb_file+" "+protein+".pqr")
    try:
        print("\t Electrostatics calculation - APBS")
        os.system("apbs "+protein+"_APBSINPUT.in")
        dx_file = protein+".pqr.dx"
        return dx_file 
    except:
        print('Error: APBSINPUT file not present')

def read_dx_file(dx_file):
    
    print("\t Creating grid from dx file")
    with open(dx_file,'r') as f:
        delta = []; potential = []
        for line in f.readlines():
            if any(x in line for x in ['#','component','attribute','positions regular']):
                continue
            elif "object 1" in line:
                content = line.split()
                nx = int(content[5])
                ny = int(content[6])
                nz = int(content[7])
            elif "origin" in line:
                content = line.split()
                x0 = float(content[1])
                y0 = float(content[2])
                z0 = float(content[3])
            elif "delta" in line:
                content = line.split()
                delta.append([content[1],content[2],content[3]])
            elif "object 2" in line or "object 3" in line:
                continue
            else:
                content = line.split()
                for i in range(len(content)):
                    potential.append(float(content[i]))

    dx = float(delta[0][0])
    dy = float(delta[1][1])
    dz = float(delta[2][2])

    print("*** grid parameters ***")
    print("nx ",nx)
    print("ny ",ny)
    print("nz ",nz)
    print("origin ",[x0,y0,z0])
    print("dx ",dx)
    print("dy ",dy)
    print("dz ",dz)

    x_array = np.linspace(x0, x0+nx*dx, nx)
    y_array = np.linspace(y0, y0+ny*dy, ny)
    z_array = np.linspace(z0, z0+nz*dz, nz)
    yv, xv, zv = np.meshgrid(y_array, x_array, z_array)

    combined_x_y_z_arrays = np.dstack([xv.ravel(), yv.ravel(), zv.ravel()])[0]	
    return combined_x_y_z_arrays, potential

def calculate_sasa(structure):
    sasa_md = md.shrake_rupley(structure,mode='residue',get_mapping=True)
    residues = [int(str(residue)[3:]) for residue in structure.topology.residues]
    sasa = {}
    for i in range(len(residues)):
        sasa[int(residues[i])] = float(sasa_md[0][0][i])
    return sasa

def calculate_difference_in_sasa(pdb1,pdb2,threshold):
    structure1 = md.load_pdb(pdb1)
    structure2 = md.load_pdb(pdb2)
    
    sasa1 = calculate_sasa(structure1)
    sasa2 = calculate_sasa(structure2)
    common_residues = set(list(sasa1.keys())).intersection(list(sasa2.keys()))

    residues = []; diff_sasa = []
    total_surf = 0; part_surf = 0
    for residue in common_residues:
        diff_sasa1_sasa2 = sasa1[residue]-sasa2[residue]
        total_surf += diff_sasa1_sasa2
        if diff_sasa1_sasa2 >= threshold:
            residues.append(residue)
            diff_sasa.append(sasa1[residue]-sasa2[residue])
            part_surf += diff_sasa1_sasa2
    print('Total surface: ',total_surf)
    print('Partial surface: ',part_surf)
    return residues, diff_sasa

if __name__ == '__main__':
    import sys
    protein = sys.argv[1]
    pdb_file = "pdb/PDB_"+protein+"_NP.pdb"
    pdb_file2 = "pdb/PDB_"+protein+"_WP.pdb"
                
    res, dabs = calculate_difference_in_sasa(pdb_file,pdb_file2,threshold)
    
    try:
        sys.argv[2] == 'skip'
        dx_file = protein+'.pqr.dx'
    except:
        dx_file = calculate_electrostatics(protein)
    grid, potential = read_dx_file(dx_file)
    interpolator_tree = cKDTree(grid)
    res, aa, zs, ys, xs = read_z_from_pdb(pdb_file,res)
    
    d = {'res':res, 'aa':aa, 'ABS':dabs, 'x':xs, 'y':ys, 'z':zs}
    res = pd.DataFrame(d)
     
    coords = []
    with open('worms/'+protein+'_pore.pdb','r') as f: #changed here
        for line in f.readlines():
            if 'HETATM' in line:
                content = line.split()
                coords.append([float(content[5]),float(content[6]),float(content[7])])
    f.close()
    
    mid_values, IPs, sequences, center_of_mass, xr, yr = running_average(res,run_avg,shift)
    zs = []
    potentials = []; k = 1
    f = open(protein+"_res.pdb","w+")
    for point in coords:
        zs.append(point[2])        
        dist, indexes = interpolator_tree.query(point, k=[i for i in range(5)]) # add numbers here for smoother interpolation 
                                                              # (e.g.: k=[1,2,3,4] for interpolation with 4 closest neighbours)
        pnts = []; pots = []
        for index in indexes:
            pnt = [grid[index][0], grid[index][1], grid[index][2]]
            pnts.append(pnt)
            pots.append(potential[index])
        pot = np.average(pots)#, weights=1/dist)
        potentials.append(pot)
        f.write('%s%s %s %s %s%s    %s%s%s%s%s%s\n' % ("ATOM".ljust(6),
                str(k).rjust(5),
                "CA".center(4),
                "ALA".ljust(3),
                "A".rjust(1),
                str(k).rjust(4),
                str('%8.3f' % (float(point[0]))).rjust(8),
                str('%8.3f' % (float(point[1]))).rjust(8),
                str('%8.3f' % (float(point[2]))).rjust(8),
                str('%6.2f' % (float(1.00))).rjust(6),
                str('%6.2f' % (float(round(-pot,2)))).ljust(6),
                "C".rjust(12)))
        k += 1
    f.write("TER\nEND")
    f.close()
    np.savetxt(protein+"_res_ip.txt", np.array([mid_values,IPs]).T)
    np.savetxt(protein+"_res_pot.txt", np.array([zs,potentials]).T)

    fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)
 #   plt.title(protein+" pI and pore Electrostatic Potential")    
    ax1.plot(mid_values,IPs,':',color='black',label='isoelec point')
    ax2.plot(zs,potentials,'-',color='black',label='interp potential')
    ax2.set_xlabel('z',fontsize=15)
    ax1.set_ylabel('Isoelectric Point',fontsize=10)
    ax2.set_ylabel('Interpolated Potential',fontsize=10)
    plt.savefig(protein+"_pi_ep.png")
    plt.show()       
    
    
