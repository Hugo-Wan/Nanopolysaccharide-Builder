import MDAnalysis as mda
from MDAnalysis.coordinates import PDB
from MDAnalysis.core.universe import Merge
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.coordinates.PDB import PDBWriter
from MDAnalysis.core.topologyattrs import Elements
from vmd import atomsel, molecule
import math
import copy
import random
import numpy as np
import warnings
import os
import glob
import shutil
import json
import sys
import re
warnings.filterwarnings("ignore", category=UserWarning)


#a_trans = 7.784
#b_trans = 8.201
#c_trans = 10.38
#gamma_angle=96.5


try:
    a_trans = float(sys.argv[1])
    b_trans = float(sys.argv[2])
    c_trans = float(sys.argv[3])
    gamma_angle = float(sys.argv[4])
    if a_trans <= 0 or b_trans <= 0 or c_trans <= 0 or gamma_angle<=0:  
        sys.stderr.write("Error: Please provide positive value for all crystallographic parameters.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length or angle value provided. Please enter valid numeric values.\n")
    sys.exit(1)
except IndexError:
    sys.stderr.write("Error: Missing one or more input arguments. Please provide four required crystallographic values.\n")
    sys.exit(1)


try:
    c_iterations = int(sys.argv[5]) 
    if c_iterations < 0 :
        sys.stderr.write("Error: Please provide a positive integer for c repetition numbers.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length value provide. Please enter a valid numeric value.\n")
    sys.exit(1)



try:
        ##Sulfate degree ref:/10.1016/j.carbpol.2019.115292
        Sulfate_target = float(sys.argv[6])  #####
        if not (0 <= Sulfate_target <= 1.5):
            sys.stderr.write("Error: Sulfate_target should be between 0 and 1.5.\n")
            sys.exit(1)


except IndexError:
    sys.stderr.write("Error: Missing one or more input arguments. Please provide all five effective values.\n")
    sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid input. Please ensure all values are numeric and within the valid ranges.\n")
    sys.exit(1)            # pH condition

##structure input
with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
neutron_pdb_1 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'charmm36', 'chain-1.pdb')
user_defined_pdb_1 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'charmm36', 'chain-1_ud.pdb')

if os.path.exists(user_defined_pdb_1):
    unit_chain_input_file_1 = user_defined_pdb_1
else:
    unit_chain_input_file_1 = neutron_pdb_1

neutron_pdb_2 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'charmm36', 'chain-2.pdb')
user_defined_pdb_2 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'charmm36', 'chain-2_ud.pdb')

if os.path.exists(user_defined_pdb_2):
    unit_chain_input_file_2 = user_defined_pdb_2
else:
    unit_chain_input_file_2 = neutron_pdb_2

chain_1_u = mda.Universe(unit_chain_input_file_1)
chain_2_u = mda.Universe(unit_chain_input_file_2)



m_Glc=163.09316                # Glc unit mass
m_difference=101.1             # difference between Glc unit with Glc-SO3- Na(+)


#chain_1 assembly    
#------------------------------------------------------------------------------
chain_1_strand = []
for chain_1_i in range(1, c_iterations + 1):
    chain_1_u.atoms.positions += [0, 0, c_trans]
    resid_1 = (c_iterations - chain_1_i + 1) * 2 - 1
    resid_2 = (c_iterations - chain_1_i + 1) * 2
    for j, atom in enumerate(chain_1_u.atoms):
        if j < 21:
            atom.residue.resid = resid_1
        else:
            atom.residue.resid = resid_2
    chain_1_output = f"chain_1-{chain_1_i}.pdb"
    chain_1_strand.append(chain_1_output )
    with mda.Writer(chain_1_output , n_atoms=chain_1_u.atoms.n_atoms) as W:
        W.write(chain_1_u.atoms)
with open("chain_1_temp.pdb", "w") as chain_1_structure:
    for chain_1_output in reversed (chain_1_strand):
        with open(chain_1_output, "r") as chain_1_pdb_file:
            for chain_1_line in chain_1_pdb_file:
                if chain_1_line.startswith("ATOM"):
                    chain_1_structure.write(chain_1_line)
            os.remove(chain_1_output) 


 ##terminal side            
chain_1_finite_build_chain = "chain_1_temp.pdb"
chain_1_finite_u = mda.Universe(chain_1_finite_build_chain)
chain_1_first_atom = chain_1_finite_u.atoms[0]
chain_1_o4_atoms = chain_1_finite_u.select_atoms("name O4")
if not chain_1_o4_atoms:
    raise ValueError("No atoms named 'O4' found in the structure.")
chain_1_last_o4 = chain_1_o4_atoms[-1]

chain_1_new_positions = {
    'O1': chain_1_first_atom.position + [0.538, -0.449, 1.247],
    'HO1': chain_1_first_atom.position + [0.129, -0.365, 1.932],
    'HO4': chain_1_last_o4.position + [0.377, -0.192, -0.892]
}

chain_1_new_atoms = {}
for name, pos in chain_1_new_positions.items():
    chain_1_base_atom = chain_1_first_atom if name in ['O1', 'HO1'] else chain_1_last_o4
    chain_1_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
    chain_1_new_uni.add_TopologyAttr('name', [name])
    chain_1_new_uni.add_TopologyAttr('type', [chain_1_base_atom.type])
    chain_1_new_uni.add_TopologyAttr('resname', [chain_1_base_atom.resname])
    chain_1_new_uni.add_TopologyAttr('resid', [chain_1_base_atom.resid])
    chain_1_new_uni.add_TopologyAttr('segid', ['0']) 
    chain_1_new_uni.atoms.positions = [pos]
    chain_1_new_atoms[name] = chain_1_new_uni

chain_1_combined = Merge(chain_1_finite_u.atoms[:chain_1_last_o4.index + 1], chain_1_new_atoms['HO4'].atoms, chain_1_finite_u.atoms[chain_1_last_o4.index + 1:])
chain_1_combined = Merge(chain_1_combined.atoms[:2], chain_1_new_atoms['O1'].atoms, chain_1_combined.atoms[2:])
chain_1_combined = Merge(chain_1_combined.atoms[:3], chain_1_new_atoms['HO1'].atoms, chain_1_combined.atoms[3:])



for atom in chain_1_combined.atoms:
    atom.segment.segid = '0'  
    elements = [atom.name[0] for atom in chain_1_combined.atoms]
    chain_1_combined.add_TopologyAttr(Elements(elements))
chain_1_combined.atoms.write("chain_1_temp.2.pdb")
#------------------------------------------------------------------------------



#chain_1 assembly    
#------------------------------------------------------------------------------
chain_2_strand = []
for chain_2_i in range(1, c_iterations + 1):
    chain_2_u.atoms.positions += [0, 0, c_trans]
    resid_1 = (c_iterations - chain_2_i + 1) * 2 - 1
    resid_2 = (c_iterations - chain_2_i + 1) * 2
    for j, atom in enumerate(chain_2_u.atoms):
        if j < 21:
            atom.residue.resid = resid_1
        else:
            atom.residue.resid = resid_2
    chain_2_output = f"chain_2-{chain_2_i}.pdb"
    chain_2_strand.append(chain_2_output )
    with mda.Writer(chain_2_output , n_atoms=chain_2_u.atoms.n_atoms) as W:
        W.write(chain_2_u.atoms)
with open("chain_2_temp.pdb", "w") as chain_2_structure:
    for chain_2_output in reversed (chain_2_strand):
        with open(chain_2_output, "r") as chain_2_pdb_file:
            for chain_2_line in chain_2_pdb_file:
                if chain_2_line.startswith("ATOM"):
                    chain_2_structure.write(chain_2_line)
            os.remove(chain_2_output) 


 ##terminal side            
chain_2_finite_build_chain = "chain_2_temp.pdb"
chain_2_finite_u = mda.Universe(chain_2_finite_build_chain)
chain_2_first_atom = chain_2_finite_u.atoms[0]
chain_2_o4_atoms = chain_2_finite_u.select_atoms("name O4")
if not chain_2_o4_atoms:
    raise ValueError("No atoms named 'O4' found in the structure.")
chain_2_last_o4 = chain_2_o4_atoms[-1]

chain_2_new_positions = {
    'O1': chain_2_first_atom.position + [0.538, -0.449, 1.247],
    'HO1': chain_2_first_atom.position + [0.129, -0.365, 1.932],
    'HO4': chain_2_last_o4.position + [0.377, -0.192, -0.892]
}

chain_2_new_atoms = {}
for name, pos in chain_2_new_positions.items():
    chain_2_base_atom = chain_2_first_atom if name in ['O1', 'HO1'] else chain_2_last_o4
    chain_2_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
    chain_2_new_uni.add_TopologyAttr('name', [name])
    chain_2_new_uni.add_TopologyAttr('type', [chain_2_base_atom.type])
    chain_2_new_uni.add_TopologyAttr('resname', [chain_2_base_atom.resname])
    chain_2_new_uni.add_TopologyAttr('resid', [chain_2_base_atom.resid])
    chain_2_new_uni.add_TopologyAttr('segid', ['0']) 
    chain_2_new_uni.atoms.positions = [pos]
    chain_2_new_atoms[name] = chain_2_new_uni

chain_2_combined = Merge(chain_2_finite_u.atoms[:chain_2_last_o4.index + 1], chain_2_new_atoms['HO4'].atoms, chain_2_finite_u.atoms[chain_2_last_o4.index + 1:])
chain_2_combined = Merge(chain_2_combined.atoms[:2], chain_2_new_atoms['O1'].atoms, chain_2_combined.atoms[2:])
chain_2_combined = Merge(chain_2_combined.atoms[:3], chain_2_new_atoms['HO1'].atoms, chain_2_combined.atoms[3:])



for atom in chain_2_combined.atoms:
    atom.segment.segid = '0'  
    elements = [atom.name[0] for atom in chain_2_combined.atoms]
    chain_2_combined.add_TopologyAttr(Elements(elements))
chain_2_combined.atoms.write("chain_2_temp.2.pdb")
#------------------------------------------------------------------------------
#
##------------------------------------------------------------------------------
###translate parameter
tilt_angle_rad= (gamma_angle - 90) * math.pi/180
angle_rad = (gamma_angle - 90) * math.pi / 180
cosA = math.cos(angle_rad)
sinA = math.sin(angle_rad)
cosB = math.cos(angle_rad - tilt_angle_rad)
sinB = math.sin(angle_rad - tilt_angle_rad)
a_par_transverse_move= a_trans * sinB 
a_par_vertical_move= -1 * a_trans * cosB 
b_par_transverse_move= b_trans * cosA 
b_par_vertical_move= -1 * b_trans * sinA  



#chain_1_1st_layer assembly  
chain_1_1st_layer_input_file = "chain_1_temp.2.pdb"
chain_1_1st_layer_u = mda.Universe(chain_1_1st_layer_input_file)
chain_1_1st_layer = []  
for i in range(1, 5):
    chain_1_1st_layer_u.atoms.positions += [ b_par_vertical_move, b_par_transverse_move, 0]

    segid_increment_value = 1 
    for atom in chain_1_1st_layer_u.atoms:
        numeric_part = atom.segid.replace('','')
        new_numeric_part = int(numeric_part) + segid_increment_value
    new_segid = f"{new_numeric_part}"
    atom.segment.segid = new_segid

    chain_1_1st_layer_output = f"chain_1_1st_layer_{i}.pdb"
    chain_1_1st_layer.append(chain_1_1st_layer_output)
    with mda.Writer(chain_1_1st_layer_output, n_atoms=chain_1_1st_layer_u.atoms.n_atoms) as W:
        W.write(chain_1_1st_layer_u.atoms)
with open("chain_1_1st_layer_temp.pdb", "w") as layer_file:
    for chain_1_1st_layer_output in chain_1_1st_layer:
        with open(chain_1_1st_layer_output, "r") as chain_1_1st_layer_pdb_file:
            for chain_1_1st_layer_line in chain_1_1st_layer_pdb_file:
                if chain_1_1st_layer_line.startswith("ATOM"):
                    layer_file.write(chain_1_1st_layer_line)
            os.remove(chain_1_1st_layer_output) 


#chain_1_2nd_layer assembly  
chain_1_2nd_layer_input_file = "chain_1_temp.2.pdb"
chain_1_2nd_layer_u = mda.Universe(chain_1_1st_layer_input_file)
chain_1_2nd_layer = [] 

for j in range(1, 3):
    k=j+1
    chain_1_2nd_layer_u = mda.Universe(chain_1_2nd_layer_input_file)
    translation_vector = [a_par_vertical_move + b_par_vertical_move * k,
                          a_par_transverse_move + b_par_transverse_move * k, 0]
    chain_1_2nd_layer_u.atoms.positions += translation_vector
    for atom in chain_1_2nd_layer_u.atoms:
        new_numeric_part = 4 + j  
        atom.segment.segid = str(new_numeric_part)  
    chain_1_2nd_layer_output = f"chain_1_2nd_layer_{j}.pdb"
    chain_1_2nd_layer.append(chain_1_2nd_layer_output)
    with mda.Writer(chain_1_2nd_layer_output, n_atoms=chain_1_2nd_layer_u.atoms.n_atoms) as W:
        W.write(chain_1_2nd_layer_u.atoms)
with open("chain_1_2nd_layer_temp.pdb", "w") as layer_file:
    for chain_1_2nd_layer_output in chain_1_2nd_layer:
        with open(chain_1_2nd_layer_output, "r") as chain_1_2nd_layer_pdb_file:
            for line in chain_1_2nd_layer_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
            os.remove(chain_1_2nd_layer_output)


#chain_1_3rd_layer assembly  
chain_1_3rd_layer_input_file = "chain_1_temp.2.pdb"
chain_1_3rd_layer_u = mda.Universe(chain_1_1st_layer_input_file)
chain_1_3rd_layer = [] 

for l in range(1, 4):  
    m = l - 1  
    chain_1_3rd_layer_u = mda.Universe(chain_1_3rd_layer_input_file)
    translation_vector = [-a_par_vertical_move + b_par_vertical_move * l,
                          -a_par_transverse_move + b_par_transverse_move * l, 0]

    chain_1_3rd_layer_u.atoms.positions += translation_vector

    for atom in chain_1_3rd_layer_u.atoms:
        new_numeric_part = 6 + l  
        atom.segment.segid = str(new_numeric_part)  
    chain_1_3rd_layer_output = f"chain_1_3rd_layer_{l}.pdb"
    chain_1_3rd_layer.append(chain_1_3rd_layer_output)
    with mda.Writer(chain_1_3rd_layer_output, n_atoms=chain_1_3rd_layer_u.atoms.n_atoms) as W:
        W.write(chain_1_3rd_layer_u.atoms)
with open("chain_1_3rd_layer_temp.pdb", "w") as layer_file:
    for chain_1_3rd_layer_output in chain_1_3rd_layer:
        with open(chain_1_3rd_layer_output, "r") as chain_1_3rd_layer_pdb_file:
            for line in chain_1_3rd_layer_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
            os.remove(chain_1_3rd_layer_output)




#chain_2_1st_layer assembly  
chain_2_1st_layer_input_file = "chain_2_temp.2.pdb"
chain_2_1st_layer_u = mda.Universe(chain_2_1st_layer_input_file)
chain_2_1st_layer = []  

for i in range(1, 5):  
    chain_2_1st_layer_u = mda.Universe(chain_2_1st_layer_input_file)
    translation_vector = [b_par_vertical_move * i,
                          b_par_transverse_move * i, 0]

    chain_2_1st_layer_u.atoms.positions += translation_vector

    for atom in chain_2_1st_layer_u.atoms:
        new_numeric_part = 9 + i  
        atom.segment.segid = str(new_numeric_part)  
    chain_2_1st_layer_output = f"chain_2_1st_layer_{i}.pdb"
    chain_2_1st_layer.append(chain_2_1st_layer_output)
    with mda.Writer(chain_2_1st_layer_output, n_atoms=chain_2_1st_layer_u.atoms.n_atoms) as W:
        W.write(chain_2_1st_layer_u.atoms)
with open("chain_2_1st_layer_temp.pdb", "w") as layer_file:
    for chain_2_1st_layer_output in chain_2_1st_layer:
        with open(chain_2_1st_layer_output, "r") as chain_2_1st_layer_pdb_file:
            for line in chain_2_1st_layer_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
            os.remove(chain_2_1st_layer_output)




#chain_2_2nd_layer assembly  
chain_2_2nd_layer_input_file = "chain_2_temp.2.pdb"
chain_2_2nd_layer_u = mda.Universe(chain_2_1st_layer_input_file)
chain_2_2nd_layer = [] 

for j in range(1, 4):
    k=j+1
    chain_2_2nd_layer_u = mda.Universe(chain_2_2nd_layer_input_file)
    translation_vector = [a_par_vertical_move + b_par_vertical_move * k,
                          a_par_transverse_move + b_par_transverse_move * k, 0]
    chain_2_2nd_layer_u.atoms.positions += translation_vector
    for atom in chain_2_2nd_layer_u.atoms:
        new_numeric_part = 13 + j  
        atom.segment.segid = str(new_numeric_part)  
    chain_2_2nd_layer_output = f"chain_2_2nd_layer_{j}.pdb"
    chain_2_2nd_layer.append(chain_2_2nd_layer_output)
    with mda.Writer(chain_2_2nd_layer_output, n_atoms=chain_2_2nd_layer_u.atoms.n_atoms) as W:
        W.write(chain_2_2nd_layer_u.atoms)
with open("chain_2_2nd_layer_temp.pdb", "w") as layer_file:
    for chain_2_2nd_layer_output in chain_2_2nd_layer:
        with open(chain_2_2nd_layer_output, "r") as chain_2_2nd_layer_pdb_file:
            for line in chain_2_2nd_layer_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
            os.remove(chain_2_2nd_layer_output)


#chain_2_3rd_layer assembly  
chain_2_3rd_layer_input_file = "chain_2_temp.2.pdb"
chain_2_3rd_layer_u = mda.Universe(chain_2_1st_layer_input_file)
chain_2_3rd_layer = [] 

for l in range(1, 3):  
    m = l + 1  
    chain_2_3rd_layer_u = mda.Universe(chain_2_3rd_layer_input_file)
    translation_vector = [-a_par_vertical_move + b_par_vertical_move * m,
                          -a_par_transverse_move + b_par_transverse_move * m, 0]

    chain_2_3rd_layer_u.atoms.positions += translation_vector

    for atom in chain_2_3rd_layer_u.atoms:
        new_numeric_part = 16 + l  
        atom.segment.segid = str(new_numeric_part)  
    chain_2_3rd_layer_output = f"chain_2_3rd_layer_{l}.pdb"
    chain_2_3rd_layer.append(chain_2_3rd_layer_output)
    with mda.Writer(chain_2_3rd_layer_output, n_atoms=chain_2_3rd_layer_u.atoms.n_atoms) as W:
        W.write(chain_2_3rd_layer_u.atoms)
with open("chain_2_3rd_layer_temp.pdb", "w") as layer_file:
    for chain_2_3rd_layer_output in chain_2_3rd_layer:
        with open(chain_2_3rd_layer_output, "r") as chain_2_3rd_layer_pdb_file:
            for line in chain_2_3rd_layer_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
            os.remove(chain_2_3rd_layer_output)


###assemble
def assemble_pdbs(output_file, input_files):
    with open(output_file, 'w') as outfile:
        for filename in input_files:
            with open(filename, 'r') as infile:
                for line in infile:
                    if line.startswith("ATOM"):
                        outfile.write(line)
            outfile.write('TER\n')  

chain_temp_pdb_files = [
    "chain_1_1st_layer_temp.pdb", "chain_1_2nd_layer_temp.pdb", "chain_1_3rd_layer_temp.pdb",
    "chain_2_1st_layer_temp.pdb", "chain_2_2nd_layer_temp.pdb", "chain_2_3rd_layer_temp.pdb"
]

unit_temp_pdb = "unit_temp.pdb"
assemble_pdbs(unit_temp_pdb, chain_temp_pdb_files)




##crystallographic planes defined
planes_110_left = [0,4,13]
planes_110_right = [8,12,17]


planes_1_10_left = [6, 9, 16]
planes_1_10_right = [3,5,15]

planes_110_count=len(planes_110_left)+len(planes_110_right)
planes_1_10_count=len(planes_1_10_left)+len(planes_1_10_right)

chain_u = mda.Universe("unit_temp.pdb")
residue_ids = [residue.resid for residue in chain_u.residues]
maximum_resid_num=residue_ids[-1]
total_residues=len(residue_ids)


ds_reverse= (1/(m_Glc*Sulfate_target*(10**-3))-(m_difference/m_Glc) )  ##(1/ds)
ds=1/ds_reverse
n_so3=round(total_residues*ds)
#print('degree of substitution', ds)

if n_so3 < 0 :
    n_so3=total_residues
#print('carboxylation_number',n_so3)

def adjust_for_sulfate_limit(initial_carboxylation_num):
    initial_sulfate_num = n_so3
    sulfate_number_limit = c_iterations * planes_110_count +  c_iterations * planes_1_10_count
    new_sulfate_number = initial_carboxylation_num
    if initial_sulfate_num > sulfate_number_limit:
        new_sulfate_number = sulfate_number_limit
    return new_sulfate_number

sulfate_num = adjust_for_sulfate_limit(n_so3)
#print(sulfate_num)


if sulfate_num > 0:
    left_110_count=0
    right_110_count=0
    left_1_10_count=0
    right_1_10_count=0
    for i in range(1, sulfate_num + 1):
        # Determine the target based on i % 4
        if i % 4 == 0:
            left_110_count += 1
        elif i % 4 == 1:
            right_110_count += 1
        elif i % 4 == 2:
            left_1_10_count += 1
        elif i % 4 == 3:
            right_1_10_count += 1
    
    max_110_count_left = c_iterations * len(planes_110_left)
    max_110_count_right = c_iterations * len(planes_110_right)

    max_1_10_count_left = c_iterations * len(planes_1_10_left)
    max_1_10_count_right = c_iterations * len(planes_1_10_right)

    left_110_count_update=left_110_count
    right_110_count_update=right_110_count
    left_1_10_count_update=left_1_10_count
    right_1_10_count_update=right_1_10_count

    if left_110_count > max_110_count_left:
        left_110_count_update = max_110_count_left
    
    if right_110_count > max_110_count_right:
        right_110_count_update = max_110_count_right
    
    if left_1_10_count > max_1_10_count_left:
        left_1_10_count_update = max_1_10_count_left
    
    if right_1_10_count > max_1_10_count_right:
        right_1_10_count_update = max_1_10_count_right    



    new_sulfate_number_update=left_110_count_update + right_110_count_update + left_1_10_count_update + right_1_10_count_update 
    #print(new_sulfate_number_update)
    #print(total_residues)
    actual_ds=(new_sulfate_number_update)/(total_residues)
    actual_sulfate_content=Sulfate_target*(actual_ds/ds)
    #print(actual_ds)
    #print(actual_sulfate_content)
    actual_sulfate_rounded=round(actual_sulfate_content,4)
    #print('surface charge density', actual_sulfate_rounded, 'mmol/g')
#
    even_numbers = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 0]
    odd_numbers = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 1]
    all_combinations_110_left = [(value1, value2) for value1 in planes_110_left  for value2 in odd_numbers]
    all_combinations_110_right = [(value1, value2) for value1 in planes_110_right  for value2 in even_numbers ]
    random.shuffle(all_combinations_110_left)
    random.shuffle(all_combinations_110_right)
    sel_110_left_array = []
    sel_110_right_array = []



    while len(sel_110_left_array) < left_110_count_update:
        if not all_combinations_110_left:  
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_110_left = random.choice(all_combinations_110_left)
        sel_110_left_array.append(selected_combination_110_left)
        all_combinations_110_left.remove(selected_combination_110_left)
    
    while len(sel_110_right_array) < right_110_count_update:
        if not all_combinations_110_right:  # Check if the list is empty to avoid IndexError
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_110_right = random.choice(all_combinations_110_right)
        sel_110_right_array.append(selected_combination_110_right)
        all_combinations_110_right.remove(selected_combination_110_right)



    all_combinations_1_10_left = [(value1, value2) for value1 in planes_1_10_left  for value2 in odd_numbers]
    all_combinations_1_10_right = [(value1, value2) for value1 in planes_1_10_right  for value2 in even_numbers ]
    random.shuffle(all_combinations_1_10_left)
    random.shuffle(all_combinations_1_10_right)
    sel_1_10_left_array = []
    sel_1_10_right_array = []


    while len(sel_1_10_left_array) < left_1_10_count_update:
        if not all_combinations_1_10_left:  
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_1_10_left = random.choice(all_combinations_1_10_left)
        sel_1_10_left_array.append(selected_combination_1_10_left)
        all_combinations_1_10_left.remove(selected_combination_1_10_left)
    
    while len(sel_1_10_right_array) < right_1_10_count_update:
        if not all_combinations_1_10_right:  # Check if the list is empty to avoid IndexError
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_1_10_right = random.choice(all_combinations_1_10_right)
        sel_1_10_right_array.append(selected_combination_1_10_right)
        all_combinations_1_10_right.remove(selected_combination_1_10_right)



    def get_fragment_indices(mol_id, sel):
        min_max_chain_indices = []
        for fragment_number, resid in sel:
            selection_query = f"fragment {fragment_number} and resid {resid}"
            selection = atomsel(selection_query, molid=mol_id)
            indices = selection.index
            if indices:
                min_index = min(indices)
                max_index = max(indices)
                min_max_chain_indices.append((min_index, max_index))
            else:
                min_max_chain_indices.append((None, None))
        return min_max_chain_indices
    
    origin_glc_chains = molecule.load("pdb", "unit_temp.pdb")
    left_glc_110_candidates_indices   = get_fragment_indices(origin_glc_chains, sel_110_left_array)   
    right_glc_110_candidates_indices  = get_fragment_indices(origin_glc_chains, sel_110_right_array)   
    left_glc_1_10_candidates_indices  = get_fragment_indices(origin_glc_chains, sel_1_10_left_array)   
    right_glc_1_10_candidates_indices = get_fragment_indices(origin_glc_chains, sel_1_10_right_array)  
    molecule.delete(origin_glc_chains) 

    def remove_atoms_based_on_indices(universe, index_ranges):
        all_removed_indices = []
    
        for index_range in index_ranges:
            for start_idx, end_idx in index_range:
                atom_group = universe.atoms[start_idx:end_idx + 1]
                for residue in atom_group.residues:
                    if residue.resid == 1:
                        indices_to_remove = residue.atoms[19:23].indices
                    elif residue.resid == maximum_resid_num:
                        indices_to_remove = residue.atoms[18:22].indices
                    else:
                        indices_to_remove = residue.atoms[17:21].indices
                    all_removed_indices.extend(indices_to_remove)
        unique_indices = sorted(set(all_removed_indices))
        indices_str = ' '.join(map(str, unique_indices))
        universe = universe.select_atoms(f"not index {indices_str}")
    
        return universe
    
    all_indices = [left_glc_110_candidates_indices, right_glc_110_candidates_indices, 
                   left_glc_1_10_candidates_indices, right_glc_1_10_candidates_indices]
    ###chain_u is the unmodified structure, with name of "alpha-chitin.temp.pdb"
    primary_hydroxyl_remove_universe = remove_atoms_based_on_indices(chain_u, all_indices)
    primary_hydroxyl_remove_universe.atoms.write("primary_hydroxyl_remove.temp.pdb")

    def get_fragment(mol_id, indices_i, indices_j):
        selection_query = f"index {indices_i} to {indices_j}"
        selection = atomsel(selection_query, molid=mol_id)
        fragment_num_indice = selection.fragment
        unique_fragments = set(fragment_num_indice)
        if len(unique_fragments) == 1:
                return unique_fragments.pop()   
        else:
            raise ValueError("Multiple unique fragments found in selection")
   
    primary_hydroxyl_remove_chains = molecule.load("pdb", "primary_hydroxyl_remove.temp.pdb")
    left_so3_candidates_indices   = get_fragment_indices(primary_hydroxyl_remove_chains, sel_110_left_array)   
    right_so3_candidates_indices  = get_fragment_indices(primary_hydroxyl_remove_chains, sel_110_right_array)   
    upper_so3_candidates_indices  = get_fragment_indices(primary_hydroxyl_remove_chains, sel_1_10_left_array)   
    bottom_so3_candidates_indices = get_fragment_indices(primary_hydroxyl_remove_chains, sel_1_10_right_array)  
    primary_hydroxyl_remove_chains_u = mda.Universe("primary_hydroxyl_remove.temp.pdb")
    so3_atoms = list(primary_hydroxyl_remove_chains_u.atoms)#
#

    def add_so3_based_on_indices(universe, indices_lists):
        all_new_atoms = []  # This will store all newly created atoms across all modifications
        for index_ranges in indices_lists:
            for start_idx, end_idx in index_ranges:
                residues = universe.select_atoms(f"index {start_idx}:{end_idx}").residues
                fragment_number = get_fragment(primary_hydroxyl_remove_chains, start_idx, end_idx)
                #print("staring from to the end:",start_idx, end_idx)
                #print("fragment number is :", fragment_number )
                #print("residue is : ", residues)
                for residue in residues:
                    #print("residue is : ", residue)
                    residue.resname = 'BGLS'
                    if residue.resid == 1 :
                       base_atom_index = 18   
                    elif residue.resid == maximum_resid_num :
                       base_atom_index = 17  
                    else :
                       base_atom_index = 16   
                    base_atom = residue.atoms[base_atom_index]
                    so3_new_positions = {} 
                    # Define new positions based on residue ID's odd/even status
                    if residue.resid % 2 == 1  and fragment_number in planes_110_left:  # Odd residue ID
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([ 0.889,    0.251,  -0.289]),
                            'H62':  base_atom.position + np.array([ 0.081,   -0.511,   0.821]),
                            'O6':   base_atom.position + np.array([ 0.478,   -1.235,  -0.500]),
                            'S6':   base_atom.position + np.array([ 1.141,   -2.655,  -1.078]),
                            'OS62': base_atom.position + np.array([ 1.791,   -3.523,   0.193]),
                            'OS63': base_atom.position + np.array([ 2.343,   -2.299,  -2.181]),
                            'OS64': base_atom.position + np.array([-0.047,   -3.562,  -1.823])
                        }
                    elif residue.resid % 2 == 0  and fragment_number in planes_110_right:  # Even residue ID
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([-0.381,  -0.045,  -0.888]),
                            'H62':  base_atom.position + np.array([-0.729,   0.050,   0.640]),                            
                            'O6':   base_atom.position + np.array([ 0.320,   1.379,  -0.012]),
                            'S6':   base_atom.position + np.array([ 0.697,   3.006,  -0.025]),
                            'OS62': base_atom.position + np.array([ 1.161,   3.482,   1.508]),
                            'OS63': base_atom.position + np.array([ 1.952,   3.278,  -1.094]),
                            'OS64': base_atom.position + np.array([-0.640,   3.884,  -0.501])
                        }
    
                    elif residue.resid % 2 == 1  and fragment_number in planes_1_10_left:  # Odd residue ID
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([ 0.889,    0.251,  -0.289]),
                            'H62':  base_atom.position + np.array([ 0.081,   -0.511,   0.821]),
                            'O6':   base_atom.position + np.array([ 0.478,   -1.235,  -0.500]),
                            'S6':   base_atom.position + np.array([ 1.141,   -2.655,  -1.078]),
                            'OS62': base_atom.position + np.array([ 1.791,   -3.523,   0.193]),
                            'OS63': base_atom.position + np.array([ 2.343,   -2.299,  -2.181]),
                            'OS64': base_atom.position + np.array([-0.047,   -3.562,  -1.823])
                        }
                    elif residue.resid % 2 == 0  and fragment_number in planes_1_10_right:  # Even residue ID
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([-0.381,  -0.045,  -0.888]),
                            'H62':  base_atom.position + np.array([-0.729,   0.050,   0.640]),                            
                            'O6':   base_atom.position + np.array([ 0.320,   1.379,  -0.012]),
                            'S6':   base_atom.position + np.array([ 0.697,   3.006,  -0.025]),
                            'OS62': base_atom.position + np.array([ 1.161,   3.482,   1.508]),
                            'OS63': base_atom.position + np.array([ 1.952,   3.278,  -1.094]),
                            'OS64': base_atom.position + np.array([-0.640,   3.884,  -0.501])
                        }
    

                    insert_pos = so3_atoms.index(base_atom) + 1
                    for name, pos in so3_new_positions.items():
                        so3_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
                        so3_new_uni.add_TopologyAttr('name', [name])
                        so3_new_uni.add_TopologyAttr('type', [base_atom.type])
                        so3_new_uni.add_TopologyAttr('resname', [base_atom.resname])
                        so3_new_uni.add_TopologyAttr('resid', [residue.resid])
                        so3_new_uni.add_TopologyAttr('segid', [base_atom.segment.segid])
                        so3_new_uni.add_TopologyAttr('chainIDs', ['X'])
                        so3_new_uni.atoms.positions = [pos]
                        so3_new_atom = so3_new_uni.atoms[0]
    
                        so3_atoms.insert(insert_pos, so3_new_atom)
                        insert_pos += 1  # Update position for the next atom
        # Create a new universe with all atoms including the added ones
        new_universe = mda.Merge(*[mda.AtomGroup([atom]) for atom in so3_atoms])
        return new_universe
    
    all_indices = [left_so3_candidates_indices, right_so3_candidates_indices, upper_so3_candidates_indices, bottom_so3_candidates_indices]
    so3_negative_u = add_so3_based_on_indices(primary_hydroxyl_remove_chains_u, all_indices)
    so3_negative_u.atoms.write("so3_negative-modified.temp.pdb")
    molecule.delete(primary_hydroxyl_remove_chains)#
#
#
#    
    u_final = mda.Universe("so3_negative-modified.temp.pdb")
    box_x = 50
    box_y = 50
    box_z = c_iterations * c_trans + 20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    u_final.atoms.translate(translation_vector)
    with mda.Writer(f"charmm36-cellulose-Ibeta-18fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
        
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    
    print(f"Half-Ester sulfate content: {actual_sulfate_content:.2f}, Degree of sulfate: {actual_ds:.4f}")


elif sulfate_num == 0:
    u_final = mda.Universe("unit_temp.pdb")
    box_x = 50
    box_y = 50
    box_z = c_iterations * c_trans + 20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    with mda.Writer(f"charmm36-cellulose-Ibeta-18fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
    

    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    actual_ds=0.0    
    actual_sulfate_content=0.0
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"Half-Ester sulfate content: {actual_sulfate_content:.2f}, Degree of sulfate: {actual_ds:.4f}")

