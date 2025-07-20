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
    a_iterations = int(sys.argv[5]) 
    b_iterations = int(sys.argv[6]) 
    c_iterations = int(sys.argv[7]) 
    if c_iterations < 0 or b_iterations < 0 or a_iterations < 0:
        sys.stderr.write("Error: Please provide a positive integer for all repetition numbers.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length value provide. Please enter a valid numeric value.\n")
    sys.exit(1)


try:
        ##Sulfate degree ref:/10.1016/j.carbpol.2019.115292
        Sulfate_target = float(sys.argv[8])  #####
        if not (0 <= Sulfate_target <= 1.5):
            sys.stderr.write("Error: Sulfate_target should be between 0 and 1.5.\n")
            sys.exit(1)


except IndexError:
    sys.stderr.write("Error: Missing one or more input arguments. Please provide all five effective values.\n")
    sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid input. Please ensure all values are numeric and within the valid ranges.\n")
    sys.exit(1)            # pH condition



with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
neutron_pdb_1 = os.path.join(main_folder_path, 'structure', 'cellulose_II', 'charmm36', 'chain-1.pdb')
user_defined_pdb_1 = os.path.join(main_folder_path, 'structure', 'cellulose_II', 'charmm36', 'chain-1_ud.pdb')

if os.path.exists(user_defined_pdb_1):
    unit_chain_input_file_1 = user_defined_pdb_1
else:
    unit_chain_input_file_1 = neutron_pdb_1

neutron_pdb_2 = os.path.join(main_folder_path, 'structure', 'cellulose_II', 'charmm36', 'chain-2.pdb')
user_defined_pdb_2 = os.path.join(main_folder_path, 'structure', 'cellulose_II', 'charmm36', 'chain-2_ud.pdb')

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
    resid_1 = chain_1_i  * 2 - 1
    resid_2 = chain_1_i  * 2   

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
    for chain_1_output in (chain_1_strand):
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
    'O1': chain_1_first_atom.position + [-0.436,-0.505,-1.265],
    'HO1': chain_1_first_atom.position + [-0.176,0.104,-1.96],
    'HO4': chain_1_last_o4.position + [0.677,-0.27,0.625]  
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



#chain_2 assembly    
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
    'O1': chain_2_first_atom.position + [-0.081,0.726,1.23],
    'HO1': chain_2_first_atom.position + [-0.119, 0.109, 1.964],
    'HO4': chain_2_last_o4.position + [0.845, 0.243, -0.386]
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
##------------------------------------------------------------------------------
##translate parameter

angle_rad = (gamma_angle - 90) * math.pi / 180
cosA = math.cos(angle_rad)
sinA = math.sin(angle_rad)

b_par_transverse_move= -b_trans * sinA
b_par_vertical_move= 1 * b_trans * cosA


#chain_1_layer assembly  
chain_1_1st_layer_input_file = "chain_1_temp.2.pdb"
chain_1_1st_layer_u = mda.Universe(chain_1_1st_layer_input_file)
chain_1_1st_layer = []  

for i in range(1, a_iterations+1):
    chain_1_1st_layer_u.atoms.positions += [a_trans, 0, 0]

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

#chain_1_1ayer_temp=mda.Universe("chain_1_1st_layer_temp.pdb")        
#for i, atom in enumerate(chain_1_1ayer_temp.atoms, start=1):
#    atom.id = i
#with mda.Writer("chain_1_1st_layer_temp.pdb", n_atoms=chain_1_1ayer_temp.atoms.n_atoms) as W:
#    W.write(chain_1_1ayer_temp.atoms)

#chain_1_sheet assembly  
chain_1_sheet_input_file = "chain_1_1st_layer_temp.pdb"
chain_1_sheet = [] 

chain_1_base_universe = mda.Universe(chain_1_sheet_input_file)

for j in range(1, b_iterations+1):
    chain_1_sheet_u = chain_1_base_universe.copy()
    translation_vector = [b_par_transverse_move * j  ,
                          b_par_vertical_move * j  , 0]
    chain_1_sheet_u.atoms.positions += translation_vector

    segid_increment_value = (j-1)* a_iterations 
    for segid, group in chain_1_sheet_u.atoms.groupby('segids').items():
        new_segid = str(int(segid) + segid_increment_value)
        for atom in group:
            atom.segment.segid = new_segid

    chain_1_sheet_output = f"chain_1_sheet_{j}.pdb"
    chain_1_sheet.append(chain_1_sheet_output)
    with mda.Writer(chain_1_sheet_output, n_atoms=chain_1_sheet_u.atoms.n_atoms) as W:
        W.write(chain_1_sheet_u.atoms)
with open("chain_1_sheet_temp.pdb", "w") as layer_file:
    for chain_1_sheet_output in chain_1_sheet:
        with open(chain_1_sheet_output, "r") as chain_1_sheet_pdb_file:
            for line in chain_1_sheet_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
                elif line.startswith("TER"): 
                    layer_file.write("TER\n")
        os.remove(chain_1_sheet_output)  


#chain_2_layer assembly

#chain_2_1st_layer assembly  
chain_2_1st_layer_input_file = "chain_2_temp.2.pdb"
chain_2_1st_layer = []  

for i in range(1, a_iterations+1):  
    chain_2_1st_layer_u = mda.Universe(chain_2_1st_layer_input_file)
    translation_vector = [a_trans * i,
                            0, 0]

    chain_2_1st_layer_u.atoms.positions += translation_vector

    for atom in chain_2_1st_layer_u.atoms:
        new_numeric_part = b_iterations * a_iterations + i
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


#chain_2_sheet assembly  
chain_2_sheet_input_file = "chain_2_1st_layer_temp.pdb"
chain_2_sheet = [] 

chain_2_base_universe = mda.Universe(chain_2_sheet_input_file)

for j in range(1, b_iterations+1):
    chain_2_sheet_u = chain_2_base_universe.copy()

    translation_vector = [b_par_transverse_move * j  ,
                          b_par_vertical_move * j  , 0]    
    chain_2_sheet_u.atoms.positions += translation_vector

    segid_increment_value = a_iterations * (j-1)
    for segid, group in chain_2_sheet_u.atoms.groupby('segids').items():
        new_segid = str(int(segid) + segid_increment_value)
        for atom in group:
            atom.segment.segid = new_segid

    chain_2_sheet_output = f"chain_2_sheet_{j}.pdb"
    chain_2_sheet.append(chain_2_sheet_output)
    with mda.Writer(chain_2_sheet_output, n_atoms=chain_2_sheet_u.atoms.n_atoms) as W:
        W.write(chain_2_sheet_u.atoms)
with open("chain_2_sheet_temp.pdb", "w") as layer_file:
    for chain_2_sheet_output in chain_2_sheet:
        with open(chain_2_sheet_output, "r") as chain_2_sheet_pdb_file:
            for line in chain_2_sheet_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
                elif line.startswith("TER"): 
                    layer_file.write("TER\n")
        os.remove(chain_2_sheet_output)  


def assemble_pdbs(output_file, input_files):
    with open(output_file, 'w') as outfile:
        for filename in input_files:
            with open(filename, 'r') as infile:
                for line in infile:
                    if line.startswith("ATOM"):
                        outfile.write(line)
            outfile.write('TER\n')  

chain_temp_pdb_files = ["chain_1_sheet_temp.pdb", "chain_2_sheet_temp.pdb"]

unit_temp_pdb = "unit_temp.pdb"
assemble_pdbs(unit_temp_pdb, chain_temp_pdb_files)

##crystallographic planes defined
planes_100_left = []
planes_100_right = []


planes_010_upper = []
planes_010_bottom = []
planes_100_left_chain_1 = []
planes_100_left_chain_2 = []
planes_100_right_chain_1 = []
planes_100_right_chain_2 = []

range1_start = 1
range1_end = a_iterations-2
range2_end = int((a_iterations * b_iterations)*2 -2)
range2_start =int(range2_end - a_iterations + 3)

planes_010_bottom.extend([range1_start, range1_end])
planes_010_upper.extend([range2_start, range2_end])


for m in range(1, b_iterations*2+1):
    value1 = (m-1)*(a_iterations)
    value2 = m*(a_iterations)-1
    planes_100_left.extend([value1])
    planes_100_right.extend([value2]) 
    if m <= b_iterations:  
        planes_100_left_chain_1.extend([value1])
        planes_100_right_chain_1.extend([value2])
    else:  
        planes_100_left_chain_2.extend([value1])
        planes_100_right_chain_2.extend([value2])


#print("001 upper", planes_001_upper)
#print("001 bottom", planes_001_bottom)
#print("100 right", planes_100_right)
#print("100 left", planes_100_left)
#print("100 right_chain1", planes_100_right_chain_1)
#print("100 left_chain1", planes_100_left_chain_1 )
#print("100 right_chain2", planes_100_right_chain_2)
#print("100 left_chain2", planes_100_left_chain_2 )


planes_010_count=(planes_010_upper[-1] - planes_010_upper[0] + 1) + (planes_010_bottom[-1] - planes_010_bottom[0] + 1)
planes_100_count=len(planes_100_left)+len(planes_100_right)




#--------------so3 modified step ------------------------------------------
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
#print('sulfate_num',n_so3)

#


sulfate_num = n_so3


if sulfate_num > 0:
    left_100_count=0
    right_100_count=0
    upper_010_count=0
    bottom_010_count=0
    for i in range(1, sulfate_num + 1):
        # Determine the target based on i % 4
        if i % 4 == 0:
            bottom_010_count += 1
        elif i % 4 == 1:
            upper_010_count += 1
        elif i % 4 == 2:
            left_100_count += 1
        elif i % 4 == 3:
            right_100_count += 1


    planes_010_upper_count  = (planes_010_upper[-1] - planes_010_upper[0] + 1) 
    planes_010_bottom_count = (planes_010_bottom[-1] - planes_010_bottom[0] + 1)
    planes_100_left_count   =  len(planes_100_left)
    planes_100_right_count  =  len(planes_100_right)

    max_100_count_left = c_iterations * planes_100_left_count 
    max_100_count_right = c_iterations * planes_100_right_count
    

    max_010_count_upper = c_iterations  * planes_010_upper_count
    max_010_count_bottom = c_iterations * planes_010_bottom_count


    left_100_count_update = left_100_count
    right_100_count_update = right_100_count
    
    upper_010_count_update = upper_010_count
    bottom_010_count_update = bottom_010_count

    if left_100_count > max_100_count_left:
        left_100_count_update = max_100_count_left
    
    if right_100_count > max_100_count_right:
        right_100_count_update = max_100_count_right
    
    if upper_010_count > max_010_count_upper:
        upper_010_count_update = max_010_count_upper
    
    if bottom_010_count > max_010_count_bottom:
        bottom_010_count_update = max_010_count_bottom    
        
    new_sulfate_number_update=upper_010_count_update + bottom_010_count_update + left_100_count_update + right_100_count_update 
    actual_ds=(new_sulfate_number_update)/(total_residues)
    actual_sulfate_content=Sulfate_target*(actual_ds/ds)
    #print(actual_ds)
    #print(actual_sulfate_content)
    actual_sulfate_rounded=round(actual_sulfate_content,4)
    #print('surface charge density', actual_sulfate_rounded, 'mmol/g')



    even_numbers_100 = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 0]
    odd_numbers_100 = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 1]

    all_combinations_100_left  = [(value1, value2) for value1 in planes_100_left for value2 in even_numbers_100 ]
    all_combinations_100_right = [(value1, value2) for value1 in planes_100_right  for value2 in odd_numbers_100 ]
    random.shuffle(all_combinations_100_left)
    random.shuffle(all_combinations_100_right)
    sel_100_left_array = []
    sel_100_right_array = []

    while len(sel_100_left_array) < left_100_count_update:
        if not all_combinations_100_left:  
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_100_left = random.choice(all_combinations_100_left)
        sel_100_left_array.append(selected_combination_100_left)
        all_combinations_100_left.remove(selected_combination_100_left)
    
    while len(sel_100_right_array) < right_100_count_update:
        if not all_combinations_100_right:  # Check if the list is empty to avoid IndexError
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_100_right = random.choice(all_combinations_100_right)
        sel_100_right_array.append(selected_combination_100_right)
        all_combinations_100_right.remove(selected_combination_100_right)
        
    
   
    numbers_010_odd = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 1]
    numbers_010_even = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 0]
    all_combinations_010_upper  = [(value1, value2) for value1 in range(min(planes_010_upper),  max(planes_010_upper)+ 1)  for value2 in numbers_010_even ]
    all_combinations_010_bottom = [(value1, value2) for value1 in range(min(planes_010_bottom), max(planes_010_bottom) + 1)  for value2 in numbers_010_odd ]

    random.shuffle(all_combinations_010_upper)
    random.shuffle(all_combinations_010_bottom)
    sel_010_upper_array = []
    sel_010_bottom_array = []
        
    while len(sel_010_upper_array) < upper_010_count_update:
        if not all_combinations_010_upper:  
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_010_upper = random.choice(all_combinations_010_upper)
        sel_010_upper_array.append(selected_combination_010_upper)
        all_combinations_010_upper.remove(selected_combination_010_upper)

    
    while len(sel_010_bottom_array) < bottom_010_count_update:
        if not all_combinations_010_bottom:  # Check if the list is empty to avoid IndexError
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_010_bottom = random.choice(all_combinations_010_bottom)
        sel_010_bottom_array.append(selected_combination_010_bottom)
        all_combinations_010_bottom.remove(selected_combination_010_bottom)

    #print("100 bottom array residue:", sel_100_bottom_array)
    #print("100 upper array residue:", sel_100_upper_array )



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
    left_glc_100_candidates_indices   = get_fragment_indices(origin_glc_chains, sel_100_left_array)    
    right_glc_100_candidates_indices  = get_fragment_indices(origin_glc_chains, sel_100_right_array)     
    upper_glc_010_candidates_indices  = get_fragment_indices(origin_glc_chains, sel_010_upper_array)     
    bottom_glc_010_candidates_indices = get_fragment_indices(origin_glc_chains, sel_010_bottom_array)    
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
    
    all_indices = [left_glc_100_candidates_indices, right_glc_100_candidates_indices, 
                   upper_glc_010_candidates_indices, bottom_glc_010_candidates_indices]
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
    left_so3_candidates_indices   = get_fragment_indices(primary_hydroxyl_remove_chains, sel_100_left_array)    
    right_so3_candidates_indices  = get_fragment_indices(primary_hydroxyl_remove_chains, sel_100_right_array)   
    upper_so3_candidates_indices  = get_fragment_indices(primary_hydroxyl_remove_chains, sel_010_upper_array)   
    bottom_so3_candidates_indices = get_fragment_indices(primary_hydroxyl_remove_chains, sel_010_bottom_array)  
    primary_hydroxyl_remove_chains_u = mda.Universe("primary_hydroxyl_remove.temp.pdb")
    so3_atoms = list(primary_hydroxyl_remove_chains_u.atoms)#



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
                    if residue.resid % 2 == 1  and fragment_number in range(min(planes_010_upper), max(planes_010_upper)+1):
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([ 0.632,   -0.219,  0.834999999999994]),
                            'H62':  base_atom.position + np.array([ 0.632,    0.179,  -0.844999999999999]),
                            'O6':   base_atom.position + np.array([-0.754,   -1.182,  -0.280000000000001]),
                            'S6':   base_atom.position + np.array([ 0.279,   -2.490,  -0.391000000000005]),
                            'OS62': base_atom.position + np.array([-0.522,   -3.874,  0.088999999999999]),
                            'OS63': base_atom.position + np.array([ 0.792,   -2.662,  -1.97000000000001]),
                            'OS64': base_atom.position + np.array([ 1.600,   -2.240,  0.599999999999994])
                        }
                    elif residue.resid % 2 == 0 and fragment_number in range(min(planes_010_upper), max(planes_010_upper)+1):
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([-0.632,   0.219,  0.835000000000001]),
                            'H62':  base_atom.position + np.array([-0.632,  -0.179,  -0.844999999999999]),
                            'O6':   base_atom.position + np.array([ 0.354,   1.182,  -0.280000000000001]),
                            'S6':   base_atom.position + np.array([-0.253,   2.513,  -0.222000000000001]),
                            'OS62': base_atom.position + np.array([-1.494,   2.225,  0.856999999999999]),
                            'OS63': base_atom.position + np.array([-0.889,   2.796,  -1.74]),
                            'OS64': base_atom.position + np.array([ 0.618,   3.849,  0.274000000000001])
                        }
    
                    elif residue.resid % 2 == 1  and fragment_number in range(min(planes_010_bottom), max(planes_010_bottom)+1):
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([ 0.632,   -0.219,  0.834999999999994]),
                            'H62':  base_atom.position + np.array([ 0.632,    0.179,  -0.844999999999999]),
                            'O6':   base_atom.position + np.array([-0.754,   -1.182,  -0.280000000000001]),
                            'S6':   base_atom.position + np.array([ 0.279,   -2.490,  -0.391000000000005]),
                            'OS62': base_atom.position + np.array([-0.522,   -3.874,  0.088999999999999]),
                            'OS63': base_atom.position + np.array([ 0.792,   -2.662,  -1.97000000000001]),
                            'OS64': base_atom.position + np.array([ 1.600,   -2.240,  0.599999999999994])
                        }
                    elif residue.resid % 2 == 0 and fragment_number in range(min(planes_010_bottom), max(planes_010_bottom)+1):
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([-0.632,   0.219,  0.835000000000001]),
                            'H62':  base_atom.position + np.array([-0.632,  -0.179,  -0.844999999999999]),
                            'O6':   base_atom.position + np.array([ 0.754,   1.182,  -0.280000000000001]),
                            'S6':   base_atom.position + np.array([-0.253,   2.513,  -0.222000000000001]),
                            'OS62': base_atom.position + np.array([-1.494,   2.225,  0.856999999999999]),
                            'OS63': base_atom.position + np.array([-0.889,   2.796,  -1.74]),
                            'OS64': base_atom.position + np.array([ 0.618,   3.849,  0.274000000000001])
                        }
    
                    elif fragment_number in planes_100_right:  # Odd residue ID
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([-0.131,	-0.715,	 0.784]),
                            'H62':  base_atom.position + np.array([ 0.308,	-0.503,	-0.893]),
                            'O6':   base_atom.position + np.array([ 0.998,	 0.949,	 0.381]),
                            'S6':   base_atom.position + np.array([ 2.485,	 0.192,	 0.466]),
                            'OS62': base_atom.position + np.array([ 3.352,	 0.511,	-0.925]),
                            'OS63': base_atom.position + np.array([ 3.327,	 0.763,	 1.790]),
                            'OS64': base_atom.position + np.array([ 2.262,	-1.456,	 0.616])
                        }
                    elif fragment_number in planes_100_left:  # Even residue ID
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([ 0.131,	 1.049,	-0.167]),
                            'H62':  base_atom.position + np.array([-0.308,	-0.164,	 1.011]),
                            'O6':   base_atom.position + np.array([-0.999,	-0.504,	-0.891]),
                            'S6':   base_atom.position + np.array([-2.484,	-0.493,	-0.127]),
                            'OS62': base_atom.position + np.array([-3.518,	 0.541,	-0.933]),
                            'OS63': base_atom.position + np.array([-2.295,	 0.026,	 1.449]),
                            'OS64': base_atom.position + np.array([-3.123,	-2.036,	-0.132])
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


    u_final = mda.Universe("so3_negative-modified.temp.pdb")
    box_x = a_iterations * a_trans+20
    box_y = b_iterations * b_trans+20
    box_z = c_iterations * c_trans+20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z,90,90,90]
    translation_vector = center_of_box - center_of_mass
    u_final.atoms.translate(translation_vector)
    with mda.Writer(f"charmm36-cellulose-II-fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
        
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    
    print(f"Half-Ester sulfate content: {actual_sulfate_content:.2f}, Degree of sulfate: {actual_ds:.4f}")


elif sulfate_num == 0:
    u_final = mda.Universe("unit_temp.pdb")
    box_x = a_iterations * a_trans+20
    box_y = b_iterations * b_trans+20
    box_z = c_iterations * c_trans+20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z,90,90,90]
    translation_vector = center_of_box - center_of_mass
    with mda.Writer(f"charmm36-cellulose-II-fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
    

    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    actual_ds=0.0    
    actual_sulfate_content=0.0
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"Half-Ester sulfate content: {actual_sulfate_content:.2f}, Degree of sulfate: {actual_ds:.4f}")
