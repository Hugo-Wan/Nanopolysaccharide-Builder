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



#-------------------------------------------- modified parameter  --------------------------------------------------
try:
    a_trans = float(sys.argv[1])
    b_trans = float(sys.argv[2])
    c_trans = float(sys.argv[3])
    if a_trans < 4 or b_trans < 15 or c_trans < 10:
        sys.stderr.write("Error: Please provide an effective values for all crystallographic parameters.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length value provided. Please enter valid numeric values.\n")
    sys.exit(1)

try:
    a_iterations = int(sys.argv[4])
    b_iterations = int(sys.argv[5])
    c_iterations = int(sys.argv[6])
    
    if a_iterations <= 0 or b_iterations <= 0 or c_iterations <= 0:
        sys.stderr.write("Error: Please provide positive integers for all repetition numbers.\n")
        sys.exit(1)
    ##deacetylation degree ref:10.1021/acs.jchemed.7b00902J
    DDA_target = float(sys.argv[7])
    if not (0 <= DDA_target < 1):
        sys.stderr.write("Error: Deacetylation degree (DDA) must be between 0 and 1.\n") 
        sys.exit(1)
    pH = float(sys.argv[8])
    if not (0 <= pH <= 14):
        sys.stderr.write("Error: pH must be between 0 and 14.\n")
        sys.exit(1)
    #print(f"Validated input values:\n  a_iterations: {a_iterations}, b_iterations: {b_iterations}, c_iterations: {c_iterations}\n  DDA_target: {DDA_target}, pH: {pH}")

except IndexError:
    sys.stderr.write("Error: Missing one or more input arguments. Please provide all five effective values.\n")
    sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid input. Please ensure all values are numeric and within the valid ranges.\n")
    sys.exit(1)  

pKa=6.3                        #pka of chitosan  ref:doi.org/10.1016/j.foodhyd.2022.108383
m_GlcNAc = 204.20              # GlcNAc unit mass
m_Glc = 179.17                 # Glc unit mass
#-------------------------------------------- modified parameter  --------------------------------------------------


with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
neutron_pdb_1 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36', 'A_configuration', 'left-unit.pdb')
user_defined_pdb_1 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36', 'A_configuration', 'left-unit_ud.pdb')

if os.path.exists(user_defined_pdb_1):
    unit_chain_input_file_1 = user_defined_pdb_1
else:
    unit_chain_input_file_1 = neutron_pdb_1

neutron_pdb_2 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36', 'A_configuration', 'right-unit.pdb')
user_defined_pdb_2 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36', 'A_configuration', 'right-unit_ud.pdb')

if os.path.exists(user_defined_pdb_2):
    unit_chain_input_file_2 = user_defined_pdb_2
else:
    unit_chain_input_file_2 = neutron_pdb_2

l_u = mda.Universe(unit_chain_input_file_1)
r_u = mda.Universe(unit_chain_input_file_2)



 #left_chain assembly    
# ------------------------------------------------   
left_chain = []
for l_i in range(1, c_iterations + 1):
    l_u.atoms.positions += [0, 0, c_trans]
    resid_1 = (c_iterations - l_i + 1) * 2 - 1
    resid_2 = (c_iterations - l_i + 1) * 2
    for j, atom in enumerate(l_u.atoms):
        if j < 27:
            atom.residue.resid = resid_1
        else:
            atom.residue.resid = resid_2
    left_chain_output = f"l-{l_i}.pdb"
    left_chain.append(left_chain_output)
    with mda.Writer(left_chain_output, n_atoms=l_u.atoms.n_atoms) as W:
        W.write(l_u.atoms)
with open("left-chain_temp.pdb", "w") as l_chain:
    for l_chain_output in reversed(left_chain):
        with open(l_chain_output, "r") as lef_chain_pdb_file:
            for left_chain_line in lef_chain_pdb_file:
                if left_chain_line.startswith("ATOM"):
                    l_chain.write(left_chain_line)
            os.remove(l_chain_output) 

 ##terminal side            
left_finite_build_chain = "left-chain_temp.pdb"
left_finite_u = mda.Universe(left_finite_build_chain)
left_first_atom = left_finite_u.atoms[0]
o4_atoms = left_finite_u.select_atoms("name O4")
if not o4_atoms:
    raise ValueError("No atoms named 'O4' found in the structure.")
last_o4 = o4_atoms[-1]

new_positions = {
    'O1': left_first_atom.position + [-0.538, -0.449, 1.247],
    'HO1': left_first_atom.position + [0.129, -0.365, 1.932],
    'HO4': last_o4.position + [0.377, -0.192, -0.892]
}

left_new_atoms = {}
for name, pos in new_positions.items():
    left_base_atom = left_first_atom if name in ['O1', 'HO1'] else last_o4
    left_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
    left_new_uni.add_TopologyAttr('name', [name])
    left_new_uni.add_TopologyAttr('type', [left_base_atom.type])
    left_new_uni.add_TopologyAttr('resname', [left_base_atom.resname])
    left_new_uni.add_TopologyAttr('resid', [left_base_atom.resid])
    left_new_uni.add_TopologyAttr('chainIDs', ['N'])
    left_new_uni.add_TopologyAttr('segid', ['0']) 
    left_new_uni.atoms.positions = [pos]
    left_new_atoms[name] = left_new_uni

left_combined = Merge(left_finite_u.atoms[:last_o4.index + 1], left_new_atoms['HO4'].atoms, left_finite_u.atoms[last_o4.index + 1:])
left_combined = Merge(left_combined.atoms[:2], left_new_atoms['O1'].atoms, left_combined.atoms[2:])
left_combined = Merge(left_combined.atoms[:3], left_new_atoms['HO1'].atoms, left_combined.atoms[3:])

for atom in left_combined.atoms:
    atom.segment.segid = '0'  
left_combined.atoms.write("left-chain_temp.2.pdb")
#------------------------------------------------------------------------------




#right_chain assembly      
#------------------------------------------------------------------------------ 
right_chain = []
for r_i in range(1, c_iterations + 1):
    r_u.atoms.positions += [0, 0, c_trans]
    resid_1 = r_i  * 2 - 1
    resid_2 = r_i  * 2
    for j, atom in enumerate(r_u.atoms):
        if j < 27:
            atom.residue.resid = resid_1
        else:
            atom.residue.resid = resid_2
    right_chain_output = f"r-{r_i}.pdb"
    right_chain.append(right_chain_output)
    with mda.Writer(right_chain_output, n_atoms=r_u.atoms.n_atoms) as W:
        W.write(r_u.atoms)
with open("right-chain_temp.pdb", "w") as chain:
    for r_chain_output in right_chain:
         with open(r_chain_output, "r") as right_chain_pdb_file:
            for right_chain_line in right_chain_pdb_file:
                if right_chain_line.startswith("ATOM"):
                    chain.write(right_chain_line)
            os.remove(r_chain_output)

##terminal side 
right_finite_build_chain = "right-chain_temp.pdb"
right_finite_u = mda.Universe(right_finite_build_chain)
right_first_atom = right_finite_u.atoms[0]
o4_atoms = right_finite_u.select_atoms("name O4")
if not o4_atoms:
    raise ValueError("No atoms named 'O4' found in the structure.")
last_o4 = o4_atoms[-1]

new_positions = {
    'O1': right_first_atom.position - [-0.538, -0.449, 1.247],
    'HO1': right_first_atom.position - [0.129, -0.365, 1.932],
    'HO4': last_o4.position - [0.377, -0.192, -0.892]
}
right_side_start_segid=a_iterations
right_new_atoms = {}
for name, pos in new_positions.items():
    right_base_atom = right_first_atom if name in ['O1', 'HO1'] else last_o4
    right_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
    right_new_uni.add_TopologyAttr('name', [name])
    right_new_uni.add_TopologyAttr('type', [right_base_atom.type])
    right_new_uni.add_TopologyAttr('resname', [right_base_atom.resname])
    right_new_uni.add_TopologyAttr('resid', [right_base_atom.resid])
    right_new_uni.add_TopologyAttr('chainIDs', ['N']) 
    right_new_uni.add_TopologyAttr('segid', ['CARB']) 
    right_new_uni.atoms.positions = [pos]
    right_new_atoms[name] = right_new_uni

right_combined = Merge(right_finite_u.atoms[:last_o4.index + 1], right_new_atoms['HO4'].atoms, right_finite_u.atoms[last_o4.index + 1:])
right_combined = Merge(right_combined.atoms[:2], right_new_atoms['O1'].atoms, right_combined.atoms[2:])
right_combined = Merge(right_combined.atoms[:3], right_new_atoms['HO1'].atoms, right_combined.atoms[3:])

for atom in right_combined.atoms:
    atom.segment.segid = f'{right_side_start_segid}' 
right_combined.atoms.write("right-chain_temp.2.pdb")


#------------------------------------------------------------------------------


 #left_layer assembly  
left_layer_input_file = "left-chain_temp.2.pdb"
layer_l_u = mda.Universe(left_layer_input_file)
left_layer = []  
for i in range(1, a_iterations + 1):
    layer_l_u.atoms.positions += [a_trans, 0, 0]

    segid_increment_value = 1 
    for atom in layer_l_u.atoms:
        numeric_part = atom.segid.replace('','')
        new_numeric_part = int(numeric_part) + segid_increment_value
    new_segid = f"{new_numeric_part}"
    atom.segment.segid = new_segid

    left_layer_output = f"left_layer_{i}.pdb"
    left_layer.append(left_layer_output)
    with mda.Writer(left_layer_output, n_atoms=layer_l_u.atoms.n_atoms) as W:
        W.write(layer_l_u.atoms)
with open("left-layer_temp.pdb", "w") as layer_file:
    for left_layer_output in left_layer:
        with open(left_layer_output, "r") as left_layer_pdb_file:
            for left_layer_line in left_layer_pdb_file:
                if left_layer_line.startswith("ATOM"):
                    layer_file.write(left_layer_line)
            os.remove(left_layer_output) 


 #right_layer assembly    
right_layer_input_file = "right-chain_temp.2.pdb"
layer_r_u = mda.Universe(right_layer_input_file)
right_layer = []  
for i in range(1, a_iterations + 1):
    layer_r_u .atoms.positions += [a_trans, 0, 0]
    #j=i+a_iterations
    segid_increment_value = 1
    for atom in layer_r_u.atoms:
        numeric_part = atom.segid.replace('','')
        new_numeric_part = int(numeric_part) + segid_increment_value
    new_segid = f"{new_numeric_part}"
    atom.segment.segid = new_segid


    right_layer_output = f"right_layer_{i}.pdb"
    right_layer.append(right_layer_output)
    with mda.Writer(right_layer_output, n_atoms=layer_r_u .atoms.n_atoms) as W:
        W.write(layer_r_u .atoms)
with open("right-layer_temp.pdb", "w") as layer_file:
    for right_layer_output in right_layer:
        with open(right_layer_output, "r") as right_layer_pdb_file:
            for right_layer_line in right_layer_pdb_file:
                if right_layer_line.startswith("ATOM"):
                    layer_file.write(right_layer_line)
            os.remove(right_layer_output) 


 #unit_structure assembly  
left_layer_input_file = "left-layer_temp.pdb"
right_layer_input_file = "right-layer_temp.pdb"
unit_output_file = "unit.pdb" 
with open(unit_output_file, "w") as unit_file:
    with open("left-layer_temp.pdb", "r") as left_file:
        for line in left_file:
            if line.startswith("ATOM"):
                unit_file.write(line)
    with open("right-layer_temp.pdb", "r") as right_file:
        for line in right_file:
            if line.startswith("ATOM"):
                unit_file.write(line)

unit_input_file="unit.pdb"
crystal_structure = mda.Universe(unit_input_file)
crystal_structure_files = [] 

for crystal_i in range(1, b_iterations + 1):
    crystal_structure.atoms.positions += [0, -b_trans, 0]
    crystal_structure_output = f"crystal_{crystal_i}.pdb"
    crystal_structure_files.append(crystal_structure_output)
    with mda.Writer(crystal_structure_output, n_atoms=crystal_structure.atoms.n_atoms) as writer:
        writer.write(crystal_structure.atoms)


for file_index, file_name in enumerate(crystal_structure_files):
    temp_structure = mda.Universe(file_name)
    segid_increment_value = file_index * 2 * a_iterations
    for segment in temp_structure.segments:
        matches = re.findall(r'\d+', segment.segid)
        if matches:
            numeric_part = int(matches[0])
            new_numeric_part = numeric_part + segid_increment_value
            new_segid = f"{new_numeric_part}"
            for atom in segment.atoms:
                atom.segment.segid = new_segid
    with mda.Writer(file_name, n_atoms=temp_structure.atoms.n_atoms) as writer:
        writer.write(temp_structure.atoms)
with open("alpha-chitin-A-temp.pdb", "w") as combined_file:
    for file_name in crystal_structure_files:
        with open(file_name, "r") as file:
            for line in file:
                if line.startswith("ATOM"):
                    combined_file.write(line)
        os.remove(file_name) 

os.remove(unit_input_file) 


#-------------------------------------------------- native chitin built ----------------------------------------------



##crystallographic planes defined
planes_010_left = []
planes_010_right = []

planes_100_even = []
planes_100_odd = []

planes_100_upper = []
planes_100_bottom = []



###m is odd is for right-chain, even is for left-chain structure
for m in range(1, b_iterations*2-1):
    value1 = m*a_iterations
    value2 = (m+1)*a_iterations - 1
    planes_100_bottom.extend([value1])
    planes_100_upper.extend([value2]) 
    if m % 2 == 0:
        planes_100_even.extend([value1, value2])
    else:
        planes_100_odd.extend([value1, value2])

range1_start = 0
range1_end = a_iterations-1
range2_start = b_iterations*2*a_iterations-a_iterations
range2_end = b_iterations*2*a_iterations-1

planes_010_left.extend([range1_start, range1_end])
planes_010_right.extend([range2_start, range2_end])



planes_010_count=(planes_010_left[-1] - planes_010_left[0] + 1) + (planes_010_right[-1] - planes_010_right[0] + 1)
planes_100_count=len(planes_100_bottom)+len(planes_100_upper)
modified_chains_count=planes_010_count+planes_100_count


#print(planes_100_count)
#print("100_planes_even:",   planes_100_even)
#print("100_planes_odd:",    planes_100_odd)
#print("100_planes_upper:",  planes_100_upper)
#print("100_planes_bottom:", planes_100_bottom)
#
#print("010_planes_r:", planes_010_right)
#print("010_planes_l:", planes_010_left)


##------------------------------------------------- nh3 modified step -------------------------------------------------------
#
#
############NH3 located step

chain_u = mda.Universe("alpha-chitin-A-temp.pdb")
residue_ids = [residue.resid for residue in chain_u.residues]
total_residues=len(residue_ids)
#
##print("Left Modified Chains:", left_modified_chains)
##print("Right Modified Chains:", right_modified_chains)
#
#
###NHx group evenly assigned to each chain
##define the dda; how many resiudes should be modified:
def calculate_dda(x):
    numerator = (total_residues - x) * m_Glc
    denominator = numerator + x * m_GlcNAc
    return numerator / denominator if denominator != 0 else 0

#def adjust_closest_x_for_positive_integer(closest_x, total_residues):
#    remainder = (total_residues - closest_x) % modified_chains_count
#    if remainder != 0:
#        closest_x += remainder
#    return max(0, closest_x)

def adjust_for_amination_limit(closest_x, total_residues):
    amination_number = total_residues - closest_x
    amination_number_limit=c_iterations*planes_010_count + 2*c_iterations*planes_100_count
    if amination_number > amination_number_limit:
        # Reduce amination_number to 2 * c * a_iterations
        new_amination_number = amination_number_limit
        closest_x = total_residues - new_amination_number
    return closest_x

closest_x = None
closest_dda = float('inf')
for x in range(total_residues + 1):
    dda = calculate_dda(x)
    if abs(dda - DDA_target) < abs(closest_dda - DDA_target):
        closest_x = x
        closest_dda = dda

closest_x = adjust_for_amination_limit(closest_x, total_residues)


if DDA_target == 0:
    closest_x = total_residues

aminiation_number = total_residues - closest_x
#



if aminiation_number > 0:
    
    
    left_010_count=0
    upper_100_count=0

    for i in range(1, aminiation_number + 1):
        # Determine the target based on i % 4
        if i % 4 == 1:
            upper_100_count += 1
        elif i % 4 == 2:
            left_010_count += 1

    
    new_aminiation_number=aminiation_number
    
    planes_010_left_count=(planes_010_left[-1] - planes_010_left[0] + 1)
    planes_100_upper_count=len(planes_100_upper)

    
    
    
    max_010_count_left = c_iterations * planes_010_left_count 

    
    
    max_100_count_upper = c_iterations * 2 * planes_100_upper_count

    
    
    left_010_count_update = left_010_count

    
    upper_100_count_update = upper_100_count

        
    if left_010_count > max_010_count_left:
        left_010_count_update = max_010_count_left
    

    if upper_100_count > max_100_count_upper:
        upper_100_count_update = max_100_count_upper
    

    
    new_aminiation_number_update=upper_100_count_update + left_010_count_update 
    actual_dda=(new_aminiation_number_update*179.17)/(new_aminiation_number_update*179.17 + (total_residues-new_aminiation_number_update)*204.20)
    #print(f"Actual DDA is : {actual_dda}")
    #print(f"aminiation_number: {new_aminiation_number_update}")
    #print(f"aminiation_number_per_chain: {modified_num_per_chain}")
    actual_dda_rounded=round(actual_dda,4)
    #print(actual_dda_rounded)
    
    ##010 planes selection
    even_numbers_010 = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 0]
    odd_numbers_010 = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 1]
    
    all_combinations_010_left = [(value1, value2) for value1 in range(min(planes_010_left), max(planes_010_left) + 1)  for value2 in even_numbers_010 ]
    all_combinations_010_right = [(value1, value2) for value1 in range(min(planes_010_right), max(planes_010_right) + 1)  for value2 in odd_numbers_010 ]
    #print(len(all_combinations_010_left))
    random.shuffle(all_combinations_010_left)
    random.shuffle(all_combinations_010_right)
    sel_010_left_array = []
    sel_010_right_array = []
    
    
    while len(sel_010_left_array) < left_010_count_update:
        if not all_combinations_010_left:  
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_010_left = random.choice(all_combinations_010_left)
        sel_010_left_array.append(selected_combination_010_left)
        all_combinations_010_left.remove(selected_combination_010_left)
    
    

    
   
    ###100 plane selection
    #planes_100_upper = []
    #planes_100_bottom = []
    numbers_100 = [x for x in range(1, 2 * c_iterations + 1)]
    
    all_combinations_100_upper = [(value1, value2) for value1 in planes_100_upper  for value2 in numbers_100 ]
    #print(all_combinations_100_upper)
    all_combinations_100_bottom = [(value1, value2) for value1 in planes_100_bottom  for value2 in numbers_100 ]
    #print(all_combinations_100_bottom)
    
    
    #print(len(all_combinations_010_left))
    random.shuffle(all_combinations_100_upper)
    random.shuffle(all_combinations_100_bottom)
    sel_100_upper_array = []
    sel_100_bottom_array = []
    
    while len(sel_100_upper_array) < upper_100_count_update:
        if not all_combinations_100_upper:  
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_100_upper = random.choice(all_combinations_100_upper)
        sel_100_upper_array.append(selected_combination_100_upper)
        all_combinations_100_upper.remove(selected_combination_100_upper)
    #print("upper array containing unique small arrays:", sel_100_upper_array)
    #print("Length of the array:", len(sel_100_upper_array) )
    


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
    
    origin_nacetyl_chains = molecule.load("pdb", "alpha-chitin-A-temp.pdb")
    left_nacetyl_candidates_indices   = get_fragment_indices(origin_nacetyl_chains, sel_010_left_array)   
    upper_nacetyl_candidates_indices  = get_fragment_indices(origin_nacetyl_chains, sel_100_upper_array)    
    molecule.delete(origin_nacetyl_chains) 

    
    def remove_atoms_based_on_indices(universe, index_ranges):
        all_removed_indices = []
    
        for index_range in index_ranges:
            for start_idx, end_idx in index_range:
                atom_group = universe.atoms[start_idx:end_idx + 1]
                for residue in atom_group.residues:
                    if residue.resid == 1:
                        indices_to_remove = residue.atoms[10:17].indices
                    else:
                        indices_to_remove = residue.atoms[8:15].indices
                    all_removed_indices.extend(indices_to_remove)
        unique_indices = sorted(set(all_removed_indices))
        indices_str = ' '.join(map(str, unique_indices))
        universe = universe.select_atoms(f"not index {indices_str}")
    
        return universe
    
    
    
    all_indices = [left_nacetyl_candidates_indices, upper_nacetyl_candidates_indices]
    ###chain_u is the unmodified structure, with name of "alpha-chitin.temp.pdb"
    nacetyl_remove_universe = remove_atoms_based_on_indices(chain_u, all_indices)
    nacetyl_remove_universe.atoms.write("nacetyl_remove.temp.pdb")
    


    def get_fragment(mol_id, indices_i, indices_j):
        selection_query = f"index {indices_i} to {indices_j}"
        selection = atomsel(selection_query, molid=mol_id)
        fragment_num_indice = selection.fragment
        unique_fragments = set(fragment_num_indice)
        if len(unique_fragments) == 1:
                return unique_fragments.pop()   
        else:
            raise ValueError("Multiple unique fragments found in selection")
    
    nacetyl_remove_chains = molecule.load("pdb", "nacetyl_remove.temp.pdb")
    left_nh3_candidates_indices   = get_fragment_indices(nacetyl_remove_chains, sel_010_left_array)     
    upper_nh3_candidates_indices  = get_fragment_indices(nacetyl_remove_chains, sel_100_upper_array)   


    nacetyl_remove_chains_u = mda.Universe("nacetyl_remove.temp.pdb")
    nh3_atoms = list(nacetyl_remove_chains_u.atoms)

    
    #print(left_nh3_candidates_indices)
    def add_nh3_based_on_indices(universe, indices_lists):
        all_new_atoms = []  # This will store all newly created atoms across all modifications
        for index_ranges in indices_lists:
            for start_idx, end_idx in index_ranges:
                residues = universe.select_atoms(f"index {start_idx}:{end_idx}").residues
                fragment_number = get_fragment(nacetyl_remove_chains, start_idx, end_idx)
                #print("staring from to the end:",start_idx, end_idx)
                #print("fragment number is :", fragment_number )
                #print("residue is : ", residues)
                for residue in residues:
                    #print("residue is : ", residue)
                    residue.resname = 'BLNP'
                    base_atom_index = 9 if residue.resid == 1 else 7
                    base_atom = residue.atoms[base_atom_index]
                    nh3_new_positions = {} 
                    # Define new positions based on residue ID's odd/even status
                    if residue.resid % 2 == 1  and fragment_number in planes_100_even:  # Odd residue ID
                        nh3_new_positions = {
                            'HN1': base_atom.position + np.array([ 0.894,   -0.056,   0.444]),
                            'HN2': base_atom.position + np.array([-0.708,   -0.305,   0.637]),
                            'HN3': base_atom.position + np.array([-0.007,   -0.585,  -0.811])
                        }
                    elif residue.resid % 2 == 0  and fragment_number in planes_100_even:  # Even residue ID
                        nh3_new_positions = {
                            'HN1': base_atom.position + np.array([-0.806,  0.338, -0.486]),
                            'HN2': base_atom.position + np.array([ 0.791,  0.567, -0.230]),
                            'HN3': base_atom.position + np.array([-0.164,  0.040,  0.986])
                        }
    
                    elif residue.resid % 2 == 1  and fragment_number in planes_100_odd:  # Odd residue ID
                        nh3_new_positions = {
                            'HN1': base_atom.position + np.array([-0.828,   -0.326,   0.456]),
                            'HN2': base_atom.position + np.array([ 0.774,   -0.574,   0.265]),
                            'HN3': base_atom.position + np.array([-0.126,   -0.046,  -0.991])
                        }
                    elif residue.resid % 2 == 0  and fragment_number in planes_100_odd:  # Even residue ID
                        nh3_new_positions = {
                            'HN1': base_atom.position + np.array([ 0.827,  0.326,  0.458]),
                            'HN2': base_atom.position + np.array([-0.775,  0.574,  0.263]),
                            'HN3': base_atom.position + np.array([ 0.127,  0.046, -0.991])
                        }
    
                    elif fragment_number in range(min(planes_010_left), max(planes_010_left)+1):  # Odd residue ID
                        nh3_new_positions = {
                            'HN1': base_atom.position + np.array([-0.805,    0.339,  -0.487]),
                            'HN2': base_atom.position + np.array([ 0.791,    0.567,  -0.229]),
                            'HN3': base_atom.position + np.array([-0.165,    0.040,   0.986])
                        }
                    elif fragment_number in range(min(planes_010_right), max(planes_010_right)+1):  # Even residue ID
                        nh3_new_positions = {
                            'HN1': base_atom.position + np.array([-0.892, -0.056, -0.449]),
                            'HN2': base_atom.position + np.array([ 0.711, -0.307, -0.633]),
                            'HN3': base_atom.position + np.array([ 0.002, -0.584,  0.812])
                        }
                    insert_pos = nh3_atoms.index(base_atom) + 1
                    for name, pos in nh3_new_positions.items():
                        nh3_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
                        nh3_new_uni.add_TopologyAttr('name', [name])
                        nh3_new_uni.add_TopologyAttr('type', [base_atom.type])
                        nh3_new_uni.add_TopologyAttr('resname', [base_atom.resname])
                        nh3_new_uni.add_TopologyAttr('resid', [residue.resid])
                        nh3_new_uni.add_TopologyAttr('segid', [base_atom.segment.segid])
                        nh3_new_uni.add_TopologyAttr('chainIDs', ['N'])
                        nh3_new_uni.atoms.positions = [pos]
                        nh3_new_atom = nh3_new_uni.atoms[0]
    
                        nh3_atoms.insert(insert_pos, nh3_new_atom)
                        insert_pos += 1  # Update position for the next atom
        # Create a new universe with all atoms including the added ones
        new_universe = mda.Merge(*[mda.AtomGroup([atom]) for atom in nh3_atoms])
        return new_universe
    
    
    all_indices = [left_nh3_candidates_indices,  upper_nh3_candidates_indices]
    nh3_u = add_nh3_based_on_indices(nacetyl_remove_chains_u, all_indices)
    nh3_u.atoms.write("nh3-modified.temp.pdb")
    molecule.delete(nacetyl_remove_chains)


    #    
    ####------------------------------------------------------------nh2 building step----------------------------------------------------------------
    
    
    #pH calculation |  NH2_step
    ratio = 10 ** (pH - pKa)
    #print("ratio is :", ratio)
    
    #print(new_aminiation_number_update)
    y = (new_aminiation_number_update) / (1 + ratio) ###nh3 number
    #print(y)
    y_rounded = round(y) ###nh3 number for integer number
    #print(y_rounded)
    
    
    if  y_rounded == new_aminiation_number_update:
        nh3_number = new_aminiation_number_update
        nh2_number = 0
        actual_pH_rounded = pH
    elif y_rounded == 0:
        nh3_number = 0
        nh2_number = new_aminiation_number_update
        actual_pH_rounded = pH
    elif y_rounded >0 and y_rounded < new_aminiation_number_update:
        nh3_number = y_rounded
        nh2_number = int(new_aminiation_number_update - y_rounded)
        exponential_part= nh2_number/nh3_number
        actual_pH_rounded = math.log10(exponential_part) + 6.3
    
    
    #print(f"nh2 number is:{nh2_number}, nh3 number is {nh3_number}")
    #print("actual pH value is:",actual_pH_rounded )
    
    
    
    def collect_nh3_atoms_to_remove_randomly(universe, indices_groups, nh2_number):
        all_indices = []
        atoms_to_remove = []  
        for group in indices_groups:
            all_indices.extend(group)
        np.random.shuffle(all_indices)
        for start_idx, end_idx in all_indices[:nh2_number]:
            possible_residues = universe.select_atoms(f"index {start_idx}:{end_idx}").residues
            if nh2_number ==0 :
                return []
            else:
                for residue in possible_residues:
                    residue.resname = 'BLND'
                    if residue.resid == 1:
                        atoms_to_remove.extend(residue.atoms[13:14].indices)  # Example indices for specific atoms
                        residue.atoms[9].name = 'N2'
                        residue.atoms[10].name = 'HN21'
                        residue.atoms[11].name = 'HN22'
                    else:
                        atoms_to_remove.extend(residue.atoms[11:12].indices)
                        residue.atoms[7].name = 'N2'
                        residue.atoms[8].name = 'HN21'
                        residue.atoms[9].name = 'HN22'
    
        return atoms_to_remove
    
    def nh2_apply_modifications_and_save(universe, atoms_to_remove, output_name):
        if not atoms_to_remove:
            #print("No atoms to remove. Saving the original universe.")
            universe.atoms.write(output_name)
            return
        query = "not bynum " + " ".join(map(str, atoms_to_remove))
        remaining_atoms = universe.select_atoms(query)
        new_universe = mda.Merge(remaining_atoms)
        new_universe.atoms.write(output_name)
        #print("Modification complete. New PDB file written:", output_name)
    
    
    nh3_deprotonation_chains = molecule.load("pdb", "nh3-modified.temp.pdb")
    left_nh2_candidates_indices   = get_fragment_indices(nh3_deprotonation_chains, sel_010_left_array)   
    upper_nh2_candidates_indices  = get_fragment_indices(nh3_deprotonation_chains, sel_100_upper_array)   
    

    
    #def get_resid(mol_id, sel):
    #    resid_indices = []
    #    for indices_i, indices_j in sel:
    #        selection_query = f"index {indices_i} to {indices_j}"
    #        selection = atomsel(selection_query, molid=mol_id)
    #        resid_num_indice = selection.resid
    #        unique_resid_set = set(resid_num_indice)
    #        if len(unique_resid_set) == 1:
    #            resid_indices.append(unique_resid_set.pop())  # Add the unique resid to the list
    #        else:
    #            raise ValueError(f"Multiple unique residues found in selection from index {indices_i} to {indices_j}")
    #    return resid_indices
    
    #def get_resid(mol_id, indices_i, indices_j):
    #    selection_query = f"index {indices_i} to {indices_j}"
    #    selection = atomsel(selection_query, molid=mol_id)
    #    resid_num_indice = selection.resid
    #    unique_resid = set(resid_num_indice)
    #    if len(unique_resid) == 1:
    #            return unique_resid.pop()  # Return the single unique fragment number
    #    else:
    #        raise ValueError("Multiple unique fragments found in selection")
    
    
    
    ##apply nh2 remove
    nh2_indices_groups = [left_nh2_candidates_indices, upper_nh2_candidates_indices]
    nh3_original_chains_u = mda.Universe("nh3-modified.temp.pdb")
    nh3_atoms_to_remove = collect_nh3_atoms_to_remove_randomly(nh3_original_chains_u, nh2_indices_groups, nh2_number)
    
    nh2_apply_modifications_and_save(nh3_original_chains_u, nh3_atoms_to_remove, "nh2-temp.pdb")
    
    
    #left_nh2_candidates_resid = get_resid(nh3_deprotonation_chains, left_nh2_candidates_indices) 
    #print(f"left resid is {left_nh2_candidates_resid}")
    molecule.delete(nh3_deprotonation_chains) 
    
    u_final = mda.Universe("nh2-temp.pdb")
    box_x = a_trans * a_iterations + 20
    box_y = b_iterations * b_trans + 20
    box_z = c_iterations * c_trans + 20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    u_final.atoms.translate(translation_vector)
    with mda.Writer(f"charmm36-alpha-chitin-A-fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
        
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    
    print(f"DDA: {actual_dda_rounded:.2f}, pH: {actual_pH_rounded:.4f}, Units: {nh2_number}")


elif aminiation_number == 0:
    u_final = mda.Universe("alpha-chitin-A-temp.pdb")
    box_x = a_trans * a_iterations + 20
    box_y = b_iterations * b_trans + 20
    box_z = c_iterations * c_trans + 20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    with mda.Writer(f"charmm36-alpha-chitin-A-fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
    

    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
        
    actual_dda_rounded=0.0
    actual_pH_rounded=0.0
    y_rounded=0.0
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"DDA: {actual_dda_rounded:.2f}, pH: {actual_pH_rounded:.4f}, Units: {y_rounded}")
