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




#-------------------------------------------------- native chitin built ----------------------------------------------
with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
left_chain_input_file = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36', 'A_configuration', 'left-unit.pdb')
right_chain_input_file = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36', 'A_configuration', 'right-unit.pdb')
l_u = mda.Universe(left_chain_input_file)
r_u = mda.Universe(right_chain_input_file)


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

left_temp = "left-chain_temp.pdb"
left_temp_u=mda.Universe(left_temp)

for atom in left_temp_u.atoms:
    atom.segment.segid = '0'  
left_temp_u.atoms.write("left-chain_temp.2.pdb")
#------------------------------------------------------------------------------




 #right_chain assembly    
 # ------------------------------------------------   
right_chain = []
right_side_start_segid=a_iterations
for r_i in range(1, c_iterations + 1):
    r_u.atoms.positions += [0, 0, c_trans]
    resid_1 = r_i * 2 - 1
    resid_2 = r_i * 2
    for j, atom in enumerate(r_u.atoms):
        if j < 27:
            atom.residue.resid = resid_1
        else:
            atom.residue.resid = resid_2
    right_chain_output = f"r-{r_i}.pdb"
    right_chain.append(right_chain_output)
    with mda.Writer(right_chain_output, n_atoms=r_u.atoms.n_atoms) as W:
        W.write(r_u.atoms)
with open("right-chain_temp.pdb", "w") as chain_file:
    for r_chain_output in right_chain:
        with open(r_chain_output, "r") as right_chain_pdb_file:
            for line in right_chain_pdb_file:
                if line.startswith("ATOM"):
                    chain_file.write(line)
        os.remove(r_chain_output)


right_side_start_segid=a_iterations
right_temp = "right-chain_temp.pdb"

right_temp_u=mda.Universe(right_temp)

for atom in right_temp_u.atoms:
    atom.segment.segid = f'{right_side_start_segid}'  
right_temp_u.atoms.write("right-chain_temp.2.pdb")

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




#------------------------------------------------- nh2 modified step -------------------------------------------------------


###########NH2 located step

chain_u = mda.Universe("alpha-chitin-A-temp.pdb")
residue_ids = [residue.resid for residue in chain_u.residues]
#resid_ids = [residue.resid for residue in chain_u.resid]
total_residues=len(residue_ids)
origin_nacetyl_chains = molecule.load("pdb", "alpha-chitin-A-temp.pdb")
native_all_atoms = atomsel('all', molid=origin_nacetyl_chains)
residue_ids = native_all_atoms.residue
#min_resid = min(residue_ids, default=None)
#max_resid = max(residue_ids, default=None)
#print(f"Minimum residue ID: {min_resid}")
#print(f"Maximum residue ID: {max_resid}")
molecule.delete(origin_nacetyl_chains)

#residue_selection = atomsel(f"residue {residu_number}", molid=origin_nacetyl_chains)



def calculate_dda(x):
    numerator = (x) * m_Glc
    denominator = numerator + (total_residues - x) * m_GlcNAc
    return numerator / denominator if denominator != 0 else 0

def adjust_for_amination_limit(closest_x):
    closest_x
    if closest_x > 2 * c_iterations * ((a_iterations-2) * (2*b_iterations-2)):
        new_amination_number = 2 * c_iterations * ((a_iterations-2) * (2*b_iterations-2)) ######limitation to the all possible selected resiudes for inside deacetylation
        closest_x = new_amination_number
    return closest_x ##Glc number is closes_x

closest_x = None
closest_dda = float('inf')
for x in range(total_residues + 1):
    dda = calculate_dda(x)
    if abs(dda - DDA_target) < abs(closest_dda - DDA_target):
        closest_x = x
        closest_dda = dda
#print(f"closet_x number is {closest_x}")

closest_x = adjust_for_amination_limit(closest_x)
#print(f"total residues number is {total_residues}")
#print(f"closet_x number is {closest_x}")

if DDA_target == 0:
    closest_x = total_residues
aminiation_number = closest_x

left_fragment, right_fragment = [], []

for j in range(1, 2 * b_iterations - 1):
    start_fragment = j * a_iterations + 1
    end_fragment   = start_fragment + a_iterations - 2   # inclusive end

    if end_fragment >= start_fragment:  # guard against empty ranges
        if j % 2 == 0:   # left
            left_fragment.extend(range(start_fragment, end_fragment))   # include end
        else:            # right
            right_fragment.extend(range(start_fragment, end_fragment))  # include end

###aminiation number
#modified_num_per_chain = aminiation_number / (2 * a_iterations)
if  aminiation_number > 0:
     
    #print(f"Actual DDA is : {actual_dda}")
    #print(f"aminiation_number: {aminiation_number}")
    #print(f"aminiation_number_per_chain: {modified_num_per_chain}")
    
    #print(actual_dda_rounded)
    ##select the range

    left_residue_count=0
    right_residue_count=0

    for i in range(1, aminiation_number + 1):
        if i % 2 == 0:
            left_residue_count += 1
        elif i % 2 == 1:
            right_residue_count += 1


    left_fragment_count=len(left_fragment)
    right_fragment_count=len(right_fragment)
    max_left_residue_count = 2 * c_iterations * left_fragment_count
    max_right_residue_count = 2*  c_iterations * right_fragment_count

    left_residue_count_update=left_residue_count
    right_residue_count_update=right_residue_count

    if left_residue_count > max_left_residue_count:
       left_residue_count_update = max_left_residue_count    

    if right_residue_count > max_right_residue_count:
       right_residue_count_update = max_right_residue_count

    new_aminiation_number_update=left_residue_count_update + right_residue_count_update 
    actual_dda=(new_aminiation_number_update*179.17)/(new_aminiation_number_update*179.17 + (total_residues-new_aminiation_number_update)*204.20)
    actual_dda_rounded=round(actual_dda,5)

    numbers = [x for x in range(1, 2 * c_iterations + 1)]
    
    all_combinations_left = [(value1, value2) for value1 in left_fragment  for value2 in numbers ]
    all_combinations_right = [(value1, value2) for value1 in right_fragment  for value2 in numbers ]

    
    #print(len(all_combinations_010_left))
    random.shuffle(all_combinations_left)
    random.shuffle(all_combinations_right)
    sel_left_array = []
    sel_right_array = []
    
    while len(sel_left_array) < left_residue_count_update:
        if not all_combinations_left:  
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_left = random.choice(all_combinations_left)
        sel_left_array.append(selected_combination_left)
        all_combinations_left.remove(selected_combination_left)

    while len(sel_right_array) < right_residue_count_update:
        if not all_combinations_right:  # Check if the list is empty to avoid IndexError
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_right = random.choice(all_combinations_right)
        sel_right_array.append(selected_combination_right)
        all_combinations_right.remove(selected_combination_right)




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
    left_nacetyl_candidates_indices   = get_fragment_indices(origin_nacetyl_chains, sel_left_array)   
    right_nacetyl_candidates_indices  = get_fragment_indices(origin_nacetyl_chains, sel_right_array)   
    molecule.delete(origin_nacetyl_chains) 

    def remove_atoms_based_on_indices(universe, index_ranges):
        all_removed_indices = []
    
        for index_range in index_ranges:
            for start_idx, end_idx in index_range:
                atom_group = universe.atoms[start_idx:end_idx + 1]
                for residue in atom_group.residues:
                    indices_to_remove = residue.atoms[8:15].indices
                    all_removed_indices.extend(indices_to_remove)
        unique_indices = sorted(set(all_removed_indices))
        indices_str = ' '.join(map(str, unique_indices))
        universe = universe.select_atoms(f"not index {indices_str}")
    
        return universe
    
    
    
    all_indices = [left_nacetyl_candidates_indices, right_nacetyl_candidates_indices]
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
    left_nh2_candidates_indices   = get_fragment_indices(nacetyl_remove_chains, sel_left_array)   
    right_nh2_candidates_indices  = get_fragment_indices(nacetyl_remove_chains, sel_right_array)   


    nacetyl_remove_chains_u = mda.Universe("nacetyl_remove.temp.pdb")
    nh2_atoms = list(nacetyl_remove_chains_u.atoms)


    def add_nh2_based_on_indices(universe, indices_lists):
        all_new_atoms = []  # This will store all newly created atoms across all modifications
        for index_ranges in indices_lists:
            for start_idx, end_idx in index_ranges:
                residues = universe.select_atoms(f"index {start_idx}:{end_idx}").residues
                fragment_number = get_fragment(nacetyl_remove_chains, start_idx, end_idx)

                #print(residues_num_1)
                #print("staring from to the end:",start_idx, end_idx)
                #print("fragment number is :", fragment_number )
                #print("residue is : ", residues)
                for residue in residues:
                    residue.resname = 'BLND'
                    base_atom_index =  7
                    base_atom = residue.atoms[base_atom_index]
                    nh2_new_positions = {} 
                    #print(residues_num)
                    # Define new positions based on residue ID's odd/even status
                    if residue.resid % 2 == 1  and fragment_number in left_fragment: # Odd residue ID
                        nh2_new_positions = {
                            'HN21': base_atom.position + np.array([ 0.942, -0.324, 0.092]),
                            'HN22': base_atom.position + np.array([-0.762, -0.622, 0.177]),
                        }
                    elif residue.resid % 2 == 0  and fragment_number in left_fragment:  # Odd residue ID
                        nh2_new_positions = {
                            'HN21': base_atom.position + np.array([ 0.762,  0.622,  0.177]),
                            'HN22': base_atom.position + np.array([-0.942,  0.324,  0.092]),
                        }
    
                    elif residue.resid % 2 == 1  and fragment_number in right_fragment: # Odd residue ID
                        nh2_new_positions = {
                            'HN21': base_atom.position + np.array([ 0.762, -0.622, -0.177]),
                            'HN22': base_atom.position + np.array([-0.942, -0.324, -0.092]),
                        }
                    elif residue.resid % 2 == 0  and fragment_number in right_fragment:  # Odd residue ID
                        nh2_new_positions = {
                            'HN21': base_atom.position + np.array([ 0.942,  0.324, -0.092]),
                            'HN22': base_atom.position + np.array([-0.762,  0.622, -0.177]),
                        }
                    residue.atoms[7].name = 'N2'
                    insert_pos = nh2_atoms.index(base_atom) + 1
                    for name, pos in nh2_new_positions.items():
                        nh2_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
                        nh2_new_uni.add_TopologyAttr('name', [name])
                        nh2_new_uni.add_TopologyAttr('type', [base_atom.type])
                        nh2_new_uni.add_TopologyAttr('resname', [base_atom.resname])
                        nh2_new_uni.add_TopologyAttr('resid', [residue.resid])
                        nh2_new_uni.add_TopologyAttr('segid', [base_atom.segment.segid])
                        nh2_new_uni.add_TopologyAttr('chainIDs', ['N'])
                        nh2_new_uni.atoms.positions = [pos]
                        nh2_new_atom = nh2_new_uni.atoms[0]
    
                        nh2_atoms.insert(insert_pos, nh2_new_atom)
                        insert_pos += 1  # Update position for the next atom
        # Create a new universe with all atoms including the added ones
        new_universe = mda.Merge(*[mda.AtomGroup([atom]) for atom in nh2_atoms])
        return new_universe
    
    
    all_indices = [left_nh2_candidates_indices, right_nh2_candidates_indices]
    nh2_u = add_nh2_based_on_indices(nacetyl_remove_chains_u, all_indices)
    nh2_u.atoms.write("nh2-modified.temp.pdb")
    molecule.delete(nacetyl_remove_chains)



    u_final = mda.Universe(f"nh2-modified.temp.pdb")
    box_x = a_iterations * a_trans
    box_y = b_iterations * b_trans
    box_z = c_iterations * c_trans
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    u_final.atoms.translate(translation_vector)
    with mda.Writer(f"charmm36-alpha-chitin-A-icm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)    
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"DDA: {actual_dda_rounded:.4f}, Units: {new_aminiation_number_update}")

    #selected_numbers_length = len(selected_numbers)
    #print("Total selected numbers:", selected_numbers_length)
    #print("All selected numbers:", selected_numbers)


elif aminiation_number == 0:
    u_final = mda.Universe("alpha-chitin-A-temp.pdb")
    box_x = a_iterations * a_trans
    box_y = b_iterations * b_trans
    box_z = c_iterations * c_trans
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    u_final.atoms.translate(translation_vector)
    with mda.Writer(f"charmm36-alpha-chitin-A-icm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
    
    # Remove temporary files
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    actual_dda_rounded=0.0
    aminiation_number=0
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"DDA: {actual_dda_rounded:.2f}, Units: {aminiation_number}")
