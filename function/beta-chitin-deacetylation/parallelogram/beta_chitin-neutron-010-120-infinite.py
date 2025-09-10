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



a_trans     = 4.819 
b_trans     = 9.238976
c_trans     = 10.384001
gamma_angle = 97.156586

try:
    c_iterations = int(sys.argv[1])
    if c_iterations <= 0:
        sys.stderr.write("Error: Please provide positive integers for c repetition numbers.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length value provided. Please enter valid numeric values.\n")
    sys.exit(1)
except IndexError:
    sys.stderr.write("Error: Missing one or more input arguments. Please provide three integer values.\n")
    sys.exit(1)

try:
    height = float(sys.argv[2])
    a_iterations = int(height // a_trans)
    if a_iterations < 1:
        sys.stderr.write("Error: Please provide an larger height value.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid width value provide. Please enter a valid numeric value.\n")
    sys.exit(1)    

try:
    width = float(sys.argv[3])
    unit_A = b_trans * math.cos((gamma_angle - 90) * math.pi / 180)
    unit_B = b_trans * math.sin((gamma_angle - 90) * math.pi / 180) + a_trans
    unit_width=math.sqrt( unit_A **2 + unit_B**2 )
    b_iterations = int(width // unit_width)
    if b_iterations < 1:
        sys.stderr.write("Error: Please provide an larger width value.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid width value provide. Please enter a valid numeric value.\n")
    sys.exit(1)    


try:
        ##deacetylation degree ref:10.1021/acs.jchemed.7b00902J
        DDA_target = float(sys.argv[4])
        if not (0 <= DDA_target < 1):
            sys.stderr.write("Error: Deacetylation degree (DDA) must be between 0 and 1.\n") 
            sys.exit(1)

        pH = float(sys.argv[5])
        if not (0 <= pH <= 14):
            sys.stderr.write("Error: pH must be between 0 and 14.\n")
            sys.exit(1)
        #print(f"Validated input values:\n  a_iterations: {a_iterations}, b_iterations: {b_iterations}, c_iterations: {c_iterations}\n  DDA_target: {DDA_target}, pH: {pH}")

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
unit_chain_input_file = os.path.join(main_folder_path, 'structure', 'beta_chitin', 'charmm36', 'unit-neutron.pdb')

unit_u = mda.Universe(unit_chain_input_file)

pKa=6.3                        #pka of chitosan  ref:doi.org/10.1016/j.foodhyd.2022.108383
m_GlcNAc = 204.20              # GlcNAc unit mass
m_Glc = 179.17                 # Glc unit mass



 #chain assembly    
 #------------------------------------------------------------------------------
unit_chain = []
for unit_i in range(1, c_iterations + 1):
    unit_u.atoms.positions -= [0, 0, c_trans]
    resid_1 = unit_i  * 2 - 1
    resid_2 = unit_i  * 2
    for j, atom in enumerate(unit_u.atoms):
        if j < 27:
            atom.residue.resid = resid_1
        else:
            atom.residue.resid = resid_2
    unit_chain_output = f"unit-{unit_i}.pdb"
    unit_chain.append(unit_chain_output)
    with mda.Writer(unit_chain_output, n_atoms=unit_u.atoms.n_atoms) as W:
        W.write(unit_u.atoms)
with open("unit-chain_temp.pdb", "w") as chain:
    for unit_chain_output in unit_chain:
         with open(unit_chain_output, "r") as unit_chain_pdb_file:
            for unit_chain_line in unit_chain_pdb_file:
                if unit_chain_line.startswith("ATOM"):
                    chain.write(unit_chain_line)
            os.remove(unit_chain_output)


unit_start_segid=a_iterations
unit_combined = "unit-chain_temp.pdb"
unit_combined_temp = mda.Universe(unit_combined)
for atom in unit_combined_temp.atoms:
    atom.segment.segid = '0'  
    elements = [atom.name[0] for atom in unit_combined_temp.atoms]
    unit_combined_temp.add_TopologyAttr(Elements(elements))
unit_combined_temp.atoms.write("unit-chain_temp.2.pdb")
#------------------------------------------------------------------------------


 #layer assembly 
unit_layer_input_file = "unit-chain_temp.2.pdb"
layer_unit_u = mda.Universe(unit_layer_input_file)
unit_layer = []  
for i in range(1, a_iterations + 1):
    layer_unit_u.atoms.positions += [a_trans, 0, 0]

    segid_increment_value = 1 
    for atom in layer_unit_u.atoms:
        numeric_part = atom.segid.replace('','')
        new_numeric_part = int(numeric_part) + segid_increment_value
    new_segid = f"{new_numeric_part}"
    atom.segment.segid = new_segid

    unit_layer_output = f"unit_layer_{i}.pdb"
    unit_layer.append(unit_layer_output)
    with mda.Writer(unit_layer_output, n_atoms=layer_unit_u.atoms.n_atoms) as W:
        W.write(layer_unit_u.atoms)
with open("unit-layer_temp.pdb", "w") as layer_file:
    for unit_layer_output in unit_layer:
        with open(unit_layer_output, "r") as unit_layer_pdb_file:
            for unit_layer_line in unit_layer_pdb_file:
                if unit_layer_line.startswith("ATOM"):
                    layer_file.write(unit_layer_line)
            os.remove(unit_layer_output) 
 


unit_output_file = "unit.pdb" 
with open(unit_output_file, "w") as unit_file:
    with open("unit-layer_temp.pdb", "r") as assembled_file:
        for line in assembled_file:
            if line.startswith("ATOM"):
                unit_file.write(line)


unit_input_file="unit.pdb"
crystal_structure = mda.Universe(unit_input_file)
crystal_structure_files = [] 


angle_rad = (gamma_angle - 90) * math.pi / 180
cosA = math.cos(angle_rad)
sinA = math.sin(angle_rad)

for crystal_i in range(1, b_iterations + 1):
    crystal_structure.atoms.positions += [-b_trans * sinA +2*a_trans, b_trans * cosA , 0]
    crystal_structure_output = f"crystal_{crystal_i}.pdb"
    crystal_structure_files.append(crystal_structure_output)
    with mda.Writer(crystal_structure_output, n_atoms=crystal_structure.atoms.n_atoms) as writer:
        writer.write(crystal_structure.atoms)


for file_index, file_name in enumerate(crystal_structure_files):
    temp_structure = mda.Universe(file_name)
    segid_increment_value = file_index  * a_iterations
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
with open("unit_temp.pdb", "w") as combined_file:
    for file_name in crystal_structure_files:
        with open(file_name, "r") as file:
            for line in file:
                if line.startswith("ATOM"):
                    combined_file.write(line)
        os.remove(file_name) 

os.remove(unit_input_file) 


##crystallographic planes defined
planes_010_left = []
planes_010_right = []


planes_120_upper = []
planes_120_bottom = []

##-------------------crystallographic-plane-definitions---------------------
for m in range(1, b_iterations-1):
    value1 = m*a_iterations
    value2 = (m+1)*a_iterations - 1
    planes_120_bottom.extend([value1])
    planes_120_upper.extend([value2]) 


range1_start = 0
range1_end = a_iterations-1
range2_start = b_iterations*a_iterations-a_iterations
range2_end = b_iterations*a_iterations-1
planes_010_right.extend([range1_start, range1_end])
planes_010_left.extend([range2_start, range2_end])




planes_010_count=(planes_010_left[-1] - planes_010_left[0] + 1) + (planes_010_right[-1] - planes_010_right[0] + 1)
planes_120_count=len(planes_120_bottom)+len(planes_120_upper)
modified_chains_count=planes_010_count+planes_120_count

#print(planes_010_count)
#print(planes_120_count)
#print("120_planes_upper:",  planes_120_upper)
#print("120_planes_bottom:", planes_120_bottom)
#
#print("010_planes_r:", planes_010_right)
#print("010_planes_l:", planes_010_left)
#




##------------------------------------------------- nh3 modified step -------------------------------------------------------
#
#
############NH3 located step

chain_u = mda.Universe("unit_temp.pdb")
residue_ids = [residue.resid for residue in chain_u.residues]
total_residues=len(residue_ids)
#
def calculate_dda(x):
    numerator = (total_residues - x) * m_Glc
    denominator = numerator + x * m_GlcNAc
    return numerator / denominator if denominator != 0 else 0


def adjust_for_amination_limit(closest_x, total_residues):
    amination_number = total_residues - closest_x
    amination_number_limit=c_iterations*planes_010_count + 2*c_iterations*planes_120_count
    if amination_number > amination_number_limit:
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




if aminiation_number > 0:
      
    left_010_count=0
    upper_120_count=0
    for i in range(1, aminiation_number + 1):
        if i % 2 == 0:
            upper_120_count += 1
        elif i % 2 == 1:
            left_010_count += 1

    
    new_aminiation_number=aminiation_number

    planes_010_left_count=(planes_010_left[-1] - planes_010_left[0] + 1)
    planes_010_right_count=(planes_010_right[-1] - planes_010_right[0] + 1)
    planes_120_upper_count=len(planes_120_upper)
    planes_120_bottom_count=len(planes_120_bottom)
    
    max_010_count_left = c_iterations * planes_010_left_count 
    max_010_count_right = c_iterations * planes_010_right_count
    
    try:
        max_120_count_upper = c_iterations * 2 * planes_120_upper_count
        max_120_count_bottom = c_iterations * 2 * planes_120_bottom_count
        if max_120_count_upper == 0 or  max_120_count_bottom == 0:
            sys.stderr.write("Error: Please provide a larger width value.\n")
            sys.exit(1)
    except ValueError:
        sys.stderr.write("Error: Invalid width value provide. Please enter a larger width value.\n")
        sys.exit(1)    
    
    
    left_010_count_update = left_010_count
    
    upper_120_count_update = upper_120_count
    if left_010_count > max_010_count_left:
        left_010_count_update = max_010_count_left
    
    if upper_120_count > max_120_count_upper:
        upper_120_count_update = max_120_count_upper
        
    new_aminiation_number_update=upper_120_count_update +  left_010_count_update
    actual_dda=(new_aminiation_number_update*179.17)/(new_aminiation_number_update*179.17 + (total_residues-new_aminiation_number_update)*204.20)
    actual_dda_rounded=round(actual_dda,4)

    even_numbers_010 = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 0]
    odd_numbers_010 = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 1]
    all_combinations_010_left = [(value1, value2) for value1 in range(min(planes_010_left), max(planes_010_left) + 1)  for value2 in odd_numbers_010 ]
    all_combinations_010_right = [(value1, value2) for value1 in range(min(planes_010_right), max(planes_010_right) + 1)  for value2 in even_numbers_010 ]
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
    
    
   
    numbers_120 = [x for x in range(1, 2 * c_iterations + 1)]
    
    all_combinations_120_upper = [(value1, value2) for value1 in planes_120_upper  for value2 in numbers_120 ]
    all_combinations_120_bottom = [(value1, value2) for value1 in planes_120_bottom  for value2 in numbers_120 ]
    random.shuffle(all_combinations_120_upper)
    random.shuffle(all_combinations_120_bottom)
    sel_120_upper_array = []
    sel_120_bottom_array = []
        
    while len(sel_120_upper_array) < upper_120_count_update:
        if not all_combinations_120_upper:  
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_120_upper = random.choice(all_combinations_120_upper)
        sel_120_upper_array.append(selected_combination_120_upper)
        all_combinations_120_upper.remove(selected_combination_120_upper)


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
    
    origin_nacetyl_chains = molecule.load("pdb", "unit_temp.pdb")
    left_nacetyl_candidates_indices   = get_fragment_indices(origin_nacetyl_chains, sel_010_left_array)   
    upper_nacetyl_candidates_indices  = get_fragment_indices(origin_nacetyl_chains, sel_120_upper_array)    
    molecule.delete(origin_nacetyl_chains) 



    def remove_atoms_based_on_indices(universe, index_ranges):
        all_removed_indices = []
    
        for index_range in index_ranges:
            for start_idx, end_idx in index_range:
                #print(f"Handling indices from {start_idx} to {end_idx}")
                atom_group = universe.atoms[start_idx:end_idx + 1]
                for residue in atom_group.residues:
                    indices_to_remove = residue.atoms[8:15].indices
                    all_removed_indices.extend(indices_to_remove)
        unique_indices = sorted(set(all_removed_indices))
        indices_str = ' '.join(map(str, unique_indices))
        universe = universe.select_atoms(f"not index {indices_str}")
    
        return universe


    
    all_indices = [left_nacetyl_candidates_indices, upper_nacetyl_candidates_indices]
    ###chain_u is the unmodified structure, with name of "unit_temp.pdb"
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
    upper_nh3_candidates_indices  = get_fragment_indices(nacetyl_remove_chains, sel_120_upper_array)   


    nacetyl_remove_chains_u = mda.Universe("nacetyl_remove.temp.pdb")
    nh3_atoms = list(nacetyl_remove_chains_u.atoms)



    
    #print(left_nh3_candidates_indices)
    def add_nh3_based_on_indices(universe, indices_lists):
        all_new_atoms = [] 
        for index_ranges in indices_lists:
            for start_idx, end_idx in index_ranges:
                residues = universe.select_atoms(f"index {start_idx}:{end_idx}").residues
                fragment_number = get_fragment(nacetyl_remove_chains, start_idx, end_idx)
                for residue in residues:
                    #print("residue is : ", residue)
                    residue.resname = 'BLNP'
                    base_atom_index =  7
                    base_atom = residue.atoms[base_atom_index]
                    nh3_new_positions = {} 
                    # Define new positions based on residue ID's odd/even status
                    if residue.resid % 2 == 0  and fragment_number in planes_120_bottom:  # Odd residue ID
                        nh3_new_positions = {
                            'HN1': base_atom.position + np.array([ 0.950,   -0.211,   0.231]),
                            'HN2': base_atom.position + np.array([-0.581,   -0.159,   0.798]),
                            'HN3': base_atom.position + np.array([-0.297,   -0.588,  -0.752])
                        }
                    elif residue.resid % 2 == 1  and fragment_number in planes_120_bottom:  # Even residue ID
                        nh3_new_positions = {
                            'HN1': base_atom.position + np.array([-0.950,  0.211,  0.231]),
                            'HN2': base_atom.position + np.array([ 0.581,  0.159,  0.798]),
                            'HN3': base_atom.position + np.array([ 0.298,  0.588, -0.752])
                        }
    
                    elif residue.resid % 2 == 0  and fragment_number in planes_120_upper:  # Odd residue ID
                        nh3_new_positions = {
                            'HN1': base_atom.position + np.array([ 0.950,   -0.211,   0.231]),
                            'HN2': base_atom.position + np.array([-0.581,   -0.159,   0.798]),
                            'HN3': base_atom.position + np.array([-0.297,   -0.588,  -0.752])
                        }
                    elif residue.resid % 2 == 1  and fragment_number in planes_120_upper:  # Even residue ID
                        nh3_new_positions = {
                            'HN1': base_atom.position + np.array([-0.950,  0.211,  0.231]),
                            'HN2': base_atom.position + np.array([ 0.581,  0.159,  0.798]),
                            'HN3': base_atom.position + np.array([ 0.298,  0.588, -0.752])
                        }
    
                    elif fragment_number in range(min(planes_010_left), max(planes_010_left)+1):  # Odd residue ID
                        nh3_new_positions = {
                            'HN1': base_atom.position + np.array([-0.950,    0.211,   0.231]),
                            'HN2': base_atom.position + np.array([ 0.581,    0.159,   0.798]),
                            'HN3': base_atom.position + np.array([ 0.298,    0.588,  -0.752])
                        }
                    elif fragment_number in range(min(planes_010_right), max(planes_010_right)+1):  # Even residue ID
                        nh3_new_positions = {
                            'HN1': base_atom.position + np.array([ 0.950, -0.211,  0.231]),
                            'HN2': base_atom.position + np.array([-0.581, -0.159,  0.798]),
                            'HN3': base_atom.position + np.array([-0.298, -0.588, -0.752])
                        }

                    insert_pos = nh3_atoms.index(base_atom) + 1
                    for name, pos in nh3_new_positions.items():
                        nh3_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
                        nh3_new_uni.add_TopologyAttr('name', [name])
                        nh3_new_uni.add_TopologyAttr('type', [base_atom.type])
                        nh3_new_uni.add_TopologyAttr('resname', [base_atom.resname])
                        nh3_new_uni.add_TopologyAttr('resid', [residue.resid])
                        nh3_new_uni.add_TopologyAttr('segid', [base_atom.segment.segid])
                        nh3_new_uni.add_TopologyAttr('chainIDs', ['X'])
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
    upper_nh2_candidates_indices  = get_fragment_indices(nh3_deprotonation_chains, sel_120_upper_array)   

    
    ##apply nh2 remove
    nh2_indices_groups = [left_nh2_candidates_indices,  upper_nh2_candidates_indices]
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
    with mda.Writer(f"charmm36-beta-chitin-para-icm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
        
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    
    print(f"DDA: {actual_dda_rounded:.2f}, pH: {actual_pH_rounded:.4f}, Units: {nh2_number}")


    
elif aminiation_number == 0:
    u_final = mda.Universe("unit_temp.pdb")
    box_x = a_trans * a_iterations + 20
    box_y = b_iterations * b_trans + 20
    box_z = c_iterations * c_trans + 20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    with mda.Writer(f"charmm36-beta-chitin-para-icm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
    

    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
        
    actual_dda_rounded=0.0
    actual_pH_rounded=0.0
    y_rounded=0.0
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"DDA: {actual_dda_rounded:.2f}, pH: {actual_pH_rounded:.4f}, Units: {y_rounded}")

