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

a_trans = 4.749 
b_trans = 18.89
c_trans = 10.33


try:
        a_iterations = int(sys.argv[1])
        b_iterations = int(sys.argv[2])
        c_iterations = int(sys.argv[3])
        
        if a_iterations <= 0 or b_iterations <= 0 or c_iterations <= 0:
            sys.stderr.write("Error: Please provide positive integers for all repetition numbers.\n")
            sys.exit(1)

        ##deacetylation degree ref:10.1021/acs.jchemed.7b00902J
        DDA_target = float(sys.argv[4])
        if not (0 <= DDA_target < 1):
            sys.stderr.write("Error: Deacetylation degree (DDA) must be between 0 and 1.\n") 
            sys.exit(1)
        #print(f"Validated input values:\n  a_iterations: {a_iterations}, b_iterations: {b_iterations}, c_iterations: {c_iterations}\n  DDA_target: {DDA_target}, pH: {pH}")

except IndexError:
    sys.stderr.write("Error: Missing one or more input arguments. Please provide all five effective values.\n")
    sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid input. Please ensure all values are numeric and within the valid ranges.\n")
    sys.exit(1)            # pH condition

#DDA_target=0.20               

pKa=6.3                        #pka of chitosan  ref:doi.org/10.1016/j.foodhyd.2022.108383
m_GlcNAc = 204.20              # GlcNAc unit mass
m_Glc = 179.17                 # Glc unit mass
#-------------------------------------------- modified parameter  --------------------------------------------------




#-------------------------------------------------- native chitin built ----------------------------------------------
with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
left_chain_input_file = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36',  'B_configuration', 'left-unit.pdb')
right_chain_input_file = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36', 'B_configuration', 'right-unit.pdb')
l_u = mda.Universe(left_chain_input_file)
r_u = mda.Universe(right_chain_input_file)


 #left_chain assembly    
 #------------------------------------------------------------------------------
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



#selected_atoms = chain_u.select_atoms(f"resid 1")
#print(f"Number of atoms in the selected residue range: {len(selected_atoms)}")

#print("Left Modified Chains:", left_modified_chains)
#print("Right Modified Chains:", right_modified_chains)


###NHx group evenly assigned to each chain
##define the dda; how many resiudes should be modified:
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
print(f"total residues number is {total_residues}")
#print(f"closet_x number is {closest_x}")

if DDA_target == 0:
    closest_x = total_residues
aminiation_number = closest_x

###aminiation number
#modified_num_per_chain = aminiation_number / (2 * a_iterations)
if  closest_x > 0:
    actual_dda=(aminiation_number*179.17)/(aminiation_number*179.17 + (total_residues-aminiation_number)*204.20)
    #print(f"Actual DDA is : {actual_dda}")
    #print(f"aminiation_number: {aminiation_number}")
    #print(f"aminiation_number_per_chain: {modified_num_per_chain}")
    actual_dda_rounded=round(actual_dda,5)
    #print(actual_dda_rounded)
    ##select the range

    #store the inside part of the structures
    selected_aminiation_residue_ranges = [] 
    all_possible_residue_numbers = []

    left_ranges = []
    right_ranges = []

    left_possible_residue_numbers = []  # For even-specific ranges
    right_possible_residue_numbers = []   # For odd-specific ranges


    for j in range(1, 2 * b_iterations-1):
        k=j+1
        start_residue = j*a_iterations*2*c_iterations + 2*c_iterations
        end_residue = k*a_iterations*2*c_iterations -1 - 2*c_iterations  ##c_iteration*2 representing the number of residue for one single chitin strand
    
        selected_aminiation_residue_ranges.append((start_residue, end_residue))
        #print(f"Start at {start_residue} and ends at {end_residue}")
        all_possible_residue_numbers.extend(range(start_residue, end_residue+1))

        if j % 2 == 0:  # left
            left_start_residue = start_residue  # Adjust this as necessary
            left_end_residue = end_residue  # Adjust this as necessary
            left_ranges.append((left_start_residue, left_end_residue))
            left_possible_residue_numbers.extend(range(left_start_residue, left_end_residue + 1))
            #print(left_possible_residue_numbers)
        else:  # right
            right_start_residue = start_residue  # Adjust this as necessary
            right_end_residue = end_residue  # Adjust this as necessary
            right_ranges.append((right_start_residue, right_end_residue))
            right_possible_residue_numbers.extend(range(right_start_residue, right_end_residue + 1))
            

    #print(f"All possible resiude number is {len(all_possible_residue_numbers)}")
    #print(f"left resiude number is {len(left_possible_residue_numbers)}")
    #print(f"right possible resiude number is {len(right_possible_residue_numbers)}")
    # Print the available ranges
    #print("Selected ranges:", selected_aminiation_residue_ranges)
    

    def get_residue_indices(mol_id, modified_residue_number):
        modified_residues_indices = []
        # Select atoms belonging to the specified fragment number
        residue_selection = atomsel(f"residue {modified_residue_number}", molid=mol_id)
        indices = residue_selection.index
        modified_residues_indices.extend(indices)
        return modified_residues_indices
  # 
    def n_acetyl_collect_removal_candidates(universe, start_idx, end_idx, modified_residues):
        n_acetyl_candidates_to_remove = []
        residue_selector = f"resid {modified_residues} and (bynum {start_idx}:{end_idx})"
        residues = universe.select_atoms(residue_selector).residues
        
        for residue in residues:
            if residue.resid == 1:
                n_acetyl_candidates_to_remove.extend(residue.atoms[10:17].indices)  
            else:
                n_acetyl_candidates_to_remove.extend(residue.atoms[8:15].indices) 
        return n_acetyl_candidates_to_remove    
    

    def apply_N_acetyl_removals(universe, candidates_to_remove):
        """Remove the collected candidate atoms from the universe and create a new modified universe."""
        all_indices = set(range(len(universe.atoms)))
        remaining_indices = all_indices - set(candidates_to_remove)
        remaining_atoms = universe.atoms[list(remaining_indices)]
        no_nacetyl = mda.Merge(remaining_atoms)
        return no_nacetyl
  # 

    # Selecting numbers randomly without repetition
    selected_numbers = []
    current_structure_name='alpha-chitin-A-temp.pdb'
    modified_loop_structure='alpha-chitin-A-temp-0.pdb'
    os.rename(current_structure_name, modified_loop_structure)
  # 
    try:#print(f"aminiation number is {aminiation_number}")
        for m in range(min(aminiation_number, len(all_possible_residue_numbers))):
            #print(f"Building structure step {m}")
            if all_possible_residue_numbers:  # Ensure the list is not empty
                random_residue_number = random.choice(all_possible_residue_numbers)
                
                origin_nacetyl_chains = molecule.load("pdb", modified_loop_structure)
                selected_modified_indices=get_residue_indices(origin_nacetyl_chains, random_residue_number)
                selected_modified_ndx_i = selected_modified_indices[0]
                selected_modified_ndx_j = selected_modified_indices[-1]
                molecule.delete(origin_nacetyl_chains)
                #selected_modified_atoms=chain_u.select_atoms(f"index {selected_modified_ndx_i} : {selected_modified_ndx_j}")
                #n-acetyl-removed steps
                original_chitin_crystal=mda.Universe(modified_loop_structure) 
                atom_i = original_chitin_crystal.atoms[selected_modified_ndx_i] 
                atom_j = original_chitin_crystal.atoms[selected_modified_ndx_j]
                modified_unit_reisd=atom_i.resid
                #print(modified_unit_reisd)
                selected_nacetyl_candidates = n_acetyl_collect_removal_candidates(original_chitin_crystal, selected_modified_ndx_i, selected_modified_ndx_j, modified_unit_reisd)  
                modified_universe = apply_N_acetyl_removals(original_chitin_crystal, selected_nacetyl_candidates)
                modified_universe.atoms.write(f"modified-alpha-chitin-temp-{m}.pdb")
    
    
                ##located whether this structure is left or right with nacetyls groups
                ##NH2 adding steps
                removal_nacetyl_resiudes = molecule.load("pdb", f"modified-alpha-chitin-temp-{m}.pdb")
                selected_removal_indices = get_residue_indices(removal_nacetyl_resiudes, random_residue_number)
                selected_removal_ndx_i = selected_removal_indices[0]
                selected_removal_ndx_j = selected_removal_indices[-1]
                molecule.delete(removal_nacetyl_resiudes)
                nacetyl_removal_residue_u = mda.Universe(f"modified-alpha-chitin-temp-{m}.pdb")
                nacetyl_removal_residue_atoms = list(nacetyl_removal_residue_u.atoms)
                nacetyl_removal_residue = nacetyl_removal_residue_u.select_atoms(f"index {selected_removal_ndx_i}:{selected_removal_ndx_j}").residues
                    # Process each selected residue
                for sel in nacetyl_removal_residue:
                    sel.resname = 'BLND' 
                    base_atom_index = 9 if sel .resid == 1 else 7
                    base_atom = sel .atoms[base_atom_index]
                    base_atom.name = 'N2'   
                    # Define new positions based on residue ID's odd/even status
                    if random_residue_number in left_possible_residue_numbers and sel .resid % 2 == 1:  # Odd residue ID
                        nh2_new_positions = {
                            'HN21': base_atom.position + np.array([ 0.942, -0.324, 0.092]),
                            'HN22': base_atom.position + np.array([-0.762, -0.622, 0.177]),
                        }
                    elif random_residue_number in left_possible_residue_numbers and sel .resid % 2 == 0:   
                        nh2_new_positions = {
                            'HN21': base_atom.position + np.array([ 0.762,  0.622,  0.177]),
                            'HN22': base_atom.position + np.array([-0.942,  0.324,  0.092]),
                        }
                    if random_residue_number in right_possible_residue_numbers and sel .resid % 2 == 1:  # Odd residue ID
                        nh2_new_positions = {
                            'HN21': base_atom.position + np.array([ 0.762, -0.622, -0.177]),
                            'HN22': base_atom.position + np.array([-0.942, -0.324, -0.092]),
                        }
                    elif random_residue_number in right_possible_residue_numbers and sel .resid % 2 == 0:   
                        nh2_new_positions = {
                            'HN21': base_atom.position + np.array([ 0.942,  0.324, -0.092]),
                            'HN22': base_atom.position + np.array([-0.762,  0.622, -0.177]),
                        }
                    # Insert new atoms
                    insert_pos = nacetyl_removal_residue_atoms.index(base_atom) + 1
                    for name, pos in nh2_new_positions.items():
                        nh2_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
                        nh2_new_uni.add_TopologyAttr('name', [name])
                        nh2_new_uni.add_TopologyAttr('type', [base_atom.type])
                        nh2_new_uni.add_TopologyAttr('resname', [base_atom.resname])
                        nh2_new_uni.add_TopologyAttr('resid', [sel.resid])
                        nh2_new_uni.add_TopologyAttr('segid', [base_atom.segment.segid])
                        nh2_new_uni.add_TopologyAttr('chainIDs', ['N'])
                        nh2_new_uni.atoms.positions = [pos]
                        left_nh3_new_atom = nh2_new_uni.atoms[0]
            
                        # Insert the new atom at the calculated position
                        nacetyl_removal_residue_atoms.insert(insert_pos, left_nh3_new_atom)
                        insert_pos += 1  # Update position for the next atom
                
                # Create a new universe with all atoms including the added ones
                nh2_new_universe = mda.Merge(*[mda.AtomGroup([atom]) for atom in nacetyl_removal_residue_atoms])
    
                new_round=m+1 #for the loop update and also for the last round , which is for the final structure
                nh2_new_universe.atoms.write(f"alpha-chitin-A-temp-{new_round}.pdb")
                modified_loop_structure=f'alpha-chitin-A-temp-{new_round}.pdb'
    
                ###avoid repeated to modify the same residue
                all_possible_residue_numbers.remove(random_residue_number)  # Remove the selected number to prevent re-selection
                selected_numbers.append(random_residue_number)
            else:
                print("No valid numbers to choose from or all numbers have been selected.")
            #print(f"Randomly selected number: {random_residue_number}")
    except KeyboardInterrupt:
        print("Loop interrupted by user (Ctrl+C). Exiting...")


    ###final round number is m
    final_loop_number=new_round
    u_final = mda.Universe(f"alpha-chitin-A-temp-{new_round}.pdb")
    box_x = a_iterations * a_trans + 20
    box_y = b_iterations * b_trans + 20
    box_z = c_iterations * c_trans + 20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    u_final.atoms.translate(translation_vector)
    with mda.Writer(f"charmm36-alpha-chitin-B-fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)    
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"DDA: {actual_dda_rounded:.4f}, Units: {min(aminiation_number, len(all_possible_residue_numbers))}")

    #selected_numbers_length = len(selected_numbers)
    #print("Total selected numbers:", selected_numbers_length)
    #print("All selected numbers:", selected_numbers)


elif aminiation_number == 0:
    u_final = mda.Universe("alpha-chitin-A-temp.pdb")
    box_x = a_iterations * a_trans + 20
    box_y = b_iterations * b_trans + 20
    box_z = c_iterations * c_trans + 20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    u_final.atoms.translate(translation_vector)
    with mda.Writer(f"charmm36-alpha-chitin-B-fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
    
    #complete_structure=mda.Universe("u_final-temp.pdb")
    #for atom in complete_structure.atoms:
    #    elements = [atom.name[0] for atom in complete_structure.atoms]
    #    complete_structure.add_TopologyAttr(Elements(elements)) 
    #complete_structure.atoms.write("charmm36-alpha-chitin-A-fcm.pdb")
    
    # Remove temporary files
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    actual_dda_rounded=0.0
    aminiation_number=0
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"DDA: {actual_dda_rounded:.2f}, Units: {aminiation_number}")
