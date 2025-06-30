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


        ##carboxylation degree ref:/10.1016/j.carbpol.2019.115292
        Tempo_target = float(sys.argv[7])  #####
        #if not (0 <= Tempo_target < 1):
        #    sys.stderr.write("Error: Tempo-oxidation degree  must be between 0 and 1.\n") 
        #    sys.exit(1)

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
    sys.exit(1)            # pH condition

              
pKa=3.25                             #pka of coo-
m_GlcNAc  = 204.20                   # GlcNAc unit mass
m_difference = 36.973769               # difference between GlcNAc unit with GlcNAc-cooNa
#-------------------------------------------- modified parameter  --------------------------------------------------




with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
neutron_pdb_1 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36', 'AB_configuration', 'left-unit.pdb')
user_defined_pdb_1 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36', 'AB_configuration', 'left-unit_ud.pdb')

if os.path.exists(user_defined_pdb_1):
    unit_chain_input_file_1 = user_defined_pdb_1
else:
    unit_chain_input_file_1 = neutron_pdb_1

neutron_pdb_2 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36', 'AB_configuration', 'right-unit.pdb')
user_defined_pdb_2 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36', 'AB_configuration', 'right-unit_ud.pdb')

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

##------------------------------------------------- carboxylationsd  modified step -------------------------------------------------------

chain_u = mda.Universe("alpha-chitin-A-temp.pdb")
residue_ids = [residue.resid for residue in chain_u.residues]
#print(residue_ids[-1])
maximum_resid_num=residue_ids[-1]
atom_ids = [atom.id for atom in chain_u.atoms]
total_residues=len(residue_ids)
total_atom=len(atom_ids)
#c_atoms = chain_u.select_atoms("type C")
#h_atoms = chain_u.select_atoms("type H")
#o_atoms = chain_u.select_atoms("type O")
#n_atoms = chain_u.select_atoms("type N")

#c_atom_ids = [atom.id for atom in c_atoms]
#h_atom_ids = [atom.id for atom in h_atoms]
#o_atom_ids = [atom.id for atom in o_atoms]
#n_atom_ids = [atom.id for atom in n_atoms]
#n_c=len(c_atom_ids)
#n_h=len(h_atom_ids)
#n_o=len(o_atom_ids)
#n_n=len(n_atom_ids)

##print("C Atom IDs:", len(c_atom_ids))
##print("H Atom IDs:", len(h_atom_ids))
##print("O Atom IDs:", len(o_atom_ids))
##print("N Atom IDs:", len(n_atom_ids))
##print('total atom number', total_atom)
#m_c=12.011
#m_h=1.00784
#m_o=15.999
#m_n=14.0067
#m_total = m_c * n_c + m_h * n_h + m_o * n_o + m_n * n_n
#n_total = n_c + n_h + n_o + n_n
#avogard_number=6.02214076e23
#dry_weight=(n_total/avogard_number) * m_total
#print('total weight in gram', dry_weight)
ds_reverse= (1/(m_GlcNAc*Tempo_target*(10**-3))-(m_difference/m_GlcNAc) )  ##(1/ds)
ds=1/ds_reverse
n_coo=round(total_residues*ds)
#print('degree of substitution', ds)
#print('carboxylation_number',n_coo)
#
#
def adjust_for_carboxylation_limit(initial_carboxylation_num):
    initial_carboxylation_num = n_coo
    carboxylation_number_limit = c_iterations * planes_010_count + 2 * c_iterations * planes_100_count
    new_carboxylation_number = initial_carboxylation_num
    if initial_carboxylation_num > carboxylation_number_limit:
        new_carboxylation_number = carboxylation_number_limit
    return new_carboxylation_number

carboxylation_num = adjust_for_carboxylation_limit(n_coo)
#print(carboxylation_num)


if carboxylation_num > 0:
    left_010_count=0
    right_010_count=0
    upper_100_count=0
    bottom_100_count=0
    for i in range(1, carboxylation_num + 1):
        # Determine the target based on i % 4
        if i % 2 == 0:
            upper_100_count += 1
        elif i % 2 == 1:
            left_010_count += 1

    
    new_carboxylation_num=carboxylation_num
    
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
    

        
    
    new_carboxylation_number_update=upper_100_count_update  + left_010_count_update 
    actual_ds=(new_carboxylation_number_update)/(total_residues)
    actual_carboxylate_content=Tempo_target*(actual_ds/ds)

    #print(f"Actual Ds : {actual_ds}")
    #print(f"carboxylation_number: {new_carboxylation_number_update}")
    
    actual_carboxylate_rounded=round(actual_carboxylate_content,4)
    #print('surface charge density', actual_carboxylate_rounded, 'mmol/g')
    
    ##010 planes selection
    even_numbers_010 = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 0]
    odd_numbers_010 = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 1]
    
    all_combinations_010_left = [(value1, value2) for value1 in range(min(planes_010_left), max(planes_010_left) + 1)  for value2 in odd_numbers_010 ]
    all_combinations_010_right = [(value1, value2) for value1 in range(min(planes_010_right), max(planes_010_right) + 1)  for value2 in even_numbers_010 ]
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
                        indices_to_remove = residue.atoms[25:29].indices
                    elif residue.resid == maximum_resid_num:
                        indices_to_remove = residue.atoms[24:28].indices
                    else:
                        indices_to_remove = residue.atoms[23:27].indices
                    all_removed_indices.extend(indices_to_remove)
        unique_indices = sorted(set(all_removed_indices))
        indices_str = ' '.join(map(str, unique_indices))
        universe = universe.select_atoms(f"not index {indices_str}")
    
        return universe
    
    
    
    all_indices = [left_nacetyl_candidates_indices,  
                   upper_nacetyl_candidates_indices]
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
    left_coo_candidates_indices   = get_fragment_indices(primary_hydroxyl_remove_chains, sel_010_left_array)   
 
    upper_coo_candidates_indices  = get_fragment_indices(primary_hydroxyl_remove_chains, sel_100_upper_array)   

    primary_hydroxyl_remove_chains_u = mda.Universe("primary_hydroxyl_remove.temp.pdb")
    coo_atoms = list(primary_hydroxyl_remove_chains_u.atoms)

    def add_coo_based_on_indices(universe, indices_lists):
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
                    residue.resname = 'BLCP'
                    if residue.resid == 1 :
                       base_atom_index = 24   
                    elif residue.resid == maximum_resid_num :
                       base_atom_index = 23   
                    else :
                       base_atom_index = 22   
                    base_atom = residue.atoms[base_atom_index]
                    coo_new_positions = {} 
                    # Define new positions based on residue ID's odd/even status
                    if residue.resid % 2 == 1  and fragment_number in planes_100_even:  # Odd residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.047,   0.956,  -0.187]),
                            'O62': base_atom.position + np.array([-1.367,   0.411,  -0.081])
                        }
                    elif residue.resid % 2 == 0  and fragment_number in planes_100_even:  # Even residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([-1.047,  -0.956,  -0.187]),
                            'O62': base_atom.position + np.array([ 1.367,  -0.411,  -0.081])
                        }
    
                    elif residue.resid % 2 == 1  and fragment_number in planes_100_odd:  # Odd residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.367,  0.411,   0.081]),
                            'O62': base_atom.position + np.array([-1.047,  0.956,   0.187])
                        }
                    elif residue.resid % 2 == 0  and fragment_number in planes_100_odd:  # Even residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.047,  -0.956,   0.188]),
                            'O62': base_atom.position + np.array([-1.367,  -0.411,   0.081])
                        }
    
                    elif fragment_number in range(min(planes_010_left), max(planes_010_left)+1):  # Odd residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.047,  0.956,   -0.187]),
                            'O62': base_atom.position + np.array([-1.367,  0.411,   -0.081])
                        }
                    elif fragment_number in range(min(planes_010_right), max(planes_010_right)+1):  # Even residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.047,  -0.956,   0.188]),
                            'O62': base_atom.position + np.array([-1.367,  -0.411,   0.081])
                        }
                    insert_pos = coo_atoms.index(base_atom) + 1
                    for name, pos in coo_new_positions.items():
                        coo_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
                        coo_new_uni.add_TopologyAttr('name', [name])
                        coo_new_uni.add_TopologyAttr('type', [base_atom.type])
                        coo_new_uni.add_TopologyAttr('resname', [base_atom.resname])
                        coo_new_uni.add_TopologyAttr('resid', [residue.resid])
                        coo_new_uni.add_TopologyAttr('segid', [base_atom.segment.segid])
                        coo_new_uni.add_TopologyAttr('chainIDs', ['N'])
                        coo_new_uni.atoms.positions = [pos]
                        coo_new_atom = coo_new_uni.atoms[0]
    
                        coo_atoms.insert(insert_pos, coo_new_atom)
                        insert_pos += 1  # Update position for the next atom
        # Create a new universe with all atoms including the added ones
        new_universe = mda.Merge(*[mda.AtomGroup([atom]) for atom in coo_atoms])
        return new_universe
    
    
    all_indices = [left_coo_candidates_indices,  upper_coo_candidates_indices]
    nh3_u = add_coo_based_on_indices(primary_hydroxyl_remove_chains_u, all_indices)
    nh3_u.atoms.write("coo_negative-modified.temp.pdb")
    molecule.delete(primary_hydroxyl_remove_chains)


    #    
#    ####------------------------------------------------------------coo_ph building step----------------------------------------------------------------
    
    
    #pH calculation |  NH2_step
    ratio = 10 ** (pKa-pH)
    #print("ratio is :", ratio)

    y = (new_carboxylation_number_update) / (1 + ratio) ###coo- number
    #print(y)
    y_rounded = round(y) ###nh3 number for integer number
    #print(y_rounded)
    
    
    if  y_rounded == new_carboxylation_number_update:
        coo_negative = new_carboxylation_number_update
        coo_netural = 0
        actual_pH_rounded = pH
    elif y_rounded == 0:
        coo_negative = 0
        coo_netural = new_carboxylation_number_update
        actual_pH_rounded = pH
    elif y_rounded >0 and y_rounded < new_carboxylation_number_update:
        coo_negative = y_rounded
        coo_netural = int(new_carboxylation_number_update - y_rounded)
        exponential_part= coo_netural/coo_negative
        actual_pH_rounded = pKa- math.log10(exponential_part) 
    
    
    #print(f"coo netural number is:{coo_netural}, coo negative number is {coo_negative}")
    #print("actual pH value is:",actual_pH_rounded )
    
    
    
    def collect_carboxylation_atoms_to_remove_randomly(universe, indices_groups, coo_netural):
        all_indices = []
        atoms_to_remove = []  
        for group in indices_groups:
            all_indices.extend(group)
        np.random.shuffle(all_indices)
        for start_idx, end_idx in all_indices[:coo_netural]:
            possible_residues = universe.select_atoms(f"index {start_idx}:{end_idx}").residues
            if coo_netural ==0 :
                return []
            else:
                for residue in possible_residues:
                    residue.resname = 'BLCD'
        return atoms_to_remove
    
    def coo_netural_apply_modifications_and_save(universe, atoms_to_remove, output_name):
        if not atoms_to_remove:
            #print("No atoms to remove. Saving the original universe.")
            universe.atoms.write(output_name)
            return
        query = "not bynum " + " ".join(map(str, atoms_to_remove))
        remaining_atoms = universe.select_atoms(query)
        new_universe = mda.Merge(remaining_atoms)
        new_universe.atoms.write(output_name)
        #print("Modification complete. New PDB file written:", output_name)
    
    
    coo_deprotonation_chains = molecule.load("pdb", "coo_negative-modified.temp.pdb")
    left_coo_netural_candidates_indices   = get_fragment_indices(coo_deprotonation_chains, sel_010_left_array)     
    upper_coo_netural_candidates_indices  = get_fragment_indices(coo_deprotonation_chains, sel_100_upper_array)   

    

    coo_netural_indices_groups = [left_coo_netural_candidates_indices,  upper_coo_netural_candidates_indices]
    coo_negative_chains_u = mda.Universe("coo_negative-modified.temp.pdb")
    coo_negative_to_remove = collect_carboxylation_atoms_to_remove_randomly(coo_negative_chains_u, coo_netural_indices_groups, coo_netural)
    
    coo_netural_apply_modifications_and_save(coo_negative_chains_u, coo_negative_to_remove, "coo_netural-modified.temp.pdb")
    
    molecule.delete(coo_deprotonation_chains) 
    
    u_final = mda.Universe("coo_netural-modified.temp.pdb")
    box_x = a_trans * a_iterations + 20
    box_y = b_iterations * b_trans + 20
    box_z = c_iterations * c_trans + 20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    u_final.atoms.translate(translation_vector)
    with mda.Writer(f"charmm36-alpha-chitin-AB-fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
        
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    
    print(f"Carboxylate content: {actual_carboxylate_content:.2f}, Degree of carboxylation: {actual_ds:.4f}, pH: {actual_pH_rounded:.4f}")


elif carboxylation_num == 0:
    u_final = mda.Universe("alpha-chitin-A-temp.pdb")
    box_x = a_trans * a_iterations + 20
    box_y = b_iterations * b_trans + 20
    box_z = c_iterations * c_trans + 20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    with mda.Writer(f"charmm36-alpha-chitin-AB-fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
    

    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    actual_ds=0.0    
    actual_carboxylate_content=0.0
    actual_pH_rounded=0.0
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"Carboxylate content: {actual_carboxylate_content:.2f}, Degree of carboxylation: {actual_ds:.4f}, pH: {actual_pH_rounded:.4f}")
