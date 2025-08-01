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
    if a_trans <= 0 or b_trans <= 0 or c_trans <= 0:
        sys.stderr.write("Error: Please provide positive value for all crystallographic parameters.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length value provided. Please enter valid numeric values.\n")
    sys.exit(1)
except IndexError:
    sys.stderr.write("Error: Missing one or more input arguments. Please provide three integer values.\n")
    sys.exit(1)



try:
    alpha_angle=float(sys.argv[4])
    beta_angle =float(sys.argv[5])
    gamma_angle=float(sys.argv[6])
    v=float(sys.argv[7])  ##angstrom^3
    if alpha_angle <= 0 or beta_angle <= 0 or gamma_angle <= 0 or v<=0:
        sys.stderr.write("Error: Please provide positive value for all crystallographic parameters.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid angle and volume parameters. Please enter valid numeric values.\n")
    sys.exit(1)
except IndexError:
    sys.stderr.write("Error: Missing one or more input arguments. Please provide all required angle and volume parameters.\n")
    sys.exit(1)



try:
    c_iterations = int(sys.argv[8])
    if c_iterations < 0:
        sys.stderr.write("Error: Please provide a positive integer for c repetition number.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length value provide. Please enter a valid numeric value.\n")
    sys.exit(1)


try:
        ##carboxylation degree ref:/10.1016/j.carbpol.2019.115292
        Tempo_target = float(sys.argv[9])  #####
        if not (0 <= Tempo_target <= 2):
            sys.stderr.write("Error: Tempo_target must be between 0 and 2.\n")
            sys.exit(1)
        #if not (0 <= Tempo_target < 1):
        #    sys.stderr.write("Error: Tempo-oxidation degree  must be between 0 and 1.\n") 
        #    sys.exit(1)

        pH = float(sys.argv[10])
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
neutron_pdb_1 = os.path.join(main_folder_path, 'structure', 'cellulose_I_alpha', 'charmm36', 'chain.pdb')
user_defined_pdb_1 = os.path.join(main_folder_path, 'structure', 'cellulose_I_alpha', 'charmm36', 'chain-1_ud.pdb')

if os.path.exists(user_defined_pdb_1):
    unit_chain_input_file_1 = user_defined_pdb_1
else:
    unit_chain_input_file_1 = neutron_pdb_1


chain_1_u = mda.Universe(unit_chain_input_file_1)



#chain_1 assembly    
#------------------------------------------------------------------------------
chain_1_strand = []
for chain_1_i in range(1, c_iterations + 1):
    chain_1_u.atoms.positions += [c_trans, 0, 0]
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
    'O1': chain_1_first_atom.position + [1.257, -0.337, -0.733],
    'HO1': chain_1_first_atom.position + [1.107, -0.364, -1.481],
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



##------------------------------------------------------------------------------
##
###------------------------------------------------------------------------------
####translate parameter

alpha_angle_rad = (alpha_angle - 90) * math.pi / 180
cos_alpha_angle = math.cos(alpha_angle_rad)
sin_alpha_angle = math.sin(alpha_angle_rad)

beta_angle_rad = (beta_angle - 90) * math.pi / 180
cos_beta_angle = math.cos(beta_angle_rad)
sin_beta_angle = math.sin(beta_angle_rad)

beta_angle_raw = beta_angle * math.pi / 180
sin_beta_angle_raw = math.sin(beta_angle_raw)
alpha_angle_raw = alpha_angle * math.pi / 180
sin_alpha_angle_raw = math.sin(alpha_angle_raw)


gamma_angle_rad = (gamma_angle - 90) * math.pi / 180
cos_gamma_angle = math.cos(gamma_angle_rad)
sin_gamma_angle = math.sin(gamma_angle_rad)


a_par_vertical_move_1= a_trans * cos_beta_angle 
a_par_vertical_to_screen_move_1= -1 * a_trans * sin_beta_angle   


b_par_vertical_to_screen_move_1= -1 * b_trans * sin_alpha_angle   
b_par_vertical= -1 * b_trans * sin_gamma_angle 

h=b_trans * cos_alpha_angle  
z_dir_trans_move = v/(a_trans*c_trans*sin_beta_angle_raw)

try:
    sqrt_content_1=(h)**2- z_dir_trans_move**2
    if sqrt_content_1 < 0:
        raise ValueError("Invalid crystallographic parameters or unit volume size.")

    sqrt_content_2=b_trans**2 - (h**2)
    if sqrt_content_2 < 0:
        raise ValueError("Invalid crystallographic parameters or unit volume size.")

    z_dir_vertical_move = -1 * math.sqrt(sqrt_content_1)
    z_dir_vertical_to_screen = -1 *  math.sqrt(sqrt_content_2)

except ValueError as e:
    sys.stderr.write(f"Error with crystallographic parameters and unit volume size")
    sys.exit(1)


#chain_1_1st_layer assembly  
chain_1_1st_layer_input_file = "chain_1_temp.2.pdb"

chain_1_1st_layer = []  
for i in range(1, 4):
    chain_1_1st_layer_u = mda.Universe(chain_1_1st_layer_input_file)
    translation_vector = [ a_par_vertical_to_screen_move_1 * (i+1) ,
                          a_par_vertical_move_1*(i+1),  
                           0
                         ]
    chain_1_1st_layer_u.atoms.positions += translation_vector
    segid_increment_value = i 
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


##chain_1_2nd_layer assembly  
chain_1_2nd_layer_input_file = "chain_1_temp.2.pdb"
chain_1_2nd_layer_u = mda.Universe(chain_1_2nd_layer_input_file)
chain_1_2nd_layer = [] 

for j in range(1, 5):
    chain_1_2nd_layer_u = mda.Universe(chain_1_2nd_layer_input_file)


    translation_vector = [z_dir_vertical_to_screen + a_par_vertical_to_screen_move_1 * (j) ,
                          z_dir_vertical_move + a_par_vertical_move_1*(j),  
                          z_dir_trans_move
    ]


    chain_1_2nd_layer_u.atoms.positions += translation_vector
    for atom in chain_1_2nd_layer_u.atoms:
        new_numeric_part = 3 + j  
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





##chain_1_3rd_layer assembly  
chain_1_3rd_layer_input_file = "chain_1_temp.2.pdb"
chain_1_3rd_layer_u = mda.Universe(chain_1_3rd_layer_input_file)
chain_1_3rd_layer = [] 

for j in range(1, 5):
    chain_1_3rd_layer_u = mda.Universe(chain_1_3rd_layer_input_file)
    translation_vector = [ z_dir_vertical_to_screen * 2 + a_par_vertical_to_screen_move_1 * (j) ,
                          z_dir_vertical_move * 2 + a_par_vertical_move_1*(j),  
                          z_dir_trans_move * 2
                         ]
    chain_1_3rd_layer_u.atoms.positions += translation_vector
    for atom in chain_1_3rd_layer_u.atoms:
        new_numeric_part = 7 + j  
        atom.segment.segid = str(new_numeric_part)  
    chain_1_3rd_layer_output = f"chain_1_3rd_layer_{j}.pdb"
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



##chain_1_4th_layer assembly  
chain_1_4th_layer_input_file = "chain_1_temp.2.pdb"
chain_1_4th_layer_u = mda.Universe(chain_1_4th_layer_input_file)
chain_1_4th_layer = [] 

for j in range(1, 5):
    chain_1_4th_layer_u = mda.Universe(chain_1_4th_layer_input_file)
    translation_vector = [ z_dir_vertical_to_screen * 3 + a_par_vertical_to_screen_move_1 * (j) ,
                          z_dir_vertical_move * 3 + a_par_vertical_move_1*(j),  
                          z_dir_trans_move * 3
                         ]
    chain_1_4th_layer_u.atoms.positions += translation_vector
    for atom in chain_1_4th_layer_u.atoms:
        new_numeric_part = 11 + j  
        atom.segment.segid = str(new_numeric_part)  
    chain_1_4th_layer_output = f"chain_1_4th_layer_{j}.pdb"
    chain_1_4th_layer.append(chain_1_4th_layer_output)
    with mda.Writer(chain_1_4th_layer_output, n_atoms=chain_1_4th_layer_u.atoms.n_atoms) as W:
        W.write(chain_1_4th_layer_u.atoms)
with open("chain_1_4th_layer_temp.pdb", "w") as layer_file:
    for chain_1_4th_layer_output in chain_1_4th_layer:
        with open(chain_1_4th_layer_output, "r") as chain_1_4th_layer_pdb_file:
            for line in chain_1_4th_layer_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
            os.remove(chain_1_4th_layer_output)



##chain_1_5th_layer assembly  
chain_1_5th_layer_input_file = "chain_1_temp.2.pdb"
chain_1_5th_layer_u = mda.Universe(chain_1_5th_layer_input_file)
chain_1_5th_layer = [] 

for j in range(1, 4):
    chain_1_5th_layer_u = mda.Universe(chain_1_5th_layer_input_file)
    translation_vector = [ z_dir_vertical_to_screen * 4 + a_par_vertical_to_screen_move_1 * (j) ,
                          z_dir_vertical_move * 4 + a_par_vertical_move_1*(j),  
                          z_dir_trans_move * 4
                         ]
    chain_1_5th_layer_u.atoms.positions += translation_vector
    for atom in chain_1_5th_layer_u.atoms:
        new_numeric_part = 15 + j  
        atom.segment.segid = str(new_numeric_part)  
    chain_1_5th_layer_output = f"chain_1_5th_layer_{j}.pdb"
    chain_1_5th_layer.append(chain_1_5th_layer_output)
    with mda.Writer(chain_1_5th_layer_output, n_atoms=chain_1_5th_layer_u.atoms.n_atoms) as W:
        W.write(chain_1_5th_layer_u.atoms)
with open("chain_1_5th_layer_temp.pdb", "w") as layer_file:
    for chain_1_5th_layer_output in chain_1_5th_layer:
        with open(chain_1_5th_layer_output, "r") as chain_1_5th_layer_pdb_file:
            for line in chain_1_5th_layer_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
            os.remove(chain_1_5th_layer_output)



####assemble
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
    "chain_1_4th_layer_temp.pdb", "chain_1_5th_layer_temp.pdb"]

unit_temp_pdb = "unit_temp.pdb"
assemble_pdbs(unit_temp_pdb, chain_temp_pdb_files)




pKa=3.25                        # pka of coo-
m_Glc=163.09316                 # Glc unit mass
m_difference=111                # difference between Glc unit with Glc-cooNa(3)




##crystallographic planes defined
planes_010_left=[15,16,17]
planes_010_right=[0,1,2]

planes_001_left=[3,7,11]
planes_001_right=[6,10,14]

planes_010_count=len(planes_010_left)+len(planes_010_right)
planes_001_count=len(planes_001_left)+len(planes_001_right)

chain_u = mda.Universe("unit_temp.pdb")
residue_ids = [residue.resid for residue in chain_u.residues]
maximum_resid_num=residue_ids[-1]
total_residues=len(residue_ids)


ds_reverse= (1/(m_Glc*Tempo_target*(10**-3))-(m_difference/m_Glc) )  ##(1/ds)
ds=1/ds_reverse
n_coo=round(total_residues*ds)
#print('degree of substitution', ds)

if n_coo < 0 :
    n_coo=total_residues



carboxylation_num = n_coo
#print(carboxylation_num)


if carboxylation_num > 0:
    left_010_count=0
    right_010_count=0
    left_001_count=0
    right_001_count=0
    for i in range(1, carboxylation_num + 1):
        # Determine the target based on i % 4
        if i % 4 == 0:
            left_010_count += 1
        elif i % 4 == 1:
            right_010_count += 1
        elif i % 4 == 2:
            left_001_count += 1
        elif i % 4 == 3:
            right_001_count += 1
    
    max_010_count_left = c_iterations * len(planes_010_left)
    max_010_count_right = c_iterations * len(planes_010_right)

    max_001_count_left = c_iterations * len(planes_001_left)
    max_001_count_right = c_iterations * len(planes_001_right)

    left_010_count_update=left_010_count
    right_010_count_update=right_010_count
    left_001_count_update=left_001_count
    right_001_count_update=right_001_count

    if left_010_count > max_010_count_left:
        left_010_count_update = max_010_count_left
    
    if right_010_count > max_010_count_right:
        right_010_count_update = max_010_count_right
    
    if left_001_count > max_001_count_left:
        left_001_count_update = max_001_count_left
    
    if right_001_count > max_001_count_right:
        right_001_count_update = max_001_count_right    


    new_carboxylation_number_update=left_010_count_update + right_010_count_update + left_001_count_update + right_001_count_update 
    #print(new_carboxylation_number_update)
    #print(total_residues)
    actual_ds=(new_carboxylation_number_update)/(total_residues)
    actual_carboxylate_content=Tempo_target*(actual_ds/ds)
    #print(actual_ds)
    #print(actual_carboxylate_content)
    actual_carboxylate_rounded=round(actual_carboxylate_content,4)
    #print('surface charge density', actual_carboxylate_rounded, 'mmol/g')

    even_numbers = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 0]
    odd_numbers = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 1]
    all_combinations_010_left = [(value1, value2) for value1 in planes_010_left  for value2 in odd_numbers]
    all_combinations_010_right = [(value1, value2) for value1 in planes_010_right  for value2 in even_numbers ]
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
    
    while len(sel_010_right_array) < right_010_count_update:
        if not all_combinations_010_right:  # Check if the list is empty to avoid IndexError
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_010_right = random.choice(all_combinations_010_right)
        sel_010_right_array.append(selected_combination_010_right)
        all_combinations_010_right.remove(selected_combination_010_right)

    all_combinations_001_left = [(value1, value2) for value1 in planes_001_left  for value2 in odd_numbers]
    all_combinations_001_right = [(value1, value2) for value1 in planes_001_right  for value2 in even_numbers ]
    random.shuffle(all_combinations_001_left)
    random.shuffle(all_combinations_001_right)
    sel_001_left_array = []
    sel_001_right_array = []

    while len(sel_001_left_array) < left_001_count_update:
        if not all_combinations_001_left:  
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_001_left = random.choice(all_combinations_001_left)
        sel_001_left_array.append(selected_combination_001_left)
        all_combinations_001_left.remove(selected_combination_001_left)
    
    while len(sel_001_right_array) < right_001_count_update:
        if not all_combinations_001_right:  # Check if the list is empty to avoid IndexError
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_001_right = random.choice(all_combinations_001_right)
        sel_001_right_array.append(selected_combination_001_right)
        all_combinations_001_right.remove(selected_combination_001_right)

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
    left_glc_010_candidates_indices  = get_fragment_indices(origin_glc_chains, sel_010_left_array)   
    right_glc_010_candidates_indices = get_fragment_indices(origin_glc_chains, sel_010_right_array)   
    left_glc_001_candidates_indices  = get_fragment_indices(origin_glc_chains, sel_001_left_array)   
    right_glc_001_candidates_indices = get_fragment_indices(origin_glc_chains, sel_001_right_array)  
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
    
    all_indices = [left_glc_010_candidates_indices, right_glc_010_candidates_indices, 
                   left_glc_001_candidates_indices, right_glc_001_candidates_indices]
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
    right_coo_candidates_indices  = get_fragment_indices(primary_hydroxyl_remove_chains, sel_010_right_array)   
    upper_coo_candidates_indices  = get_fragment_indices(primary_hydroxyl_remove_chains, sel_001_left_array)   
    bottom_coo_candidates_indices = get_fragment_indices(primary_hydroxyl_remove_chains, sel_001_right_array)  
    primary_hydroxyl_remove_chains_u = mda.Universe("primary_hydroxyl_remove.temp.pdb")
    coo_atoms = list(primary_hydroxyl_remove_chains_u.atoms)#

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
                    residue.resname = 'BGLA'
                    if residue.resid == 1 :
                       #print(residue.resid)
                       base_atom_index = 18   
                    elif residue.resid == maximum_resid_num :
                       #print(residue.resid)
                       base_atom_index = 17   
                    else :
                       #print(residue.resid)
                       base_atom_index = 16   
                    base_atom = residue.atoms[base_atom_index]
                    coo_new_positions = {} 
                    # Define new positions based on residue ID's odd/even status
                    if residue.resid % 2 == 1  and fragment_number in planes_001_left:  # Odd residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.345,  -0.473,  -0.114]),
                            'O62': base_atom.position + np.array([-1.093,  -0.896,  -0.216])
                        }
                    elif residue.resid % 2 == 0  and fragment_number in planes_001_right:  # Even residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.093,   0.896,  -0.216]),
                            'O62': base_atom.position + np.array([-1.345,   0.473,  -0.114])
                        }
                    elif residue.resid % 2 == 1  and fragment_number in planes_010_left:  # Odd residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.345,  -0.473,  -0.114]),
                            'O62': base_atom.position + np.array([-1.093,  -0.896,  -0.216])
                        }
                    elif residue.resid % 2 == 0  and fragment_number in planes_010_right:  # Even residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.093,   0.896,  -0.216]),
                            'O62': base_atom.position + np.array([-1.345,   0.473,  -0.114])
                        }
                    insert_pos = coo_atoms.index(base_atom) + 1
                    for name, pos in coo_new_positions.items():
                        coo_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
                        coo_new_uni.add_TopologyAttr('name', [name])
                        coo_new_uni.add_TopologyAttr('type', [base_atom.type])
                        coo_new_uni.add_TopologyAttr('resname', [base_atom.resname])
                        coo_new_uni.add_TopologyAttr('resid', [residue.resid])
                        coo_new_uni.add_TopologyAttr('segid', [base_atom.segment.segid])
                        coo_new_uni.add_TopologyAttr('chainIDs', ['X'])
                        coo_new_uni.atoms.positions = [pos]
                        coo_new_atom = coo_new_uni.atoms[0]
                        coo_atoms.insert(insert_pos, coo_new_atom)
                        insert_pos += 1  # Update position for the next atom
        new_universe = mda.Merge(*[mda.AtomGroup([atom]) for atom in coo_atoms])
        return new_universe
    

    all_indices = [left_coo_candidates_indices, right_coo_candidates_indices, upper_coo_candidates_indices, bottom_coo_candidates_indices]
    coo_negative_u = add_coo_based_on_indices(primary_hydroxyl_remove_chains_u, all_indices)
    coo_negative_u.atoms.write("coo_negative-modified.temp.pdb")
    molecule.delete(primary_hydroxyl_remove_chains)#

   
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
                    residue.resname = 'BGLD'
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
    right_coo_netural_candidates_indices  = get_fragment_indices(coo_deprotonation_chains, sel_010_right_array)   
    upper_coo_netural_candidates_indices  = get_fragment_indices(coo_deprotonation_chains, sel_001_left_array)   
    bottom_coo_netural_candidates_indices = get_fragment_indices(coo_deprotonation_chains, sel_001_right_array)  
    

    coo_netural_indices_groups = [left_coo_netural_candidates_indices, right_coo_netural_candidates_indices, upper_coo_netural_candidates_indices, bottom_coo_netural_candidates_indices]
    coo_negative_chains_u = mda.Universe("coo_negative-modified.temp.pdb")
    coo_negative_to_remove = collect_carboxylation_atoms_to_remove_randomly(coo_negative_chains_u, coo_netural_indices_groups, coo_netural)
    
    coo_netural_apply_modifications_and_save(coo_negative_chains_u, coo_negative_to_remove, "coo_netural-modified.temp.pdb")
    
    molecule.delete(coo_deprotonation_chains) 
    
    u_final = mda.Universe("coo_netural-modified.temp.pdb")
    box_x = 50
    box_y = 50
    box_z = c_iterations * c_trans + 20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    u_final.atoms.translate(translation_vector)
    with mda.Writer(f"charmm36-cellulose-Ialpha-18fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
        
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    
    print(f"Carboxylate content: {actual_carboxylate_content:.2f}, Degree of carboxylation: {actual_ds:.4f}, pH: {actual_pH_rounded:.4f}")


elif carboxylation_num == 0:
    u_final = mda.Universe("unit_temp.pdb")
    box_x = 50
    box_y = 50
    box_z = c_iterations * c_trans + 20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    with mda.Writer(f"charmm36-cellulose-Ialpha-18fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
    

    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    actual_ds=0.0    
    actual_carboxylate_content=0.0
    actual_pH_rounded=0.0
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"Carboxylate content: {actual_carboxylate_content:.2f}, Degree of carboxylation: {actual_ds:.4f}, pH: {actual_pH_rounded:.4f}")
