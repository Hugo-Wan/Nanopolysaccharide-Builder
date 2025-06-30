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


m_Glc=163.09316                # Glc unit mass
m_difference=101.1             # difference between Glc unit with Glc-SO3- Na(+)



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
        ##Sulfate degree ref:/10.1016/j.carbpol.2019.115292
        Sulfate_target = float(sys.argv[9])  #####
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
chain_1_1st_layer_u = mda.Universe(chain_1_1st_layer_input_file)
chain_1_1st_layer = []  
for i in range(1, 7):
    chain_1_1st_layer_u.atoms.positions += [ a_par_vertical_to_screen_move_1, a_par_vertical_move_1, 0]

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


##chain_1_2nd_layer assembly  
chain_1_2nd_layer_input_file = "chain_1_temp.2.pdb"
chain_1_2nd_layer_u = mda.Universe(chain_1_2nd_layer_input_file)
chain_1_2nd_layer = [] 

for j in range(1, 7):
    chain_1_2nd_layer_u = mda.Universe(chain_1_2nd_layer_input_file)


    translation_vector = [z_dir_vertical_to_screen + a_par_vertical_to_screen_move_1 * (j) ,
                          z_dir_vertical_move + a_par_vertical_move_1*(j),  
                          z_dir_trans_move
    ]


    chain_1_2nd_layer_u.atoms.positions += translation_vector
    for atom in chain_1_2nd_layer_u.atoms:
        new_numeric_part = 6 + j  
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

for j in range(1, 7):
    chain_1_3rd_layer_u = mda.Universe(chain_1_3rd_layer_input_file)
    translation_vector = [ z_dir_vertical_to_screen * 2 + a_par_vertical_to_screen_move_1 * (j) ,
                          z_dir_vertical_move * 2 + a_par_vertical_move_1*(j),  
                          z_dir_trans_move * 2
                         ]
    chain_1_3rd_layer_u.atoms.positions += translation_vector
    for atom in chain_1_3rd_layer_u.atoms:
        new_numeric_part = 12 + j  
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

for j in range(1, 6):
    chain_1_4th_layer_u = mda.Universe(chain_1_4th_layer_input_file)
    translation_vector = [ z_dir_vertical_to_screen * 3 + a_par_vertical_to_screen_move_1 * (j) ,
                          z_dir_vertical_move * 3 + a_par_vertical_move_1*(j),  
                          z_dir_trans_move * 3
                         ]
    chain_1_4th_layer_u.atoms.positions += translation_vector
    for atom in chain_1_4th_layer_u.atoms:
        new_numeric_part = 18 + j  
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

for j in range(1, 5):
    chain_1_5th_layer_u = mda.Universe(chain_1_5th_layer_input_file)
    translation_vector = [ z_dir_vertical_to_screen * 4 + a_par_vertical_to_screen_move_1 * (j) ,
                          z_dir_vertical_move * 4 + a_par_vertical_move_1*(j),  
                          z_dir_trans_move * 4
                         ]
    chain_1_5th_layer_u.atoms.positions += translation_vector
    for atom in chain_1_5th_layer_u.atoms:
        new_numeric_part = 23 + j  
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





#chain_1_6th_layer assembly  
chain_1_6th_layer_input_file = "chain_1_temp.2.pdb"
chain_1_6th_layer_u = mda.Universe(chain_1_6th_layer_input_file)
chain_1_6th_layer = [] 

for j in range(1, 6):
    chain_1_6th_layer_u = mda.Universe(chain_1_6th_layer_input_file)
    translation_vector = [ z_dir_vertical_to_screen * -1 + a_par_vertical_to_screen_move_1 * (j+1) ,
                          z_dir_vertical_move * -1 + a_par_vertical_move_1*(j+1),  
                          z_dir_trans_move * -1
                         ]
    chain_1_6th_layer_u.atoms.positions += translation_vector
    for atom in chain_1_6th_layer_u.atoms:
        new_numeric_part = 27 + j  
        atom.segment.segid = str(new_numeric_part)  
    chain_1_6th_layer_output = f"chain_1_6th_layer_{j}.pdb"
    chain_1_6th_layer.append(chain_1_6th_layer_output)
    with mda.Writer(chain_1_6th_layer_output, n_atoms=chain_1_6th_layer_u.atoms.n_atoms) as W:
        W.write(chain_1_6th_layer_u.atoms)
with open("chain_1_6th_layer_temp.pdb", "w") as layer_file:
    for chain_1_6th_layer_output in chain_1_6th_layer:
        with open(chain_1_6th_layer_output, "r") as chain_1_6th_layer_pdb_file:
            for line in chain_1_6th_layer_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
            os.remove(chain_1_6th_layer_output)





#chain_1_7th_layer assembly  
chain_1_7th_layer_input_file = "chain_1_temp.2.pdb"
chain_1_7th_layer_u = mda.Universe(chain_1_7th_layer_input_file)
chain_1_7th_layer = [] 

for j in range(1, 5):
    chain_1_7th_layer_u = mda.Universe(chain_1_7th_layer_input_file)
    translation_vector = [ z_dir_vertical_to_screen * -2 + a_par_vertical_to_screen_move_1 * (j+2) ,
                          z_dir_vertical_move * -2 + a_par_vertical_move_1*(j+2),  
                          z_dir_trans_move * -2
                         ]
    chain_1_7th_layer_u.atoms.positions += translation_vector
    for atom in chain_1_7th_layer_u.atoms:
        new_numeric_part = 32 + j  
        atom.segment.segid = str(new_numeric_part)  
    chain_1_7th_layer_output = f"chain_1_7th_layer_{j}.pdb"
    chain_1_7th_layer.append(chain_1_7th_layer_output)
    with mda.Writer(chain_1_7th_layer_output, n_atoms=chain_1_7th_layer_u.atoms.n_atoms) as W:
        W.write(chain_1_7th_layer_u.atoms)
with open("chain_1_7th_layer_temp.pdb", "w") as layer_file:
    for chain_1_7th_layer_output in chain_1_7th_layer:
        with open(chain_1_7th_layer_output, "r") as chain_1_7th_layer_pdb_file:
            for line in chain_1_7th_layer_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
            os.remove(chain_1_7th_layer_output)



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
    "chain_1_4th_layer_temp.pdb", "chain_1_5th_layer_temp.pdb",  "chain_1_6th_layer_temp.pdb", "chain_1_7th_layer_temp.pdb"]

unit_temp_pdb = "unit_temp.pdb"
assemble_pdbs(unit_temp_pdb, chain_temp_pdb_files)




pKa=3.25                        # pka of coo-
m_Glc=163.09316                 # Glc unit mass
m_difference=111                # difference between Glc unit with Glc-cooNa(3)




##crystallographic planes defined
planes_010_left=[23, 24, 25, 26]


planes_010_right=[32, 33, 34, 35]

planes_001_left=[0,6,12,18]
planes_001_right=[11, 17, 5, 31]

planes_010_count=len(planes_010_left)+len(planes_010_right)
planes_001_count=len(planes_001_left)+len(planes_001_right)

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

sulfate_num = n_so3
#print(sulfate_num)

if sulfate_num > 0:
    left_001_count=sulfate_num

    max_010_count_left = c_iterations * len(planes_010_left)
    max_010_count_right = c_iterations * len(planes_010_right)

    max_001_count_left = c_iterations * len(planes_001_left)
    max_001_count_right = c_iterations * len(planes_001_right)

    left_001_count_update=left_001_count

    if left_001_count > max_001_count_left:
        left_001_count_update = max_001_count_left
    


    new_sulfate_number_update= left_001_count_update  
    #print(new_carboxylation_number_update)
    #print(total_residues)
    #print(total_residues)
    actual_ds=(new_sulfate_number_update)/(total_residues)
    actual_sulfate_content=Sulfate_target*(actual_ds/ds)
    #print(actual_ds)
    #print(actual_sulfate_content)
    actual_sulfate_rounded=round(actual_sulfate_content,4)
    #print('surface charge density', actual_sulfate_rounded, 'mmol/g')

    even_numbers = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 0]
    odd_numbers = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 1]

    all_combinations_001_left = [(value1, value2) for value1 in planes_001_left  for value2 in odd_numbers]
    all_combinations_001_right = [(value1, value2) for value1 in planes_001_right  for value2 in even_numbers ]
    random.shuffle(all_combinations_001_left)
    random.shuffle(all_combinations_001_right)
    sel_001_left_array = []


    while len(sel_001_left_array) < left_001_count_update:
        if not all_combinations_001_left:  
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_001_left = random.choice(all_combinations_001_left)
        sel_001_left_array.append(selected_combination_001_left)
        all_combinations_001_left.remove(selected_combination_001_left)
    
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
    left_glc_001_candidates_indices  = get_fragment_indices(origin_glc_chains, sel_001_left_array)   
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
    
    all_indices = [left_glc_001_candidates_indices ]
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
    upper_so3_candidates_indices   = get_fragment_indices(primary_hydroxyl_remove_chains, sel_001_left_array)   
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
                    if fragment_number in planes_001_left:  # Odd residue ID
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([ 0.353, -0.613, -0.803000000000001]),
                            'H62':  base_atom.position + np.array([-1.069,  0.004, -0.044000000000001]),                            
                            'O6':   base_atom.position + np.array([ 0.369, -0.647, 1.22]),
                            'S6':   base_atom.position + np.array([ 0.174, -2.031, 2.134]),
                            'OS62': base_atom.position + np.array([ 1.544, -2.293, 3.05]),
                            'OS63': base_atom.position + np.array([-1.145, -1.849, 3.143]),
                            'OS64': base_atom.position + np.array([-0.076, -3.335, 1.121])
                        }
                    elif fragment_number in planes_001_right:  # Even residue ID
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([-0.889,   -0.033,  -0.202]),
                            'H62':  base_atom.position + np.array([ 0.104,    0.316,   0.911]),
                            'O6':   base_atom.position + np.array([ 0.640,    0.987,  -0.787]),
                            'S6':   base_atom.position + np.array([ 1.396,    2.155,  -1.715]),
                            'OS62': base_atom.position + np.array([ 2.791,    2.323,  -1.193]),
                            'OS63': base_atom.position + np.array([ 1.343,    1.709,  -3.323]),
                            'OS64': base_atom.position + np.array([ 0.612,    3.562,  -1.522])
                        }
                    elif fragment_number in planes_010_left:  # Odd residue ID
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([ 0.353, -0.613, -0.803000000000001]),
                            'H62':  base_atom.position + np.array([-1.069,  0.004, -0.044000000000001]),                            
                            'O6':   base_atom.position + np.array([ 0.369, -0.647, 1.22]),
                            'S6':   base_atom.position + np.array([ 0.174, -2.031, 2.134]),
                            'OS62': base_atom.position + np.array([ 1.544, -2.293, 3.05]),
                            'OS63': base_atom.position + np.array([-1.145, -1.849, 3.143]),
                            'OS64': base_atom.position + np.array([-0.076, -3.335, 1.121])
                        }
                    elif fragment_number in planes_010_right:  # Even residue ID
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([-0.889,   -0.033,  -0.202]),
                            'H62':  base_atom.position + np.array([ 0.104,    0.316,   0.911]),
                            'O6':   base_atom.position + np.array([ 0.640,    0.987,  -0.787]),
                            'S6':   base_atom.position + np.array([ 1.396,    2.155,  -1.715]),
                            'OS62': base_atom.position + np.array([ 2.791,    2.323,  -1.193]),
                            'OS63': base_atom.position + np.array([ 1.343,    1.709,  -3.323]),
                            'OS64': base_atom.position + np.array([ 0.612,    3.562,  -1.522])
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
    

    all_indices = [upper_so3_candidates_indices]
    so3_negative_u = add_so3_based_on_indices(primary_hydroxyl_remove_chains_u, all_indices)
    so3_negative_u.atoms.write("so3_negative-modified.temp.pdb")
    molecule.delete(primary_hydroxyl_remove_chains)#
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
    with mda.Writer(f"charmm36-cellulose-Ialpha-36fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
        
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    
    print(f"Sulfate content: {actual_sulfate_content:.2f}, Degree of sulfation: {actual_ds:.4f}")


elif sulfate_num == 0:
    u_final = mda.Universe("unit_temp.pdb")
    box_x = 50
    box_y = 50
    box_z = c_iterations * c_trans + 20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    with mda.Writer(f"charmm36-cellulose-Ialpha-36fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
    

    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    actual_ds=0.0    
    actual_sulfate_content=0.0
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"Carboxylate content: {actual_sulfate_content:.2f}, Degree of carboxylation: {actual_ds:.4f}")

