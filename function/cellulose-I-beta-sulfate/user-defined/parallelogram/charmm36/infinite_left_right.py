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




m_Glc=163.09316                # Glc unit mass
m_difference=101.1             # difference between Glc unit with Glc-SO3- Na(+)

try:
    c_iterations = int(sys.argv[5])
    if c_iterations < 0:
        sys.stderr.write("Error: Please provide a positive integer for c repetition number.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length value provide. Please enter a valid numeric value.\n")
    sys.exit(1)
except IndexError:
    sys.stderr.write("Error: Missing one or more input arguments. Please provide three integer values.\n")
    sys.exit(1)

    
try:
    height = float(sys.argv[6])
    num_h = int(round(height /7.8))
    if num_h < 1:
        sys.stderr.write("Error: Please provide an larger height value.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid height value provide. Please enter a valid numeric value.\n")
    sys.exit(1)    


try:
    width = float(sys.argv[7])
    num_w = int((round(width/5.3)/2))
    if num_w < 1:
        sys.stderr.write("Error: Please provide an larger width value.\n")
        sys.exit(1)


        sys.exit(1)        
except ValueError:
    sys.stderr.write("Error: Invalid width value provide. Please enter a valid numeric value.\n")
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

#DDA_target=0.20      



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

chain_1_combined=mda.Universe("chain_1_temp.pdb")
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


chain_2_combined=mda.Universe("chain_2_temp.pdb")
for atom in chain_2_combined.atoms:
    atom.segment.segid = '0'  
    elements = [atom.name[0] for atom in chain_2_combined.atoms]
    chain_2_combined.add_TopologyAttr(Elements(elements))
chain_2_combined.atoms.write("chain_2_temp.2.pdb")



#------------------------------------------------------------------------------
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


#chain_1_layer assembly  
chain_1_1st_layer_input_file = "chain_1_temp.2.pdb"
chain_2_1st_layer_input_file = "chain_2_temp.2.pdb"
chain_1_1st_layer_u = mda.Universe(chain_1_1st_layer_input_file)
chain_2_1st_layer_u = mda.Universe(chain_2_1st_layer_input_file)
chain_1_1st_layer = []  
chain_2_1st_layer = []  
##chain_1_layer assembly  
#chain_1_1st_layer_input_file = "chain_1_temp.2.pdb"
#chain_1_1st_layer_u = mda.Universe(chain_1_1st_layer_input_file)
#chain_1_1st_layer = []  
#
for i in range(1, num_w+1):
    j=i-1
    chain_1_1st_layer_u.atoms.positions += [ a_par_vertical_move + b_par_vertical_move, a_par_transverse_move+ b_par_transverse_move, 0]
    segid_increment_value = i 
    for atom in chain_1_1st_layer_u.atoms:
        numeric_part = atom.segid.replace('','')
        new_numeric_part = 0 + segid_increment_value
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


for i in range(1, num_w+1):
    j=i-1
    chain_2_1st_layer_u.atoms.positions += [ a_par_vertical_move+b_par_vertical_move, a_par_transverse_move+b_par_transverse_move, 0]
    segid_increment_value = 1 
    for atom in chain_2_1st_layer_u.atoms:
        numeric_part = atom.segid.replace('','')
        new_numeric_part = num_w + i
    new_segid = f"{new_numeric_part}"
    atom.segment.segid = new_segid

    chain_2_1st_layer_output = f"chain_2_1st_layer_{i}.pdb"
    chain_2_1st_layer.append(chain_2_1st_layer_output)
    with mda.Writer(chain_2_1st_layer_output, n_atoms=chain_2_1st_layer_u.atoms.n_atoms) as W:
        W.write(chain_2_1st_layer_u.atoms)
with open("chain_2_1st_layer_temp.pdb", "w") as layer_file:
    for chain_2_1st_layer_output in chain_2_1st_layer:
        with open(chain_2_1st_layer_output, "r") as chain_2__1st_layer_pdb_file:
            for chain_2_1st_layer_line in chain_2__1st_layer_pdb_file:
                if chain_2_1st_layer_line.startswith("ATOM"):
                    layer_file.write(chain_2_1st_layer_line)
            os.remove(chain_2_1st_layer_output) 


def assemble_pdbs(output_file, input_files):
    with open(output_file, 'w') as outfile:
        for filename in input_files:
            with open(filename, 'r') as infile:
                for line in infile:
                    if line.startswith("ATOM"):
                        outfile.write(line) 
chain_temp_pdb_files = ["chain_1_1st_layer_temp.pdb", "chain_2_1st_layer_temp.pdb"]
#
unit_temp_pdb = "layer-1-temp.pdb"
assemble_pdbs(unit_temp_pdb, chain_temp_pdb_files)   

layer_after_layer=[]
for i in range(1, num_h+1):  
    layer_after_layer_u = mda.Universe("layer-1-temp.pdb")       
    translation_vector = [a_par_vertical_move*i, a_par_transverse_move*i, 0]
    layer_after_layer_u.atoms.positions += translation_vector
    segid_increment_value = (i-1)*2*num_w
    for segid, group in layer_after_layer_u.atoms.groupby('segids').items():
        new_segid = str(int(segid) + segid_increment_value)
        for atom in group:
            atom.segment.segid = new_segid

    layer_after_layer_output = f"layer_after_{i}.pdb"
    layer_after_layer.append(layer_after_layer_output)
    with mda.Writer(layer_after_layer_output, n_atoms=layer_after_layer_u.atoms.n_atoms) as W:
        W.write(layer_after_layer_u.atoms)
with open("layer_after_temp.pdb", "w") as layer_file:
    for layer_after_layer_output in layer_after_layer:
        with open(layer_after_layer_output, "r") as layer_after_layer_pdb_file:
            for line in layer_after_layer_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
        os.remove(layer_after_layer_output)  







##crystallographic planes defined
planes_left = []
planes_right = []


##-------------------crystallographic-plane-definitions---------------------
for m in range(1, num_h + 1):
    if m == 1:
        value_1_1=0
        value1 = num_w
        value2 = num_w - 1
        planes_right.extend([value1, value_1_1])
        planes_left.extend([value2])
    elif m>1 and m<num_h:
        value1 = num_w + 2 * (m - 1) * num_w
        value2 = (num_w - 1) + 2 * (m - 1) * num_w
        value3 = num_w + 2 * (m - 1) * num_w + num_w - 1
        value4 = num_w * 2 * (m - 1)
        planes_right.extend([value1, value4])
        planes_left.extend([value2, value3])
    elif m==num_h:
        value1 = num_w + 2 * (m - 1) * num_w
        value2 = (num_w - 1) + 2 * (m - 1) * num_w
        value_2_1=value2+num_w
        planes_right.extend([value1])
        planes_left.extend([value2, value_2_1])

left_right_count=len(planes_right)+len(planes_left)


range1_start = 0
range1_end = num_w * 2 - 1
planes_110_upper = [x for x in range(range1_start, range1_end + 1) if x not in (num_w, num_w - 1)]
range2_start = num_w*2*num_h-1 - (num_w * 2 - 1)
range2_end = num_w*2*num_h-1
planes_110_bottom = [x for x in range(range2_start, range2_end + 1) if x not in (num_w*2*num_h-(num_w), num_w*2*num_h-(num_w )-1)]


#print("right", planes_right)
#print("left",  planes_left)
#print("110 upper", planes_110_upper)
#print("110 bottom", planes_110_bottom)

planes_110_count=len(planes_110_upper)+len(planes_110_bottom)

left_right_count=len(planes_right)+len(planes_left)



chain_u = mda.Universe("layer_after_temp.pdb")
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




sulfate_num = n_so3


if sulfate_num > 0:
    left_count=0
    right_count=0
    for i in range(1, sulfate_num + 1):
        # Determine the target based on i % 4
        if i % 2 == 0:
            left_count += 1
        elif i % 2 == 1:
            right_count += 1
    
    planes_110_upper_count  = len(planes_110_upper)
    planes_110_bottom_count = len(planes_110_bottom)
    planes_left_count   =  len(planes_left)
    planes_right_count  =  len(planes_right)


    max_count_left = c_iterations * planes_left_count 
    max_count_right = c_iterations * planes_right_count
    

    max_110_count_upper = c_iterations * planes_110_upper_count
    max_110_count_bottom = c_iterations * planes_110_bottom_count


    left_count_update = left_count
    right_count_update = right_count

    if left_count > max_count_left:
        left_count_update = max_count_left
    
    if right_count > max_count_right:
        right_count_update = max_count_right

    new_sulfate_number_update= left_count_update + right_count_update 
    actual_ds=(new_sulfate_number_update)/(total_residues)
    actual_sulfate_content=Sulfate_target*(actual_ds/ds)
    #print(actual_ds)
    #print(actual_sulfate_content)
    actual_sulfate_rounded=round(actual_sulfate_content,4)
    #print('surface charge density', actual_sulfate_rounded, 'mmol/g')




    even_numbers = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 0]
    odd_numbers= [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 1]
    all_combinations_left  = [(value1, value2) for value1 in planes_left for value2 in even_numbers ]
    all_combinations_right = [(value1, value2) for value1 in planes_right  for value2 in odd_numbers ]
    random.shuffle(all_combinations_left)
    random.shuffle(all_combinations_right)
    sel_left_array = []
    sel_right_array = []

    while len(sel_left_array) < left_count_update:
        if not all_combinations_left:  
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_left = random.choice(all_combinations_left)
        sel_left_array.append(selected_combination_left)
        all_combinations_left.remove(selected_combination_left)
    
    while len(sel_right_array) < right_count_update:
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
    
    origin_glc_chains = molecule.load("pdb", "layer_after_temp.pdb")
    left_glc_candidates_indices   = get_fragment_indices(origin_glc_chains, sel_left_array)    
    right_glc_candidates_indices  = get_fragment_indices(origin_glc_chains, sel_right_array)      
    molecule.delete(origin_glc_chains) 


    def remove_atoms_based_on_indices(universe, index_ranges):
        all_removed_indices = []
    
        for index_range in index_ranges:
            for start_idx, end_idx in index_range:
                atom_group = universe.atoms[start_idx:end_idx + 1]
                for residue in atom_group.residues:
                    indices_to_remove = residue.atoms[17:21].indices
                    all_removed_indices.extend(indices_to_remove)
        unique_indices = sorted(set(all_removed_indices))
        indices_str = ' '.join(map(str, unique_indices))
        universe = universe.select_atoms(f"not index {indices_str}")

        return universe

    all_indices = [left_glc_candidates_indices, right_glc_candidates_indices]

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
    left_so3_candidates_indices   = get_fragment_indices(primary_hydroxyl_remove_chains, sel_left_array)    
    right_so3_candidates_indices  = get_fragment_indices(primary_hydroxyl_remove_chains, sel_right_array)    
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
                    base_atom_index = 16   
                    base_atom = residue.atoms[base_atom_index]
                    so3_new_positions = {}
                    # Define new positions based on residue ID's odd/even status
                    if fragment_number in planes_left:  # Odd residue ID
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([-0.632, -0.207, -0.838000000000001]),
                            'H62':  base_atom.position + np.array([-0.632,  0.151, 0.850000000000001]),
                            'O6':   base_atom.position + np.array([ 0.658,  1.242, -0.263999999999996]),
                            'S6':   base_atom.position + np.array([-0.449,  2.489, -0.177999999999997]),
                            'OS62': base_atom.position + np.array([-1.658,  2.087, 0.902000000000001]),
                            'OS63': base_atom.position + np.array([ 0.319,  3.879, 0.338999999999999]),
                            'OS64': base_atom.position + np.array([-1.114,  2.748, -1.687])
                        }
                    elif fragment_number in planes_right:  # Even residue ID
                        so3_new_positions = {
                            'H61':  base_atom.position + np.array([ 0.632,   0.207, -0.838000000000001]),
                            'H62':  base_atom.position + np.array([ 0.632,  -0.151, 0.849999999999994]),
                            'O6':   base_atom.position + np.array([-0.659,  -1.241, -0.264000000000003]),
                            'S6':   base_atom.position + np.array([ 0.474,  -2.465, -0.356000000000002]),
                            'OS62': base_atom.position + np.array([ 1.000,  -2.618, -1.934]),
                            'OS63': base_atom.position + np.array([ 1.771,  -2.099, 0.629999999999995]),
                            'OS64': base_atom.position + np.array([-0.217,  -3.902, 0.141999999999996])
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
    
    
    all_indices = [left_so3_candidates_indices, right_so3_candidates_indices]
    so3_negative_u = add_so3_based_on_indices(primary_hydroxyl_remove_chains_u, all_indices)
    so3_negative_u.atoms.write("so3_negative-modified.temp.pdb")
    molecule.delete(primary_hydroxyl_remove_chains)#
#

#
    u_final = mda.Universe("so3_negative-modified.temp.pdb")
    box_x = num_h * a_trans+50
    box_y = num_h * b_trans+50
    box_z = c_iterations * c_trans
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z,90,90,90]
    translation_vector = center_of_box - center_of_mass
    u_final.atoms.translate(translation_vector)
    with mda.Writer(f"charmm36-cellulose-Ibeta-para-icm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
        
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    
    print(f"Half-Ester sulfate content: {actual_sulfate_content:.2f}, Degree of sulfate: {actual_ds:.4f}")

elif sulfate_num == 0:
    u_final = mda.Universe("layer_after_temp.pdb")
    box_x = num_h * a_trans+50
    box_y = num_h * b_trans+50
    box_z = c_iterations * c_trans
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    with mda.Writer(f"charmm36-cellulose-Ibeta-para-icm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
    


    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    actual_ds=0.0    
    actual_sulfate_content=0.0
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"Half-Ester sulfate content: {actual_sulfate_content:.2f}, Degree of sulfate: {actual_ds:.4f}")

