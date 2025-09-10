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
    gamma_angle= float(sys.argv[4])
    if a_trans <= 0 or b_trans <= 0 or c_trans <= 0:
        sys.stderr.write("Error: Please provide positive values for all crystallographic parameters.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length value provided. Please enter valid numeric values.\n")
    sys.exit(1)



try:
    c_iterations = int(sys.argv[5])
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
    height = float(sys.argv[6])
    a_iterations = int(height // a_trans)
    if a_iterations < 1:
        sys.stderr.write("Error: Please provide an larger height value.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid width value provide. Please enter a valid numeric value.\n")
    sys.exit(1)    

try:
    width = float(sys.argv[7])
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
        ##carboxylation degree ref:/10.1016/j.carbpol.2019.115292
        Tempo_target = float(sys.argv[8])  #####
        #if not (0 <= Tempo_target < 1):
        #    sys.stderr.write("Error: Tempo-oxidation degree  must be between 0 and 1.\n") 
        #    sys.exit(1)

        pH = float(sys.argv[9])
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
#DDA_target=0.20      
# 
#          

pKa=3.25                             #pka of coo-
m_GlcNAc  = 204.20                   # GlcNAc unit mass
m_difference = 36.973769               # difference between GlcNAc unit with GlcNAc-cooNa
#-------------------------------------------- modified parameter  --------------------------

with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
neutron_pdb = os.path.join(main_folder_path, 'structure', 'beta_chitin', 'charmm36', 'unit-neutron.pdb')
user_defined_pdb = os.path.join(main_folder_path, 'structure', 'beta_chitin', 'charmm36', 'unit-ud.pdb')

if os.path.exists(user_defined_pdb):
    unit_chain_input_file = user_defined_pdb
else:
    unit_chain_input_file = neutron_pdb

unit_u = mda.Universe(unit_chain_input_file)


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

chain_u = mda.Universe("unit_temp.pdb")
residue_ids = [residue.resid for residue in chain_u.residues]
maximum_resid_num=residue_ids[-1]
total_residues=len(residue_ids)
#

ds_reverse= (1/(m_GlcNAc*Tempo_target*(10**-3))-(m_difference/m_GlcNAc) )  ##(1/ds)
ds=1/ds_reverse
n_coo=round(total_residues*ds)
#print('degree of substitution', ds)
#print('carboxylation_number',n_coo)
#
#
def adjust_for_carboxylation_limit(initial_carboxylation_num):
    initial_carboxylation_num = n_coo
    carboxylation_number_limit = c_iterations * planes_010_count + 2 * c_iterations * planes_120_count
    new_carboxylation_number = initial_carboxylation_num
    if initial_carboxylation_num > carboxylation_number_limit:
        new_carboxylation_number = carboxylation_number_limit
    return new_carboxylation_number

carboxylation_num = adjust_for_carboxylation_limit(n_coo)




if carboxylation_num > 0:

    upper_120_count=carboxylation_num
    #print(upper_120_count)
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


    upper_120_count_update = upper_120_count
    #print(upper_120_count_update)
    if upper_120_count > max_120_count_upper:
        upper_120_count_update = max_120_count_upper
    #print(upper_120_count_update)
        
    new_carboxylation_number_update=upper_120_count_update 
    #print(new_carboxylation_number_update)
    actual_ds=(new_carboxylation_number_update)/(total_residues)
    actual_carboxylate_content=Tempo_target*(actual_ds/ds)


    actual_carboxylate_rounded=round(actual_carboxylate_content,4)
    #print('surface charge density', actual_carboxylate_rounded, 'mmol/g')


   
    numbers_120 = [x for x in range(1, 2 * c_iterations + 1)] 
    all_combinations_120_upper = [(value1, value2) for value1 in planes_120_upper  for value2 in numbers_120 ]
    random.shuffle(all_combinations_120_upper)
    sel_120_upper_array = []

        
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
    upper_nacetyl_candidates_indices  = get_fragment_indices(origin_nacetyl_chains, sel_120_upper_array)   

    molecule.delete(origin_nacetyl_chains) 



    
    def remove_atoms_based_on_indices(universe, index_ranges):
        all_removed_indices = []
    
        for index_range in index_ranges:
            for start_idx, end_idx in index_range:
                atom_group = universe.atoms[start_idx:end_idx + 1]
                for residue in atom_group.residues:
                    indices_to_remove = residue.atoms[23:27].indices
                    all_removed_indices.extend(indices_to_remove)
        unique_indices = sorted(set(all_removed_indices))
        indices_str = ' '.join(map(str, unique_indices))
        universe = universe.select_atoms(f"not index {indices_str}")
        return universe

    
    all_indices = [upper_nacetyl_candidates_indices]
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
    upper_coo_candidates_indices  = get_fragment_indices(primary_hydroxyl_remove_chains, sel_120_upper_array)   
    primary_hydroxyl_remove_chains_u = mda.Universe("primary_hydroxyl_remove.temp.pdb")
    coo_atoms = list(primary_hydroxyl_remove_chains_u.atoms)#
#

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
                    base_atom_index = 22   
                    base_atom = residue.atoms[base_atom_index]
                    coo_new_positions = {} 
                    # Define new positions based on residue ID's odd/even status
                    if residue.resid % 2 == 1  and fragment_number in planes_120_upper:  # Odd residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.345,  -0.473,  -0.114]),
                            'O62': base_atom.position + np.array([-1.093,  -0.896,  -0.216])
                        }
                    elif residue.resid % 2 == 0  and fragment_number in planes_120_upper:  # Even residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.093,   0.896,  -0.216]),
                            'O62': base_atom.position + np.array([-1.345,   0.473,  -0.114])
                        }
    
                    elif residue.resid % 2 == 1  and fragment_number in planes_120_bottom:  # Odd residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.345,  -0.473,  -0.114]),
                            'O62': base_atom.position + np.array([-1.093,  -0.896,  -0.216])
                        }
                    elif residue.resid % 2 == 0  and fragment_number in planes_120_bottom:  # Even residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.093,   0.896,  -0.216]),
                            'O62': base_atom.position + np.array([-1.345,   0.473,  -0.114])
                        }
    
                    elif fragment_number in range(min(planes_010_left), max(planes_010_left)+1):  # Odd residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.093,  0.896,   -0.216]),
                            'O62': base_atom.position + np.array([-1.345,  0.473,   -0.114])
                        }
                    elif fragment_number in range(min(planes_010_right), max(planes_010_right)+1):  # Even residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([-1.093,  -0.896,  -0.216]),
                            'O62': base_atom.position + np.array([ 1.345,  -0.473,  -0.114])
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
        # Create a new universe with all atoms including the added ones
        new_universe = mda.Merge(*[mda.AtomGroup([atom]) for atom in coo_atoms])
        return new_universe
    
    
    all_indices = [upper_coo_candidates_indices]
    nh3_u = add_coo_based_on_indices(primary_hydroxyl_remove_chains_u, all_indices)
    nh3_u.atoms.write("coo_negative-modified.temp.pdb")
    molecule.delete(primary_hydroxyl_remove_chains)#

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
    upper_coo_netural_candidates_indices  = get_fragment_indices(coo_deprotonation_chains, sel_120_upper_array)   


    coo_netural_indices_groups = [upper_coo_netural_candidates_indices]
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
    with mda.Writer(f"charmm36-beta-chitin-para-icm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
        
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    
    print(f"Carboxylate content: {actual_carboxylate_content:.2f}, Degree of carboxylation: {actual_ds:.4f}, pH: {actual_pH_rounded:.4f}")


elif carboxylation_num == 0:
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
    actual_ds=0.0    
    actual_carboxylate_content=0.0
    actual_pH_rounded=0.0
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"Carboxylate content: {actual_carboxylate_content:.2f}, Degree of carboxylation: {actual_ds:.4f}, pH: {actual_pH_rounded:.4f}")