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



a_trans = 7.784
b_trans = 8.201
c_trans = 10.38
gamma_angle=96.5


try:
    c_iterations = int(sys.argv[1])
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
    height = float(sys.argv[2])
    num_h = int(round(height /7.8))
    if num_h < 1:
        sys.stderr.write("Error: Please provide an larger height value.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid height value provide. Please enter a valid numeric value.\n")
    sys.exit(1)    


try:
    width = float(sys.argv[3])
    num_w = int((round(width/5.3)/2))
    if num_w < 1:
        sys.stderr.write("Error: Please provide an larger width value.\n")
        sys.exit(1)


        sys.exit(1)        
except ValueError:
    sys.stderr.write("Error: Invalid width value provide. Please enter a valid numeric value.\n")
    sys.exit(1)    


try:
    Tempo_target = float(sys.argv[4])  #####
    #if not (0 <= Tempo_target < 1):
    #    sys.stderr.write("Error: Tempo-oxidation degree  must be between 0 and 1.\n") 
    #    sys.exit(1)
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

#DDA_target=0.20      




with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
chain_1 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'charmm36', 'chain-1.pdb')
chain_2 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'charmm36', 'chain-2.pdb')

chain_1_u = mda.Universe(chain_1)
chain_2_u = mda.Universe(chain_2)

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


 ##terminal side            
chain_1_finite_build_chain = "chain_1_temp.pdb"
chain_1_finite_u = mda.Universe(chain_1_finite_build_chain)
chain_1_first_atom = chain_1_finite_u.atoms[0]
chain_1_o4_atoms = chain_1_finite_u.select_atoms("name O4")
if not chain_1_o4_atoms:
    raise ValueError("No atoms named 'O4' found in the structure.")
chain_1_last_o4 = chain_1_o4_atoms[-1]

chain_1_new_positions = {
    'O1': chain_1_first_atom.position + [0.538, -0.449, 1.247],
    'HO1': chain_1_first_atom.position + [0.129, -0.365, 1.932],
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


 ##terminal side            
chain_2_finite_build_chain = "chain_2_temp.pdb"
chain_2_finite_u = mda.Universe(chain_2_finite_build_chain)
chain_2_first_atom = chain_2_finite_u.atoms[0]
chain_2_o4_atoms = chain_2_finite_u.select_atoms("name O4")
if not chain_2_o4_atoms:
    raise ValueError("No atoms named 'O4' found in the structure.")
chain_2_last_o4 = chain_2_o4_atoms[-1]

chain_2_new_positions = {
    'O1': chain_2_first_atom.position + [0.538, -0.449, 1.247],
    'HO1': chain_2_first_atom.position + [0.129, -0.365, 1.932],
    'HO4': chain_2_last_o4.position + [0.377, -0.192, -0.892]
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
#------------------------------------------------------------------------------
###translate parameter
tilt_angle_rad= 6.5 * math.pi/180
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




####modified coo-
pKa=3.25                        # pka of coo-
m_Glc=163.09316                 # Glc unit mass
m_difference=111                # difference between Glc unit with Glc-cooNa(3)




##crystallographic planes defined



planes_left = []


##-------------------crystallographic-plane-definitions---------------------
for m in range(1, num_h + 1):
    if m == 1:
        value1 = num_w
        value2 = num_w - 1
        planes_left.extend([value2])
    elif m>1 and m<num_h:
        value1 = num_w + 2 * (m - 1) * num_w
        value2 = (num_w - 1) + 2 * (m - 1) * num_w
        value3 = num_w + 2 * (m - 1) * num_w + num_w - 1
        value4 = num_w * 2 * (m - 1)
        planes_left.extend([value2, value3])
    elif m==num_h:
        value1 = num_w + 2 * (m - 1) * num_w
        value2 = (num_w - 1) + 2 * (m - 1) * num_w
        value_2_1=value2+num_w
        planes_left.extend([value2, value_2_1])

left_right_count=len(planes_left)



#--------------coo modified step ------------------------------------------
chain_u = mda.Universe("layer_after_temp.pdb")
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


if carboxylation_num > 0:
    left_count=carboxylation_num
    planes_left_count   =  len(planes_left)



    max_count_left = c_iterations * planes_left_count 


    left_count_update = left_count

    if left_count > max_count_left:
        left_count_update = max_count_left

        
    new_carboxylation_number_update=left_count_update  
    actual_ds=(new_carboxylation_number_update)/(total_residues)
    actual_carboxylate_content=Tempo_target*(actual_ds/ds)
    actual_carboxylate_rounded=round(actual_carboxylate_content,4)




    even_numbers = [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 0]
    odd_numbers= [x for x in range(1, 2 * c_iterations + 1) if x % 2 == 1]
    all_combinations_left  = [(value1, value2) for value1 in planes_left for value2 in even_numbers ]
    random.shuffle(all_combinations_left)
    sel_left_array = []

    while len(sel_left_array) < left_count_update:
        if not all_combinations_left:  
            raise ValueError("Ran out of unique combinations before reaching the desired count.")
        selected_combination_left = random.choice(all_combinations_left)
        sel_left_array.append(selected_combination_left)
        all_combinations_left.remove(selected_combination_left)
    

    #print("100 bottom array residue:", sel_100_bottom_array)



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

    all_indices = [left_glc_candidates_indices]

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
    left_coo_candidates_indices   = get_fragment_indices(primary_hydroxyl_remove_chains, sel_left_array)    
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
                       base_atom_index = 18   
                    elif residue.resid == maximum_resid_num :
                       base_atom_index = 17   
                    else :
                       base_atom_index = 16   
                    base_atom = residue.atoms[base_atom_index]
                    coo_new_positions = {} 

                    if residue.resid % 2 == 0 and fragment_number in planes_left:  # Odd residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 0.658,  1.242,   -0.264]),
                            'O62': base_atom.position + np.array([-1.429, -0.063,    0.013])
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
    
    
    all_indices = [left_coo_candidates_indices]
    nh3_u = add_coo_based_on_indices(primary_hydroxyl_remove_chains_u, all_indices)
    nh3_u.atoms.write("coo_negative-modified.temp.pdb")
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
    left_coo_netural_candidates_indices   = get_fragment_indices(coo_deprotonation_chains, sel_left_array)    

    coo_netural_indices_groups = [left_coo_netural_candidates_indices]
    coo_negative_chains_u = mda.Universe("coo_negative-modified.temp.pdb")
    coo_negative_to_remove = collect_carboxylation_atoms_to_remove_randomly(coo_negative_chains_u, coo_netural_indices_groups, coo_netural)
    
    coo_netural_apply_modifications_and_save(coo_negative_chains_u, coo_negative_to_remove, "coo_netural-modified.temp.pdb")
    
    molecule.delete(coo_deprotonation_chains) 
    
    u_final = mda.Universe("coo_netural-modified.temp.pdb")
    box_x = num_h * a_trans+50
    box_y = num_h * b_trans+50
    box_z = c_iterations * c_trans+20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    u_final.atoms.translate(translation_vector)
    with mda.Writer(f"charmm36-cellulose-Ibeta-para-fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
        
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    
    print(f"Carboxylate content: {actual_carboxylate_content:.2f}, Degree of carboxylation: {actual_ds:.4f}, pH: {actual_pH_rounded:.4f}")


elif carboxylation_num == 0:
    u_final = mda.Universe("layer_after_temp.pdb")
    box_x = num_h * a_trans+50
    box_y = num_h * b_trans+50
    box_z = c_iterations * c_trans+20
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    with mda.Writer(f"charmm36-cellulose-Ibeta-para-fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
    

    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    actual_ds=0.0    
    actual_carboxylate_content=0.0
    actual_pH_rounded=0.0
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"Carboxylate content: {actual_carboxylate_content:.2f}, Degree of carboxylation: {actual_ds:.4f}, pH: {actual_pH_rounded:.4f}")

