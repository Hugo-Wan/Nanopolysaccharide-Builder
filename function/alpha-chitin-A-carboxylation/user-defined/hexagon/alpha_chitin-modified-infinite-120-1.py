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
    if a_trans < 4 or b_trans < 15 or c_trans < 10:
        sys.stderr.write("Error: Please provide an effective values for all crystallographic parameters.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid crystallographic value provided. Please enter valid numeric values.\n")
    sys.exit(1)



#-------------------------------------------- modified parameter  --------------------------------------------------


try:
    c_iterations = int(sys.argv[4])
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
    height = float(sys.argv[5])
    a_iterations = int(height // a_trans)
    #print(a_iterations)
    if a_iterations < 1:
        sys.stderr.write("Error: Please provide an larger height value.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid height value provide. Please enter a valid numeric value.\n")
    sys.exit(1)    


try:
    width = float(sys.argv[6])
    b_iterations = int(round(width / (b_trans * 0.5)))
    if b_iterations % 2 == 0:  
       b_iterations -= 1   

    if b_iterations < 1:
        sys.stderr.write("Error: Please provide an larger width value.\n")
        sys.exit(1)


    mid_value = (b_iterations - 1) // 2 
    #print(mid_value)
    #print(a_iterations)
    extra_check_for_a=a_iterations-mid_value
    #print(extra_check_for_a)
    if extra_check_for_a <= 0:
        sys.stderr.write("Error: Please decrease the width value or increase the height value.\n")

        sys.exit(1)        
except ValueError:
    sys.stderr.write("Error: Invalid width value provide. Please enter a valid numeric value.\n")
    sys.exit(1)    


try:
        ##carboxylation degree ref:/10.1016/j.carbpol.2019.115292
        Tempo_target = float(sys.argv[7])  #####

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

pKa=3.25                               # pka of coo-
m_GlcNAc  = 204.20                     # GlcNAc unit mass
m_difference = 36.973769               # difference between GlcNAc unit with GlcNAc-cooNa
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

left_temp = "left-chain_temp.pdb"
left_temp_u=mda.Universe(left_temp)

for atom in left_temp_u.atoms:
    atom.segment.segid = '0'  
left_temp_u.atoms.write("left-chain_temp.2.pdb")



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
 # ------------------------------------------------   
 #right_chain assembly    




#-----------------------------------------------------------------------------
#center_layer assembly


#center_layer_assembly
center_layer_input_file = "left-chain_temp.2.pdb"
layer_r_u = mda.Universe(center_layer_input_file)
center_layer = []  
for i in range(1, a_iterations + 1):
    layer_r_u .atoms.positions += [a_trans, 0, 0]
    #j=i+a_iterations
    segid_increment_value = 1
    for atom in layer_r_u.atoms:
        numeric_part = atom.segid.replace('','')
        new_numeric_part = int(numeric_part) + segid_increment_value
    new_segid = f"{new_numeric_part}"
    atom.segment.segid = new_segid

    center_layer_output = f"center_layer_{i}.pdb"
    center_layer.append(center_layer_output)
    with mda.Writer(center_layer_output, n_atoms=layer_r_u .atoms.n_atoms) as W:
        W.write(layer_r_u .atoms)
with open(f"center_layer_temp.pdb", "w") as layer_file:
    for center_layer_output in center_layer:
        with open(center_layer_output, "r") as center_layer_pdb_file:
            for center_layer_line in center_layer_pdb_file:
                if center_layer_line.startswith("ATOM"):
                    layer_file.write(center_layer_line)
            os.remove(center_layer_output) 
        

 # Use integer division

for m in range(1, mid_value+1):
    if m % 2 == 1:   #odd number
        right_layer_input_file_1 = "right-chain_temp.2.pdb"
        right_initial_trans_vector = np.array([a_trans * (m - 1)/2, -b_trans * (m - 1)/2, 0])
        right_input_u = mda.Universe(right_layer_input_file_1)
        right_input_u.atoms.positions += right_initial_trans_vector
        right_input_u.atoms.write(f"right_initial_{m}.pdb")
        layer_initial_input=f"right_initial_{m}.pdb"
        layer_r_u = mda.Universe(layer_initial_input)
        right_layer = []  
        for i in range(1, (a_iterations-m+1)):
            layer_r_u.atoms.positions += [a_trans, 0, 0]
            segid_increment_value = 1
            for atom in layer_r_u.atoms:
                numeric_part = atom.segid.replace('','')
                new_numeric_part = int(numeric_part) + segid_increment_value
            new_segid = f"{new_numeric_part}"
            atom.segment.segid = new_segid

            right_layer_output = f"right_layer_{i}.pdb"
            right_layer.append(right_layer_output)
            with mda.Writer(right_layer_output, n_atoms=layer_r_u.atoms.n_atoms) as W:
                W.write(layer_r_u.atoms)
        with open(f"right_temp_{m}.pdb", "w") as layer_file:
            for right_layer_output in right_layer:
                with open(right_layer_output, "r") as right_layer_pdb_file:
                    for right_layer_line in right_layer_pdb_file:
                        if right_layer_line.startswith("ATOM"):
                            layer_file.write(right_layer_line)
                    os.remove(right_layer_output) 
        os.remove(layer_initial_input)


        right_output_file=f"right_temp_{m}.pdb"
        right_output_u=mda.Universe(right_output_file)
        for seg in right_output_u.segments:  
            if m == 1:
                segid_increment_value = 0
            elif m > 1:
                segid_increment_value = 0
                for t in range(1, m):
                    segid_increment_value += (a_iterations - t) * 2
            new_id = str(int(seg.segid) + segid_increment_value)  
            #print(new_id)
            for atom in seg.atoms:
                atom.segment.segid = new_id
        right_output_new =f"right_updated_temp_{m}.pdb"
        with mda.Writer(right_output_new, n_atoms=right_output_u.atoms.n_atoms, reindex=True) as W:
            W.write(right_output_u.atoms) 


        ###another side the same structure   
        left_input_u=mda.Universe(right_output_new)
        left_duplicate = left_input_u.copy()
        left_translation_vector_1 = [0, m*b_trans, 0]
        left_duplicate.atoms.translate(left_translation_vector_1)

        for seg_left in left_duplicate.segments:  
            segid_increment_value_left = a_iterations - m
            original_segid = seg_left.segid
            #print(segid_increment_value_left)
            #print(original_segid)
            new_id_left = str(int(seg_left.segid) + segid_increment_value_left)  
            #print(new_id_left)
            for atom_left in seg_left.atoms:
                atom_left.segment.segid = new_id_left
        left_output =f"left_updated_temp_{m}.pdb"
        with mda.Writer(left_output, n_atoms=left_duplicate.atoms.n_atoms, reindex=True) as W:
            W.write(left_duplicate.atoms)
            
    if m % 2 == 0:   #even number
        right_layer_input_file_1 = "left-chain_temp.2.pdb"
        right_initial_trans_vector = np.array([a_trans * (m)/2, -b_trans * (m)/2, 0])
        right_input_u = mda.Universe(right_layer_input_file_1)
        right_input_u.atoms.positions += right_initial_trans_vector
        right_input_u.atoms.write(f"right_initial_{m}.pdb")
        layer_initial_input=f"right_initial_{m}.pdb"
        layer_r_u = mda.Universe(layer_initial_input)
        right_layer = []  
        for i in range(1, (a_iterations-m+1)):
            layer_r_u.atoms.positions += [a_trans, 0, 0]
            segid_increment_value = 1
            for atom in layer_r_u.atoms:
                numeric_part = atom.segid.replace('','')
                new_numeric_part = int(numeric_part) + segid_increment_value
            new_segid = f"{new_numeric_part}"
            atom.segment.segid = new_segid

            right_layer_output = f"right_layer_{i}.pdb"
            right_layer.append(right_layer_output)
            with mda.Writer(right_layer_output, n_atoms=layer_r_u.atoms.n_atoms) as W:
                W.write(layer_r_u.atoms)
        with open(f"right_temp_{m}.pdb", "w") as layer_file:
            for right_layer_output in right_layer:
                with open(right_layer_output, "r") as right_layer_pdb_file:
                    for right_layer_line in right_layer_pdb_file:
                        if right_layer_line.startswith("ATOM"):
                            layer_file.write(right_layer_line)
                    os.remove(right_layer_output) 
        os.remove(layer_initial_input)

        right_output_file=f"right_temp_{m}.pdb"
        right_output_u=mda.Universe(right_output_file)
        n=m-1
        #print(n)
        left_prev_file = f"left_updated_temp_{n}.pdb"
        left_prev_u = mda.Universe(left_prev_file)
        last_seg_id = int(left_prev_u.segments[-1].segid)
        #print(last_seg_id)
        for seg in right_output_u.segments:  
            segid_increment_value = last_seg_id 
            new_id = str(int(seg.segid) + segid_increment_value)  
            for atom in seg.atoms:
                atom.segment.segid = new_id
        right_output_new =f"right_updated_temp_{m}.pdb"
        with mda.Writer(right_output_new, n_atoms=right_output_u.atoms.n_atoms, reindex=True) as W:
            W.write(right_output_u.atoms) 

        ###another side the same structure    
        left_input_u=mda.Universe(right_output_new)
        left_duplicate = left_input_u.copy()

        left_translation_vector_2 = [0, m*b_trans, 0]
        left_duplicate.atoms.translate(left_translation_vector_2)
        for seg in left_duplicate.segments:  
            segid_increment_value = (a_iterations-m)
            new_id = str(int(seg.segid) + segid_increment_value)  
            for atom in seg.atoms:
                atom.segment.segid = new_id
        left_output =f"left_updated_temp_{m}.pdb"
        with mda.Writer(left_output, n_atoms=left_duplicate.atoms.n_atoms) as W:
            W.write(left_duplicate.atoms)


center_structure_file = f"center_layer_temp.pdb"
center_u = mda.Universe(center_structure_file)

crystal_structures = [center_u]

for w in range(1, mid_value + 1):
    first_structure =f"right_updated_temp_{w}.pdb"   # Use w instead of m
    second_structure = f"left_updated_temp_{w}.pdb" 
    first_u = mda.Universe(first_structure)
    second_u = mda.Universe(second_structure)
    combined_u = mda.Merge(first_u.atoms, second_u.atoms)
    crystal_structures.append(combined_u)

if not crystal_structures:
    raise ValueError("No structures were loaded or combined. Please check input files and parameters.")

final_structure = crystal_structures[0]

for universe in crystal_structures[1:]:
    final_structure = mda.Merge(final_structure.atoms, universe.atoms)

with mda.Writer("alpha-chitin-A-temp.pdb", n_atoms=final_structure.atoms.n_atoms, reindex=True) as W:
    W.write(final_structure.atoms)


#-------------------------------------------------- native chitin built ----------------------------------------------




#------------------------------------------------- coo modified step -------------------------------------------------------
####120 plane selection
###fragment 
# m=0
# (0, a_iteration-1), 
# m=1
# (a_iteration, a_iteration+a_iterations-(m+1)), 
# (a_iteration+ a_iterations-(m+1) + 1, a_iteration+ a_iterations-(1m+1) + 1 + a_iterations-(m+1))
# m=2
# (2 * a_iteration+1 + a_iterations-1 +1 ,                    2 * a_iteration+1 + a_iterations-1 +  a_iterations-1 ), 
# (2 * a_iteration+1 + a_iterations-1 + + a_iterations-1 +1 , 2 * a_iteration+1 + a_iterations-1 +  a_iterations-1 +1 + a_iterations-1 )
# m=3
# (3 * a_iteration+1 + a_iterations-1 +1 ,                    2 * a_iteration+1 + a_iterations-1 +  a_iterations-1 ), 
# ....
# m=m
# (a_iteration+1, a_iteration+1+a_iterations-1)
## (2m-1)*a_iterations-((m-1)**2+(m-1)) ;  2*m*a_iterations-(m**2+1) ; 2*m*a_iterations-m**2 ; (2*m+1)*a_iterations-((m+1)**2-m+1) ;  

####010 plane selection
# (2m-1)*a_iterations-((m-1)**2+(m-1)) to 2*m*a_iterations-(m**2+1);  2*m*a_iterations-m**2 to (2*m+1)*a_iterations-((m+1)**2-m+1)
# 
# 
# 

##crystallographic planes defined
planes_010_left = []
planes_010_right = []

planes_120_even = []
planes_120_odd = []
###for region categorized
planes_120_upper = []
planes_120_bottom = []
planes_120_bottom.extend([0])
planes_120_upper.extend([a_iterations-1])
planes_120_even.extend([0, a_iterations-1])

for m in range(1, mid_value):
    # Calculate values for 120_planes
    value1 = (2*m - 1) * a_iterations - ((m - 1)**2 + (m - 1))
    value2 = 2 * m * a_iterations - (m**2 + 1)
    value3 = 2 * m * a_iterations - m**2
    value4 = (2*m + 1) * a_iterations - ((m + 1)**2 - m)
    planes_120_bottom.extend([value1, value3])
    planes_120_upper.extend([value2, value4]) 
    if m % 2 == 0:
        planes_120_even.extend([value1, value2, value3, value4])
    else:
        planes_120_odd.extend([value1, value2, value3, value4])
    # Add calculated values to the 120_planes list
 

# Calculate ranges for 010_planes when m = mid_value
range1_start = (2*mid_value - 1) * a_iterations - ((mid_value - 1)**2 + (mid_value - 1))
range1_end = 2 * mid_value * a_iterations - (mid_value**2 + 1)
range2_start = 2 * mid_value * a_iterations - mid_value**2
range2_end = (2*mid_value + 1) * a_iterations - ((mid_value + 1)**2 - mid_value)

# Add ranges as tuples to the 010_planes list
planes_010_left.extend([range2_start, range2_end])
planes_010_right.extend([range1_start, range1_end])
planes_010_count=(planes_010_left[-1] - planes_010_left[0] + 1) + (planes_010_right[-1] - planes_010_right[0] + 1)
planes_120_count=len(planes_120_even)+len(planes_120_odd)

modified_chains_count=planes_010_count+planes_120_count

#print(planes_120_count)
#print("120_planes_even:", planes_120_even)
#print("120_planes_odd:", planes_120_odd)
#print("120_planes_upper:", planes_120_upper)
#print("120_planes_bottom:", planes_120_bottom)
#
#print("010_planes_r:", planes_010_right)
#print("010_planes_l:", planes_010_left)


chain_u = mda.Universe("alpha-chitin-A-temp.pdb")
residue_ids = [residue.resid for residue in chain_u.residues]
total_residues=len(residue_ids)
maximum_resid_num=residue_ids[-1]
#origin_nacetyl_chains = molecule.load("pdb", "alpha-chitin-A-temp.pdb")
#native_all_atoms = atomsel('all', molid=origin_nacetyl_chains)
#residue_ids = native_all_atoms.resid
#min_resid = min(residue_ids, default=None)
#max_resid = max(residue_ids, default=None)
#print(f"Minimum residue ID: {min_resid}")
#print(f"Maximum residue ID: {max_resid}")
#print("Left Modified Chains:", left_modified_chains)
#print("Right Modified Chains:", right_modified_chains)
ds_reverse= (1/(m_GlcNAc*Tempo_target*(10**-3))-(m_difference/m_GlcNAc) )  ##(1/ds)
ds=1/ds_reverse
n_coo=round(total_residues*ds)

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


    planes_010_left_count=(planes_010_left[-1] - planes_010_left[0] + 1)
    planes_010_right_count=(planes_010_right[-1] - planes_010_right[0] + 1)
    planes_120_upper_count=len(planes_120_upper)
    planes_120_bottom_count=len(planes_120_bottom)

    max_010_count_left = c_iterations * planes_010_left_count 
    max_010_count_right = c_iterations * planes_010_right_count

    max_120_count_upper = c_iterations * 2 * planes_120_upper_count
    max_120_count_bottom = c_iterations * 2 * planes_120_bottom_count

    upper_120_count_update = upper_120_count

    if upper_120_count > max_120_count_upper:
        upper_120_count_update = max_120_count_upper

    new_carboxylation_number_update=upper_120_count_update 
    actual_ds=(new_carboxylation_number_update)/(total_residues)
    actual_carboxylate_content=Tempo_target*(actual_ds/ds)

    
    
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

    numbers_120 = [x for x in range(1, 2 * c_iterations + 1)]
    
    all_combinations_120_upper = [(value1, value2) for value1 in planes_120_upper  for value2 in numbers_120 ]
    #print(all_combinations_120_upper)
    all_combinations_120_bottom = [(value1, value2) for value1 in planes_120_bottom  for value2 in numbers_120 ]
    #print(all_combinations_120_bottom)
    
    
    #print(len(all_combinations_010_left))
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
    
    origin_nacetyl_chains = molecule.load("pdb", "alpha-chitin-A-temp.pdb")
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
                return unique_fragments.pop()  # Return the single unique fragment number
        else:
            raise ValueError("Multiple unique fragments found in selection")
    
    primary_hydroxyl_remove_chains = molecule.load("pdb", "primary_hydroxyl_remove.temp.pdb")
    upper_coo_candidates_indices  = get_fragment_indices(primary_hydroxyl_remove_chains, sel_120_upper_array)   
 

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
                    base_atom_index = 22   
                    base_atom = residue.atoms[base_atom_index]
                    coo_new_positions = {} 
                    # Define new positions based on residue ID's odd/even status
                    if residue.resid % 2 == 1  and fragment_number in planes_120_even:  # Odd residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([1.367,   0.411, 0.081]),
                            'O62': base_atom.position + np.array([-1.046,  0.957, 0.187])
                        }
                    elif residue.resid % 2 == 0  and fragment_number in planes_120_even:  # Even residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([-1.367, -0.41,  0.080]),
                            'O62': base_atom.position + np.array([1.046,  -0.957, 0.188])
                        }
    
                    elif residue.resid % 2 == 1  and fragment_number in planes_120_odd:  # Odd residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.047, 0.956, -0.187]),
                            'O62': base_atom.position + np.array([-1.367, 0.411, -0.081])                            

                        }
                    elif residue.resid % 2 == 0  and fragment_number in planes_120_odd:  # Even residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([-1.046,  -0.957, -0.188]),
                            'O62': base_atom.position + np.array([-1.367,  -0.411,  0.081])                            
                        }
    
                    elif fragment_number in range(min(planes_010_left), max(planes_010_left)+1):  # Odd residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.047,  0.956,  -0.187]),
                            'O62': base_atom.position + np.array([-1.367,   0.411,  -0.081])
                        }
                    elif fragment_number in range(min(planes_010_right), max(planes_010_right)+1):  # Even residue ID
                        coo_new_positions = {
                            'O61': base_atom.position + np.array([ 1.367,  -0.41,   0.080]),
                            'O62': base_atom.position + np.array([-1.046,  -0.957, -0.188])
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
    
    
    all_indices = [upper_coo_candidates_indices]
    coo_negative_u = add_coo_based_on_indices(primary_hydroxyl_remove_chains_u, all_indices)
    coo_negative_u.atoms.write("coo_negative-modified.temp.pdb")
    molecule.delete(primary_hydroxyl_remove_chains)
    #    
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
    box_x = a_trans * a_iterations
    box_y = b_iterations * b_trans
    box_z = c_iterations * c_trans
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    u_final.atoms.translate(translation_vector)
    with mda.Writer(f"charmm36-alpha-chitin-A-hexagon-icm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
        
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    
    print(f"Carboxylate content: {actual_carboxylate_content:.2f}, Degree of carboxylation: {actual_ds:.4f}, pH: {actual_pH_rounded:.4f}")

elif carboxylation_num == 0:
    u_final = mda.Universe("alpha-chitin-A-temp.pdb")
    box_x = a_trans * a_iterations 
    box_y = b_iterations * b_trans 
    box_z = c_iterations * c_trans 
    center_of_mass = u_final.atoms.center_of_mass()
    center_of_box = [box_x / 2, box_y / 2, box_z / 2]
    u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
    translation_vector = center_of_box - center_of_mass
    with mda.Writer(f"charmm36-alpha-chitin-A-hexagon-icm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True, ) as W:
        W.write(u_final.atoms)
    

    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    actual_ds=0.0    
    actual_carboxylate_content=0.0
    actual_pH_rounded=0.0
    #print("Generated alpha-chitin-A fcm model with surface modification.")
    print(f"Carboxylate content: {actual_carboxylate_content:.2f}, Degree of carboxylation: {actual_ds:.4f}, pH: {actual_pH_rounded:.4f}")

