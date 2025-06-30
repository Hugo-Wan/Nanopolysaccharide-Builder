import MDAnalysis as mda
from MDAnalysis.coordinates import PDB
from MDAnalysis.core.universe import Merge
import warnings
import os
import glob
import json
import sys
import re
import numpy as np
warnings.filterwarnings("ignore", category=UserWarning)



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
    c_iterations = int(sys.argv[4])
    if c_iterations < 0:
        sys.stderr.write("Error: Please provide a positive integer for c repetition number.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length value provide. Please enter a valid numeric value.\n")
    sys.exit(1)
    
try:
    height = float(sys.argv[5])
    a_iterations = int(height // a_trans)
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



with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
neutron_pdb_1 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'glycam06', 'AB_configuration', 'left-unit.pdb')
user_defined_pdb_1 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'glycam06', 'AB_configuration', 'left-unit_ud.pdb')

if os.path.exists(user_defined_pdb_1):
    unit_chain_input_file_1 = user_defined_pdb_1
else:
    unit_chain_input_file_1 = neutron_pdb_1

neutron_pdb_2 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'glycam06', 'AB_configuration', 'right-unit.pdb')
user_defined_pdb_2 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'glycam06', 'AB_configuration', 'right-unit_ud.pdb')

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

box_x = a_trans * a_iterations + 20
box_y = b_iterations * b_trans / 2 + 20
box_z = c_iterations * c_trans 
center_of_mass = final_structure.atoms.center_of_mass()
center_of_box = [box_x / 2, box_y / 2, box_z / 2]
box_dimensions = [box_x, box_y, box_z, 90, 90, 90]
final_structure.dimensions = box_dimensions
box_translation_vector = center_of_box - center_of_mass
final_structure.atoms.translate(box_translation_vector)
with mda.Writer("glycam06-alpha-chitin-AB-hexagon-icm.pdb", n_atoms=final_structure.atoms.n_atoms, reindex=True) as W:
    W.write(final_structure.atoms)

# Cleanup any remaining temporary files
for temp_file in glob.glob("*temp*.pdb"):
    os.remove(temp_file)


#print("Generated hexagonal-shape alpha-chitin-A with infinite-chain model.")