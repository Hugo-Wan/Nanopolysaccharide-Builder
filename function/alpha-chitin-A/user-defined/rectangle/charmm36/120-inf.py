import MDAnalysis as mda
from MDAnalysis.coordinates import PDB
from MDAnalysis.core.universe import Merge
import warnings
import os
import glob
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
    sys.stderr.write("Error: Invalid length value provided. Please enter valid numeric values.\n")
    sys.exit(1)

try:
    a_iterations = int(sys.argv[4])
    b_iterations = int(sys.argv[5])
    c_iterations = int(sys.argv[6])
    if a_iterations <= 0 or b_iterations <= 0 or c_iterations <= 0:
        sys.stderr.write("Error: Please provide positive integers for all repetition numbers.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length value provided. Please enter valid numeric values.\n")
    sys.exit(1)
except IndexError:
    sys.stderr.write("Error: Missing one or more input arguments. Please provide three integer values.\n")
    sys.exit(1)


with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
left_chain_input_file = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36', 'A_configuration', 'left-unit.pdb')
right_chain_input_file = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'charmm36', 'A_configuration', 'right-unit.pdb')
l_u = mda.Universe(left_chain_input_file)
r_u = mda.Universe(right_chain_input_file)




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




 #left_layer assembly  
left_layer_input_file = "left-chain_temp.2.pdb"
layer_l_u = mda.Universe(left_layer_input_file)
left_layer = []  
for i in range(1, a_iterations + 1):
    layer_l_u.atoms.positions += [a_trans, 0, 0]

    segid_increment_value = 1 
    for atom in layer_l_u.atoms:
        numeric_part = atom.segid.replace('', '')
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
    
    segid_increment_value = 1
    for atom in layer_r_u.atoms:
        numeric_part = atom.segid.replace('', '')
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


left_layer_files = []
right_layer_files = []

left_origin_sheet_u = mda.Universe("left-layer_temp.pdb")
right_origin_sheet_u = mda.Universe("right-layer_temp.pdb")

left_sheet_u = left_origin_sheet_u.copy()
left_sheet_u.atoms.positions += [2*a_trans,  0, 0]
left_sheet_u.atoms.write(f"left-sheet_temp.pdb")

right_sheet_u = right_origin_sheet_u.copy()
#right_sheet_u.atoms.positions += [ a_trans,  b_trans, 0]
right_sheet_u.atoms.write(f"right-sheet_temp.pdb")

# Iterate and create translated structures
for i in range(1, b_iterations + 1):
    # Load universes each iteration
    left_sheet_update_u = mda.Universe("left-sheet_temp.pdb")
    right_sheet_update_u = mda.Universe("right-sheet_temp.pdb")

    # Translate atoms
    left_sheet_update_u.atoms.positions += [2*i*a_trans, i * b_trans, 0]
    right_sheet_update_u.atoms.positions += [2*i*a_trans, i * b_trans, 0]

    # Write individual translated PDBs
    left_output = f"left_layer_{i}.pdb"
    right_output = f"right_layer_{i}.pdb"
    left_layer_files.append(left_output)
    right_layer_files.append(right_output)

    with mda.Writer(left_output, left_sheet_update_u.atoms.n_atoms) as W:
        W.write(left_sheet_update_u.atoms)

    with mda.Writer(right_output, right_sheet_update_u.atoms.n_atoms) as W:
        W.write(right_sheet_update_u.atoms)

# Assemble left-layer_temp.pdb
with open("left-sheet_all_temp.pdb", "w") as left_assembled:
    for left_file in left_layer_files:
        with open(left_file, "r") as f:
            for line in f:
                if line.startswith("ATOM"):
                    left_assembled.write(line)
        os.remove(left_file)  # Clean up individual pdb file

# Assemble right-layer_temp.pdb
with open("right-sheet_all_temp.pdb", "w") as right_assembled:
    for right_file in right_layer_files:
        with open(right_file, "r") as f:
            for line in f:
                if line.startswith("ATOM"):
                    right_assembled.write(line)
        os.remove(right_file)  # Clean up individual pdb file

# Merge all layers into one final PDB file
final_files = ["left-sheet_all_temp.pdb", "right-sheet_all_temp.pdb"]

with open("all_temp.pdb", "w") as final_pdb:
    for fname in final_files:
        with open(fname, "r") as infile:
            final_pdb.write(infile.read() + "\n")


final_structure=mda.Universe("all_temp.pdb")
with mda.Writer("charmm36-alpha-chitin-A-icm.pdb", n_atoms=final_structure.atoms.n_atoms, reindex=True) as W:
    W.write(final_structure.atoms)
    
# Cleanup any remaining temporary files
for temp_file in glob.glob("*temp*.pdb"):
    os.remove(temp_file)

print("Generated alpha-chitin-A with infinite-chain model.")