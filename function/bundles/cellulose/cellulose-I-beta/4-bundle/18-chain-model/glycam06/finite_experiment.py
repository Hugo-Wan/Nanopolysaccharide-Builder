import MDAnalysis as mda
from MDAnalysis.coordinates import PDB
from MDAnalysis.core.universe import Merge
from MDAnalysis.core.topologyattrs import Elements
from vmd import atomsel, molecule
from MDAnalysis.analysis import distances
import warnings
import math
import os
import glob
import json
import sys
import re
import numpy as np
warnings.filterwarnings("ignore", category=UserWarning)
from psfgen import PsfGen
from numpy.linalg import svd
from numpy.linalg import norm


a_trans = 7.784
b_trans = 8.201
c_trans = 10.38
gamma_angle=96.5



try:
    length = float(sys.argv[1]) 
    c_iterations = int(length // c_trans)
    if c_iterations < 1:
        sys.stderr.write("Error: Please provide an larger length value.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length value provide. Please enter a valid numeric value.\n")
    sys.exit(1)
    
try:
    d_up_center = float(sys.argv[2])

    if d_up_center < 0 :
        raise ValueError("Inter-crystal distance (100) must be equivalent or greater than 0.")
    
    d_center_ur = float(sys.argv[3])

    if d_center_ur < 0 :
        raise ValueError("Inter-crystal distance (1-10) must be equivalent or greater than 0.")

    d_center_ul = float(sys.argv[4])
    if  d_center_ul < 0:
        raise ValueError("Inter-crystal distance (110) must be equivalent or greater than 0.")

except ValueError as e:
    sys.stderr.write(f"Error: {str(e)}\n")
    sys.exit(1)



with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
chain_1 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'glycam06', 'chain-1-finite.pdb')
chain_2 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'glycam06', 'chain-2-finite.pdb')

chain_1_u = mda.Universe(chain_1)
chain_2_u = mda.Universe(chain_2)



#chain_1 assembly    
#------------------------------------------------------------------------------
chain_1_strand = []
for chain_1_i in range(1, c_iterations + 1):
    chain_1_u.atoms.positions += [0, 0, c_trans]
    resid_1 = (c_iterations - chain_1_i + 1) * 2 
    resid_2 = (c_iterations - chain_1_i + 1) * 2 + 1
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
chain_1_h4_atoms = chain_1_finite_u.select_atoms("name H4")
if not chain_1_o4_atoms:
    raise ValueError("No atoms named 'O4' found in the structure.")
chain_1_last_o4 = chain_1_o4_atoms[-1]
chain_1_last_h4 = chain_1_h4_atoms[-1]


chain_1_new_positions_ROH = {
    'O1': chain_1_first_atom.position + [0.538, -0.449, 1.247],
    'HO1': chain_1_first_atom.position + [0.129, -0.365, 1.932],
}

chain_1_new_positions_0GB= {
    'O4':chain_1_last_o4.position,
    'H4O': chain_1_last_o4.position + [0.377, -0.192, -0.892]
}


chain_1_new_atoms_ROH = {}
for name, pos in chain_1_new_positions_ROH.items():
    if chain_1_first_atom.name in ['O1', 'HO1']:
        chain_1_base_atom = chain_1_first_atom
    else:
        chain_1_base_atom = chain_1_finite_u.atoms[0]  
    chain_1_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
    chain_1_new_uni.add_TopologyAttr('name', [name])
    chain_1_new_uni.add_TopologyAttr('type', [chain_1_base_atom.type])
    chain_1_new_uni.add_TopologyAttr('resname', ['ROH'])
    chain_1_new_uni.add_TopologyAttr('resid', ['1'])
    chain_1_new_uni.add_TopologyAttr('segid', ['0'])
    chain_1_new_uni.atoms.positions = [pos]
    chain_1_new_atoms_ROH[name] = chain_1_new_uni

chain_1_new_atoms_0GB = {}
for name, pos in chain_1_new_positions_0GB.items():
    chain_1_base_atom = chain_1_last_o4  # Using the last O4 for properties
    chain_1_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
    chain_1_new_uni.add_TopologyAttr('name', [name])
    chain_1_new_uni.add_TopologyAttr('type', [chain_1_base_atom.type])
    chain_1_new_uni.add_TopologyAttr('resname', ['4GB'])  # Resname set to '0GB' as specified
    chain_1_new_uni.add_TopologyAttr('resid', [chain_1_base_atom.resid])
    chain_1_new_uni.add_TopologyAttr('segid', ['0'])
    chain_1_new_uni.atoms.positions = [pos]
    chain_1_new_atoms_0GB[name] = chain_1_new_uni


chain_1_combined = Merge(chain_1_finite_u.atoms[:chain_1_last_h4.index + 1], chain_1_new_atoms_0GB['O4'].atoms, chain_1_finite_u.atoms[chain_1_last_h4.index + 1:])
chain_1_combined = Merge(chain_1_combined.atoms[:chain_1_last_h4.index + 2], chain_1_new_atoms_0GB['H4O'].atoms, chain_1_combined.atoms[chain_1_last_h4.index + 2:])

chain_1_combined = Merge(chain_1_new_atoms_ROH['HO1'].atoms, chain_1_combined.atoms)
chain_1_combined = Merge(chain_1_combined.atoms[:1], chain_1_new_atoms_ROH['O1'].atoms, chain_1_combined.atoms[1:])
chain_1_all_o4_atoms = chain_1_combined.select_atoms("name O4")
chain_1_last_o4_to_remove = chain_1_all_o4_atoms[-1]
chain_1_atoms_to_keep = chain_1_combined.select_atoms("not (index %d)" % chain_1_last_o4_to_remove.index)
chain_1_atoms_to_keep.write("chain_1_temp.2.pdb")


chain_1_final=mda.Universe("chain_1_temp.2.pdb")
for atom in chain_1_final.residues[-3:].atoms:  
    atom.residue.resname = '0GB'



for atom in chain_1_final.atoms:
    atom.segment.segid = '0'  
    elements = [atom.name[0] for atom in chain_1_final.atoms]
    chain_1_final.add_TopologyAttr(Elements(elements))
chain_1_final.atoms.write("chain_1_temp.2.pdb")
#------------------------------------------------------------------------------



#chain_2 assembly    
#------------------------------------------------------------------------------
chain_2_strand = []
for chain_2_i in range(1, c_iterations + 1):
    chain_2_u.atoms.positions += [0, 0, c_trans]
    resid_1 = (c_iterations - chain_2_i + 1) * 2 
    resid_2 = (c_iterations - chain_2_i + 1) * 2 + 1
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
chain_2_h4_atoms = chain_2_finite_u.select_atoms("name H4")
if not chain_2_o4_atoms:
    raise ValueError("No atoms named 'O4' found in the structure.")
chain_2_last_o4 = chain_2_o4_atoms[-1]
chain_2_last_h4 = chain_2_h4_atoms[-1]

chain_2_new_positions_ROH = {
    'O1': chain_2_first_atom.position + [0.538, -0.449, 1.247],
    'HO1': chain_2_first_atom.position + [0.129, -0.365, 1.932],
}

chain_2_new_positions_0GB= {
    'O4':chain_2_last_o4.position,
    'H4O': chain_2_last_o4.position + [0.377, -0.192, -0.892]
}


chain_2_new_atoms_ROH = {}
for name, pos in chain_2_new_positions_ROH.items():
    if chain_2_first_atom.name in ['O1', 'HO1']:
        chain_2_base_atom = chain_2_first_atom
    else:
        chain_2_base_atom = chain_2_finite_u.atoms[0]  
    chain_2_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
    chain_2_new_uni.add_TopologyAttr('name', [name])
    chain_2_new_uni.add_TopologyAttr('type', [chain_2_base_atom.type])
    chain_2_new_uni.add_TopologyAttr('resname', ['ROH'])
    chain_2_new_uni.add_TopologyAttr('resid', ['1'])
    chain_2_new_uni.add_TopologyAttr('segid', ['0'])
    chain_2_new_uni.atoms.positions = [pos]
    chain_2_new_atoms_ROH[name] = chain_2_new_uni

chain_2_new_atoms_0GB = {}
for name, pos in chain_2_new_positions_0GB.items():
    chain_2_base_atom = chain_2_last_o4  # Using the last O4 for properties
    chain_2_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
    chain_2_new_uni.add_TopologyAttr('name', [name])
    chain_2_new_uni.add_TopologyAttr('type', [chain_2_base_atom.type])
    chain_2_new_uni.add_TopologyAttr('resname', ['4GB'])  # Resname set to '0GB' as specified
    chain_2_new_uni.add_TopologyAttr('resid', [chain_2_base_atom.resid])
    chain_2_new_uni.add_TopologyAttr('segid', ['0'])
    chain_2_new_uni.atoms.positions = [pos]
    chain_2_new_atoms_0GB[name] = chain_2_new_uni


chain_2_combined = Merge(chain_2_finite_u.atoms[:chain_2_last_h4.index + 1], chain_2_new_atoms_0GB['O4'].atoms, chain_2_finite_u.atoms[chain_2_last_h4.index + 1:])
chain_2_combined = Merge(chain_2_combined.atoms[:chain_2_last_h4.index + 2], chain_2_new_atoms_0GB['H4O'].atoms, chain_2_combined.atoms[chain_2_last_h4.index + 2:])

chain_2_combined = Merge(chain_2_new_atoms_ROH['HO1'].atoms, chain_2_combined.atoms)
chain_2_combined = Merge(chain_2_combined.atoms[:1], chain_2_new_atoms_ROH['O1'].atoms, chain_2_combined.atoms[1:])

chain_2_all_o4_atoms = chain_2_combined.select_atoms("name O4")
chain_2_last_o4_to_remove = chain_2_all_o4_atoms[-1]
chain_2_atoms_to_keep = chain_2_combined.select_atoms("not (index %d)" % chain_2_last_o4_to_remove.index)
chain_2_atoms_to_keep.write("chain_2_temp.2.pdb")


chain_2_final=mda.Universe("chain_2_temp.2.pdb")
for atom in chain_2_final.residues[-3:].atoms:  
    atom.residue.resname = '0GB'



for atom in chain_2_final.atoms:
    atom.segment.segid = '0'  
    elements = [atom.name[0] for atom in chain_2_final.atoms]
    chain_2_final.add_TopologyAttr(Elements(elements))
chain_2_final.atoms.write("chain_2_temp.2.pdb")

#------------------------------------------------------------------------------
#
##------------------------------------------------------------------------------
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


#chain_1_1st_layer assembly  
chain_1_1st_layer_input_file = "chain_1_temp.2.pdb"
chain_1_1st_layer_u = mda.Universe(chain_1_1st_layer_input_file)
chain_1_1st_layer = []  
for i in range(1, 5):
    chain_1_1st_layer_u.atoms.positions += [ b_par_vertical_move, b_par_transverse_move, 0]

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


#chain_1_2nd_layer assembly  
chain_1_2nd_layer_input_file = "chain_1_temp.2.pdb"
chain_1_2nd_layer_u = mda.Universe(chain_1_1st_layer_input_file)
chain_1_2nd_layer = [] 

for j in range(1, 3):
    k=j+1
    chain_1_2nd_layer_u = mda.Universe(chain_1_2nd_layer_input_file)
    translation_vector = [a_par_vertical_move + b_par_vertical_move * k,
                          a_par_transverse_move + b_par_transverse_move * k, 0]
    chain_1_2nd_layer_u.atoms.positions += translation_vector
    for atom in chain_1_2nd_layer_u.atoms:
        new_numeric_part = 4 + j  
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


#chain_1_3rd_layer assembly  
chain_1_3rd_layer_input_file = "chain_1_temp.2.pdb"
chain_1_3rd_layer_u = mda.Universe(chain_1_1st_layer_input_file)
chain_1_3rd_layer = [] 

for l in range(1, 4):  
    m = l - 1  
    chain_1_3rd_layer_u = mda.Universe(chain_1_3rd_layer_input_file)
    translation_vector = [-a_par_vertical_move + b_par_vertical_move * l,
                          -a_par_transverse_move + b_par_transverse_move * l, 0]

    chain_1_3rd_layer_u.atoms.positions += translation_vector

    for atom in chain_1_3rd_layer_u.atoms:
        new_numeric_part = 6 + l  
        atom.segment.segid = str(new_numeric_part)  
    chain_1_3rd_layer_output = f"chain_1_3rd_layer_{l}.pdb"
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




#chain_2_1st_layer assembly  
chain_2_1st_layer_input_file = "chain_2_temp.2.pdb"
chain_2_1st_layer_u = mda.Universe(chain_2_1st_layer_input_file)
chain_2_1st_layer = []  

for i in range(1, 5):  
    chain_2_1st_layer_u = mda.Universe(chain_2_1st_layer_input_file)
    translation_vector = [b_par_vertical_move * i,
                          b_par_transverse_move * i, 0]

    chain_2_1st_layer_u.atoms.positions += translation_vector

    for atom in chain_2_1st_layer_u.atoms:
        new_numeric_part = 9 + i  
        atom.segment.segid = str(new_numeric_part)  
    chain_2_1st_layer_output = f"chain_2_1st_layer_{i}.pdb"
    chain_2_1st_layer.append(chain_2_1st_layer_output)
    with mda.Writer(chain_2_1st_layer_output, n_atoms=chain_2_1st_layer_u.atoms.n_atoms) as W:
        W.write(chain_2_1st_layer_u.atoms)
with open("chain_2_1st_layer_temp.pdb", "w") as layer_file:
    for chain_2_1st_layer_output in chain_2_1st_layer:
        with open(chain_2_1st_layer_output, "r") as chain_2_1st_layer_pdb_file:
            for line in chain_2_1st_layer_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
            os.remove(chain_2_1st_layer_output)




#chain_2_2nd_layer assembly  
chain_2_2nd_layer_input_file = "chain_2_temp.2.pdb"
chain_2_2nd_layer_u = mda.Universe(chain_2_1st_layer_input_file)
chain_2_2nd_layer = [] 

for j in range(1, 4):
    k=j+1
    chain_2_2nd_layer_u = mda.Universe(chain_2_2nd_layer_input_file)
    translation_vector = [a_par_vertical_move + b_par_vertical_move * k,
                          a_par_transverse_move + b_par_transverse_move * k, 0]
    chain_2_2nd_layer_u.atoms.positions += translation_vector
    for atom in chain_2_2nd_layer_u.atoms:
        new_numeric_part = 13 + j  
        atom.segment.segid = str(new_numeric_part)  
    chain_2_2nd_layer_output = f"chain_2_2nd_layer_{j}.pdb"
    chain_2_2nd_layer.append(chain_2_2nd_layer_output)
    with mda.Writer(chain_2_2nd_layer_output, n_atoms=chain_2_2nd_layer_u.atoms.n_atoms) as W:
        W.write(chain_2_2nd_layer_u.atoms)
with open("chain_2_2nd_layer_temp.pdb", "w") as layer_file:
    for chain_2_2nd_layer_output in chain_2_2nd_layer:
        with open(chain_2_2nd_layer_output, "r") as chain_2_2nd_layer_pdb_file:
            for line in chain_2_2nd_layer_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
            os.remove(chain_2_2nd_layer_output)


#chain_2_3rd_layer assembly  
chain_2_3rd_layer_input_file = "chain_2_temp.2.pdb"
chain_2_3rd_layer_u = mda.Universe(chain_2_1st_layer_input_file)
chain_2_3rd_layer = [] 

for l in range(1, 3):  
    m = l + 1  
    chain_2_3rd_layer_u = mda.Universe(chain_2_3rd_layer_input_file)
    translation_vector = [-a_par_vertical_move + b_par_vertical_move * m,
                          -a_par_transverse_move + b_par_transverse_move * m, 0]

    chain_2_3rd_layer_u.atoms.positions += translation_vector

    for atom in chain_2_3rd_layer_u.atoms:
        new_numeric_part = 16 + l  
        atom.segment.segid = str(new_numeric_part)  
    chain_2_3rd_layer_output = f"chain_2_3rd_layer_{l}.pdb"
    chain_2_3rd_layer.append(chain_2_3rd_layer_output)
    with mda.Writer(chain_2_3rd_layer_output, n_atoms=chain_2_3rd_layer_u.atoms.n_atoms) as W:
        W.write(chain_2_3rd_layer_u.atoms)
with open("chain_2_3rd_layer_temp.pdb", "w") as layer_file:
    for chain_2_3rd_layer_output in chain_2_3rd_layer:
        with open(chain_2_3rd_layer_output, "r") as chain_2_3rd_layer_pdb_file:
            for line in chain_2_3rd_layer_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
            os.remove(chain_2_3rd_layer_output)


###assemble
def rotate_z(coordinates, angle_degrees):
    """Rotate coordinates around the z-axis by the specified angle in degrees."""
    angle_radians = np.radians(angle_degrees)
    rotation_matrix = np.array([
        [np.cos(angle_radians), -np.sin(angle_radians), 0],
        [np.sin(angle_radians), np.cos(angle_radians), 0],
        [0, 0, 1]
    ])
    return np.dot(coordinates, rotation_matrix.T)

def assemble_pdbs(output_file, input_files):
    # Initialize a universe from the first file to start with
    universe = mda.Universe(input_files[0])
    
    # Load additional files into the same universe as new segments
    for pdb_file in input_files[1:]:
        additional = mda.Universe(pdb_file)
        universe = mda.Merge(universe.atoms, additional.atoms)
    
    # Write combined file
    universe.atoms.write(output_file)
    return universe

# Define file names
input_files = [
    "chain_1_1st_layer_temp.pdb", "chain_1_2nd_layer_temp.pdb", "chain_1_3rd_layer_temp.pdb",
    "chain_2_1st_layer_temp.pdb", "chain_2_2nd_layer_temp.pdb", "chain_2_3rd_layer_temp.pdb"
]
combined_pdb = "combined_temp.pdb"

# Assemble all PDBs into one
combined_universe = assemble_pdbs(combined_pdb, input_files)

##up-center
#move 3a
step_1_move_vertical= (3* a_trans) + d_up_center  #+ 1 * b_trans * sinA
step_1_move_transverse=  0



###center-upper-right
step_3_move_vertical_coef= 2* a_trans   + 3 * b_trans * sinA
step_3_move_transverse_coef = -3 * b_trans * cosA

step_3_move_vertical = step_3_move_vertical_coef
step_3_move_transverse = step_3_move_transverse_coef
if d_center_ur > 0:

    step_3_ratio = abs(step_3_move_vertical_coef / step_3_move_transverse_coef)
    current_total_distance = math.sqrt(step_3_move_vertical_coef**2 + step_3_move_transverse_coef**2)
    desired_total_distance = current_total_distance + d_center_ur
    step_3_theta = math.atan(step_3_ratio)
    
    # Calculate the changes needed for the movements
    step3_v_add = (desired_total_distance * math.sin(step_3_theta)) - abs(step_3_move_vertical_coef)
    step3_t_add = (desired_total_distance * math.cos(step_3_theta)) - abs(step_3_move_transverse_coef)


    step_3_move_vertical = abs(step_3_move_vertical_coef) + step3_v_add
    step_3_move_transverse = abs(step_3_move_transverse_coef) + step3_t_add
    step_3_move_vertical *= math.copysign(1, step_3_move_vertical_coef)
    step_3_move_transverse *= math.copysign(1, step_3_move_transverse_coef)


###center-upper-left
step_5_move_vertical_coef= 1* a_trans  - 3 * b_trans * sinA
step_5_move_transverse_coef= 3 * b_trans * cosA

step_5_move_vertical = step_5_move_vertical_coef
step_5_move_transverse = step_5_move_transverse_coef

if d_center_ul > 0:
    step_5_ratio = abs(step_5_move_vertical_coef / step_5_move_transverse_coef)
    current_total_distance = math.sqrt(step_5_move_vertical_coef**2 + step_5_move_transverse_coef**2)
    desired_total_distance = current_total_distance + d_center_ul
    step_5_theta = math.atan(step_5_ratio)
    
    # Calculate the changes needed for the movements
    step3_v_add = (desired_total_distance * math.sin(step_5_theta)) - abs(step_5_move_vertical_coef)
    step3_t_add = (desired_total_distance * math.cos(step_5_theta)) - abs(step_5_move_transverse_coef)


    step_5_move_vertical = abs(step_5_move_vertical_coef) + step3_v_add
    step_5_move_transverse = abs(step_5_move_transverse_coef) + step3_t_add
    step_5_move_vertical *= math.copysign(1, step_5_move_vertical_coef)
    step_5_move_transverse *= math.copysign(1, step_5_move_transverse_coef)





##structure-1 100-100 (up)
u_original = mda.Universe("combined_temp.pdb")
translation_vector_1 = [step_1_move_vertical, step_1_move_transverse , 0]
u_duplicate_1 = u_original.copy()
u_duplicate_1.atoms.translate(translation_vector_1)


##structure-3 110-110 (right-upper)
u_original = mda.Universe("combined_temp.pdb")
translation_vector_3 = [step_3_move_vertical, step_3_move_transverse , 0]
u_duplicate_3 = u_original.copy()
u_duplicate_3.atoms.translate(translation_vector_3)


##structure-5 110-110 (up-left)
u_original = mda.Universe("combined_temp.pdb")
translation_vector_5 = [step_5_move_vertical, step_5_move_transverse , 0]
u_duplicate_5 = u_original.copy()
u_duplicate_5.atoms.translate(translation_vector_5)




duplicates = [u_duplicate_1,u_duplicate_3,  u_duplicate_5]
increments = [18,36,54]  # Define increment values

for u_duplicate, increment in zip(duplicates, increments):
    for seg in u_duplicate.segments:
        try:
            new_id = str(int(seg.segid) + increment)  
        except ValueError:
            new_id = f"{seg.segid}_{increment}"  
        for atom in seg.atoms:
            atom.segment.segid = new_id
u_combined = mda.Merge(u_original.atoms, u_duplicate_1.atoms,  u_duplicate_3.atoms, u_duplicate_5.atoms)




# Save the combined structure to a new PDB file
with mda.Writer("glycam06-c1b-4-bundles-fcm.pdb", n_atoms=u_combined.atoms.n_atoms, reindex=True) as W:
    W.write(u_combined.atoms)

# Cleanup temporary files
for temp_file in glob.glob("*temp*.pdb"):
    os.remove(temp_file)
