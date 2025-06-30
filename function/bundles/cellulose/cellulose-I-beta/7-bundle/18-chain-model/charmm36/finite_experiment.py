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
    d_down_center = float(sys.argv[3])
    if d_up_center < 0 or d_down_center < 0:
        raise ValueError("Inter-crystal distance (100) must be equivalent or greater than 0.")
    
    d_center_ur = float(sys.argv[4])
    d_center_dl = float(sys.argv[7])
    if d_center_ur < 0 or d_center_dl < 0:
        raise ValueError("Inter-crystal distance (1-10) must be equivalent or greater than 0.")

    d_center_dr = float(sys.argv[5])
    d_center_ul = float(sys.argv[6])
    if d_center_dr < 0 or d_center_ul < 0:
        raise ValueError("Inter-crystal distance (110) must be equivalent or greater than 0.")

except ValueError as e:
    sys.stderr.write(f"Error: {str(e)}\n")
    sys.exit(1)



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
    'HO1': chain_1_first_atom.position + [-0.129, -0.365, 1.932],
    'HO4': chain_1_last_o4.position + [-0.617, -0.192, -0.892]
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
    'HO1': chain_2_first_atom.position + [-0.129, -0.365, 1.932],
    'HO4': chain_2_last_o4.position + [-0.617, -0.192, -0.892]
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

###move_1=24.31345 + up_down_add   ###up and  down nanostructures ###24.31345 make sure the minimum distance between H (closet) is equivalent to 2.4 angstrom (vdw radii for H)
##up -down
step_2_move_vertical= 3* -a_trans - d_up_center  #- 1 * b_trans * sinA
step_2_move_transverse=   0 #-2* b_trans * cosA


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

###center-down-right
step_4_move_vertical_coef= -1* a_trans   + 3 * b_trans * sinA
step_4_move_transverse_coef= -3 * b_trans * cosA

step_4_move_vertical = step_4_move_vertical_coef
step_4_move_transverse = step_4_move_transverse_coef

if d_center_dr > 0:
    step_4_ratio = abs(step_4_move_vertical_coef / step_4_move_transverse_coef)
    current_total_distance = math.sqrt(step_4_move_vertical_coef**2 + step_4_move_transverse_coef**2)
    desired_total_distance = current_total_distance + d_center_dr
    step_4_theta = math.atan(step_4_ratio)
    
    # Calculate the changes needed for the movements
    step3_v_add = (desired_total_distance * math.sin(step_4_theta)) - abs(step_4_move_vertical_coef)
    step3_t_add = (desired_total_distance * math.cos(step_4_theta)) - abs(step_4_move_transverse_coef)


    step_4_move_vertical = abs(step_4_move_vertical_coef) + step3_v_add
    step_4_move_transverse = abs(step_4_move_transverse_coef) + step3_t_add
    step_4_move_vertical *= math.copysign(1, step_4_move_vertical_coef)
    step_4_move_transverse *= math.copysign(1, step_4_move_transverse_coef)


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



###center-bottom-left
step_6_move_vertical_coef = -2* a_trans  - 3 * b_trans * sinA
step_6_move_transverse_coef = 3 * b_trans * cosA

step_6_move_vertical = step_6_move_vertical_coef
step_6_move_transverse = step_6_move_transverse_coef

if d_center_dl > 0:
    step_6_ratio = abs(step_6_move_vertical_coef / step_6_move_transverse_coef)
    current_total_distance = math.sqrt(step_6_move_vertical_coef**2 + step_6_move_transverse_coef**2)
    desired_total_distance = current_total_distance + d_center_dl
    step_6_theta = math.atan(step_6_ratio)
    
    # Calculate the changes needed for the movements
    step3_v_add = (desired_total_distance * math.sin(step_6_theta)) - abs(step_6_move_vertical_coef)
    step3_t_add = (desired_total_distance * math.cos(step_6_theta)) - abs(step_6_move_transverse_coef)


    step_6_move_vertical = abs(step_6_move_vertical_coef) + step3_v_add
    step_6_move_transverse = abs(step_6_move_transverse_coef) + step3_t_add
    step_6_move_vertical *= math.copysign(1, step_6_move_vertical_coef)
    step_6_move_transverse *= math.copysign(1, step_6_move_transverse_coef)



##structure-1 100-100 (up)
u_original = mda.Universe("combined_temp.pdb")
translation_vector_1 = [step_1_move_vertical, step_1_move_transverse , 0]
u_duplicate_1 = u_original.copy()
u_duplicate_1.atoms.translate(translation_vector_1)

##structure-2 100-100 (down)
u_original = mda.Universe("combined_temp.pdb")
translation_vector_2 = [step_2_move_vertical, step_2_move_transverse , 0]
u_duplicate_2 = u_original.copy()
u_duplicate_2.atoms.translate(translation_vector_2)

##structure-3 110-110 (right-upper)
u_original = mda.Universe("combined_temp.pdb")
translation_vector_3 = [step_3_move_vertical, step_3_move_transverse , 0]
u_duplicate_3 = u_original.copy()
u_duplicate_3.atoms.translate(translation_vector_3)

##structure-4 110-110 (down-right)
u_original = mda.Universe("combined_temp.pdb")
translation_vector_4 = [step_4_move_vertical, step_4_move_transverse , 0]
u_duplicate_4 = u_original.copy()
u_duplicate_4.atoms.translate(translation_vector_4)



##structure-5 110-110 (up-left)
u_original = mda.Universe("combined_temp.pdb")
translation_vector_5 = [step_5_move_vertical, step_5_move_transverse , 0]
u_duplicate_5 = u_original.copy()
u_duplicate_5.atoms.translate(translation_vector_5)


##structure-6 110-110 (up-left)
u_original = mda.Universe("combined_temp.pdb")
translation_vector_6 = [step_6_move_vertical, step_6_move_transverse , 0]
u_duplicate_6 = u_original.copy()
u_duplicate_6.atoms.translate(translation_vector_6)



duplicates = [u_duplicate_1, u_duplicate_2, u_duplicate_3,  u_duplicate_4,  u_duplicate_5, u_duplicate_6]
increments = [18,36,54,72,90, 108]  # Define increment values

for u_duplicate, increment in zip(duplicates, increments):
    for seg in u_duplicate.segments:
        try:
            new_id = str(int(seg.segid) + increment)  
        except ValueError:
            new_id = f"{seg.segid}_{increment}"  
        for atom in seg.atoms:
            atom.segment.segid = new_id
u_combined = mda.Merge(u_original.atoms, u_duplicate_1.atoms, u_duplicate_2.atoms, u_duplicate_3.atoms,u_duplicate_4.atoms, u_duplicate_5.atoms, u_duplicate_6.atoms)




# Save the combined structure to a new PDB file
with mda.Writer("charmm36-c1b-7-bundles-fcm.pdb", n_atoms=u_combined.atoms.n_atoms, reindex=True) as W:
    W.write(u_combined.atoms)

# Cleanup temporary files
for temp_file in glob.glob("*temp*.pdb"):
    os.remove(temp_file)


#structure_input_file="charmm36-cellulose-Ibeta-18fcm.pdb"
#structure = molecule.load("pdb", structure_input_file)
#
#def get_H_indices(mol_id, start_fragment,end_fragment):
#    structure_index = []
#    fragment_selection = atomsel(f"fragment {start_fragment} {end_fragment} and name H1 H2 H3 H4 H5", molid=mol_id)
#    indices = fragment_selection.index
#    structure_index.extend(indices)
#    return structure_index
#
#def get_structure_indices(mol_id, start_fragment,end_fragment):
#    structure_index = []
#    fragment_selection = atomsel(f"fragment {start_fragment} {end_fragment}", molid=mol_id)
#    indices = fragment_selection.index
#    structure_index.extend(indices)
#    return structure_index
#
#
#selected_indices_1 = get_H_indices(structure, int(16), int(17))
#selected_indices_2 = get_H_indices(structure, int(22), int(23))
#selected_indices_3 = get_H_indices(structure, int(4), int(5))
#selected_indices_4 = get_H_indices(structure, int(52), int(53))
#
#selected_cel1 = get_structure_indices(structure, int(0), int(17))
#selected_cel2 = get_structure_indices(structure, int(18), int(35))
#selected_cel3 = get_structure_indices(structure, int(36), int(53))
#
#
#selected_layer1 = get_structure_indices(structure, int(16), int(17))
#selected_layer2 = get_structure_indices(structure, int(22), int(23))
#
#
#
#cel=mda.Universe("charmm36-cellulose-Ibeta-18fcm.pdb")
#
#
#def calculate_minimum_distance_and_atoms(cel, indices1, indices2):
#    """Calculate the minimum distance between two sets of atom indices in the XY plane and identify the atoms involved."""
#    group1 = cel.atoms[indices1]
#    group2 = cel.atoms[indices2]
#    # Set Z component to zero to ignore it
#    pos1 = np.copy(group1.positions)
#    pos2 = np.copy(group2.positions)
#    pos1[:, 2] = 0
#    pos2[:, 2] = 0
#    dist_matrix = distances.distance_array(pos1, pos2)
#    min_dist = dist_matrix.min()
#    min_idx = np.unravel_index(np.argmin(dist_matrix), dist_matrix.shape)
#    atom1 = group1[min_idx[0]]
#    atom2 = group2[min_idx[1]]
#    return min_dist, atom1, atom2
#
#def calculate_com_distance(cel, indices1, indices2):
#    """Calculate the minimum distance between two sets of atom indices and identify the atoms involved."""
#    group1 = cel.atoms[indices1]
#    group2 = cel.atoms[indices2]
# 
#    com1 = group1.center_of_geometry()
#    com2 = group2.center_of_geometry()
#    com_distance = distances.distance_array(com1.reshape(1, 3), com2.reshape(1, 3))[0][0]
#    return com_distance
#
#
#def calculate_com(cel, indices1, indices2):
#    """Calculate the minimum distance between two sets of atom indices and identify the atoms involved."""
#    group1 = cel.atoms[indices1]
#    group2 = cel.atoms[indices2]
#    com1 = group1.center_of_geometry()
#    com2 = group2.center_of_geometry()
#    return com1, com2
#
#
#
##up and center  100-100
#min_dist, atom1, atom2 = calculate_minimum_distance_and_atoms(cel, selected_indices_1, selected_indices_2)
#print(f"Up and center minimum distance: {min_dist} Å between {atom1.name} (id: {atom1.id}) and {atom2.name} (id: {atom2.id})")
#com_dist = calculate_com_distance(cel, selected_cel1, selected_cel2)
#print(f"Up and center geo-center distance: {com_dist} Å ")
#
#
##down and center  100-100
#min_dist2, atom3, atom4 = calculate_minimum_distance_and_atoms(cel, selected_indices_3, selected_indices_4)
#print(f"down and center minimum distance: {min_dist2} Å between {atom3.name} (id: {atom3.id}) and {atom4.name} (id: {atom4.id})")
#com_dist = calculate_com_distance(cel, selected_cel1, selected_cel3)
#print(f"down and center geo-center distance: {com_dist} Å ")

#
#com = calculate_com(cel, selected_layer1, selected_layer2)
#print(f"{com} ")





#u_final = mda.Universe("unit_temp.pdb")
#box_x = 50 
#box_y = 50 
#box_z = c_iterations * c_trans + 20
#center_of_mass = u_final.atoms.center_of_mass()
#center_of_box = [box_x / 2, box_y / 2, box_z / 2]
#u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
#translation_vector = center_of_box - center_of_mass
#u_final.atoms.translate(translation_vector)