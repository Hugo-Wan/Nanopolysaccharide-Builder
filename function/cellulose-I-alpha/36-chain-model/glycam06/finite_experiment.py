import MDAnalysis as mda
from MDAnalysis.coordinates import PDB
from MDAnalysis.core.universe import Merge
from MDAnalysis.core.topologyattrs import Elements
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

a_trans = 6.717
b_trans = 5.962
c_trans = 10.40
alpha_angle=118.08
beta_angle=114.80
gamma_angle=80.37
v=333.3374   ##angstrom^3

try:
    c_iterations = int(sys.argv[1])
    if c_iterations < 0:
        sys.stderr.write("Error: Please provide a positive integer for c repetition number.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length value provide. Please enter a valid numeric value.\n")
    sys.exit(1)



with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
chain_1 = os.path.join(main_folder_path, 'structure', 'cellulose_I_alpha', 'glycam06', 'chain-finite.pdb')

chain_1_u = mda.Universe(chain_1)


#chain_1 assembly    
#------------------------------------------------------------------------------
chain_1_strand = []
for chain_1_i in range(1, c_iterations + 1):
    chain_1_u.atoms.positions += [c_trans, 0, 0]
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
    'O1': chain_1_first_atom.position +  [1.291, -0.166, -0.593],
    'HO1': chain_1_first_atom.position + [1.969, -0.022,  0.060],
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
#
##------------------------------------------------------------------------------
###translate parameter

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
z_dir_vertical_move = -1 * math.sqrt((h)**2- z_dir_trans_move**2)
z_dir_vertical_to_screen = -1 *  math.sqrt(b_trans**2 - (h**2))


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

u_final = mda.Universe("unit_temp.pdb")
#u_final.atoms.positions = np.dot(u_final.atoms.positions - u_final.atoms.center_of_mass(), rot_matrix) + u_final.atoms.center_of_mass()
#center_of_mass = u_final.atoms.center_of_mass()
#
#box_z = 80
#box_y = 80
#box_x = c_iterations * c_trans + 30
#center_of_box = [box_x/2, box_y /2, box_z/2]
#
#translation_vector = center_of_box - center_of_mass
#u_final.atoms.translate(translation_vector)
#u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
with mda.Writer("glycam06-cellulose-Ialpha-36fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True) as W:
    W.write(u_final.atoms)
for temp_file in glob.glob("*temp*.pdb"):
    os.remove(temp_file)
