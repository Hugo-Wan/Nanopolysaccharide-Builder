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
    a_iterations = int(sys.argv[8]) 
    b_iterations = int(sys.argv[9]) 
    c_iterations = int(sys.argv[10]) 
    if c_iterations < 0 or b_iterations < 0 or a_iterations < 0:
        sys.stderr.write("Error: Please provide a positive integer for all repetition numbers.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length value provide. Please enter a valid numeric value.\n")
    sys.exit(1)


with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
neutron_pdb_1 = os.path.join(main_folder_path, 'structure', 'cellulose_I_alpha', 'glycam06', 'chain.pdb')
user_defined_pdb_1 = os.path.join(main_folder_path, 'structure', 'cellulose_I_alpha', 'glycam06', 'chain_ud.pdb')

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

chain_1_combined=mda.Universe("chain_1_temp.pdb")
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


try:
    h = b_trans * cos_alpha_angle  
    z_dir_trans_move = v / (a_trans * c_trans * sin_beta_angle_raw)

    value_to_be_squared_1 = (h)**2 - z_dir_trans_move**2
    if value_to_be_squared_1 < 0:
        raise ValueError("Error with a, c, beta angle parameter or single cell volume.")

    z_dir_vertical_move = -1 * math.sqrt(value_to_be_squared_1)
    value_to_be_squared_2 = b_trans**2 - (h)**2
    if value_to_be_squared_2 < 0:
        raise ValueError("Error with a, b, or alpha angle parameter.")

    z_dir_vertical_to_screen = -1 * math.sqrt(value_to_be_squared_2)
except ValueError as e:
    for temp_file in glob.glob("*temp*.pdb"):
        os.remove(temp_file)
    print(str(e), file=sys.stderr)
    sys.exit(1)



#chain_1_1st_layer assembly  
chain_1_1st_layer_input_file = "chain_1_temp.2.pdb"
chain_1_1st_layer_u = mda.Universe(chain_1_1st_layer_input_file)
chain_1_1st_layer = []  
for i in range(1, a_iterations+1):
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



chain_1_sheet_input_file = "chain_1_1st_layer_temp.pdb"
chain_1_sheet = [] 

chain_1_base_universe = mda.Universe(chain_1_sheet_input_file)

for j in range(1, b_iterations+1):
    chain_1_sheet_u = chain_1_base_universe.copy()

    translation_vector = [z_dir_vertical_to_screen*(j) ,
                          z_dir_vertical_move*(j) ,
                          z_dir_trans_move*(j)
    ]

    chain_1_sheet_u.atoms.positions += translation_vector

    segid_increment_value = (j-1)* a_iterations 
    for segid, group in chain_1_sheet_u.atoms.groupby('segids').items():
        new_segid = str(int(segid) + segid_increment_value)
        for atom in group:
            atom.segment.segid = new_segid

    chain_1_sheet_output = f"chain_1_sheet_{j}.pdb"
    chain_1_sheet.append(chain_1_sheet_output)
    with mda.Writer(chain_1_sheet_output, n_atoms=chain_1_sheet_u.atoms.n_atoms) as W:
        W.write(chain_1_sheet_u.atoms)
with open("chain_1_sheet_temp.pdb", "w") as layer_file:
    for chain_1_sheet_output in chain_1_sheet:
        with open(chain_1_sheet_output, "r") as chain_1_sheet_pdb_file:
            for line in chain_1_sheet_pdb_file:
                if line.startswith("ATOM"):
                    layer_file.write(line)
                elif line.startswith("TER"): 
                    layer_file.write("TER\n")
        os.remove(chain_1_sheet_output)  




####assemble
def assemble_pdbs(output_file, input_files):
    with open(output_file, 'w') as outfile:
        for filename in input_files:
            with open(filename, 'r') as infile:
                for line in infile:
                    if line.startswith("ATOM"):
                        outfile.write(line)
            outfile.write('TER\n')  

chain_temp_pdb_files = ["chain_1_sheet_temp.pdb"]

unit_temp_pdb = "unit_temp.pdb"
assemble_pdbs(unit_temp_pdb, chain_temp_pdb_files)


u_final = mda.Universe("unit_temp.pdb")
center_of_mass = u_final.atoms.center_of_mass()

box_z = b_iterations * b_trans
box_y = a_iterations * a_trans
box_x = c_iterations * c_trans



gamma_degrees_rad = (180-gamma_angle) * math.pi / 180
alpha_degrees_rad = (180-alpha_angle) * math.pi / 180
beta_degrees_rad = (180-beta_angle) * math.pi / 180
cos_alpha_degrees_rad=math.cos(alpha_degrees_rad)
cos_gamma_degrees_rad=math.sin(gamma_degrees_rad)
cos_beta_degrees_rad=math.tan(beta_degrees_rad)

#+box_z*cos_gamma_degrees_rad  -+
center_of_box = [(box_x-box_z*cos_beta_degrees_rad)/2, (box_y)/2, box_z/2]

translation_vector = center_of_box - center_of_mass
u_final.atoms.translate(translation_vector)
u_final.dimensions = [ box_x, box_y, box_z, gamma_angle, alpha_angle, beta_angle]
with mda.Writer("glycam06-cellulose-Ialpha-icm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True) as W:
    W.write(u_final.atoms)
for temp_file in glob.glob("*temp*.pdb"):
    os.remove(temp_file)

