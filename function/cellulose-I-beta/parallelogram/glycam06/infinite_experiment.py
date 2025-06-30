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




with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
chain_1 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'glycam06', 'chain-1.pdb')
chain_2 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'glycam06', 'chain-2.pdb')

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


#
#
u_final = mda.Universe("layer_after_temp.pdb")
box_x = num_h * 6 + 20
box_y = num_w * 6 + 20
box_z = c_iterations * c_trans
center_of_mass = u_final.atoms.center_of_mass()
center_of_box = [box_x / 2, box_y / 2, box_z / 2]
u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
translation_vector = center_of_box - center_of_mass
u_final.atoms.translate(translation_vector)
with mda.Writer("glycam06-cellulose-Ibeta-para-icm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True) as W:
    W.write(u_final.atoms)

# Cleanup any remaining temporary files
for temp_file in glob.glob("*temp*.pdb"):
    os.remove(temp_file)

