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
    gamma_angle = float(sys.argv[4])
    if a_trans <= 0 or b_trans <= 0 or c_trans <= 0 or gamma_angle<=0:  
        sys.stderr.write("Error: Please provide positive value for all crystallographic parameters.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid length or angle value provided. Please enter valid numeric values.\n")
    sys.exit(1)
except IndexError:
    sys.stderr.write("Error: Missing one or more input arguments. Please provide four required crystallographic values.\n")
    sys.exit(1)





try:
    c_iterations = int(sys.argv[5])
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
    height = float(sys.argv[6])
    num_h = int(round(height /7.8))
    if num_h < 1:
        sys.stderr.write("Error: Please provide an larger height value.\n")
        sys.exit(1)
except ValueError:
    sys.stderr.write("Error: Invalid height value provide. Please enter a valid numeric value.\n")
    sys.exit(1)    


try:
    width = float(sys.argv[7])
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
neutron_pdb_1 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'glycam06', 'chain-1-finite.pdb')
user_defined_pdb_1 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'glycam06', 'chain-1-finite_ud.pdb')

if os.path.exists(user_defined_pdb_1):
    unit_chain_input_file_1 = user_defined_pdb_1
else:
    unit_chain_input_file_1 = neutron_pdb_1

neutron_pdb_2 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'glycam06', 'chain-2-finite.pdb')
user_defined_pdb_2 = os.path.join(main_folder_path, 'structure', 'cellulose_I_beta', 'glycam06', 'chain-2-finite_ud.pdb')

if os.path.exists(user_defined_pdb_2):
    unit_chain_input_file_2 = user_defined_pdb_2
else:
    unit_chain_input_file_2 = neutron_pdb_2

chain_1_u = mda.Universe(unit_chain_input_file_1)
chain_2_u = mda.Universe(unit_chain_input_file_2)



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
#------------------------------------------------------------------------------
###translate parameter
tilt_angle_rad= (gamma_angle - 90) * math.pi/180
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
box_z = c_iterations * c_trans + 20
center_of_mass = u_final.atoms.center_of_mass()
center_of_box = [box_x / 2, box_y / 2, box_z / 2]
u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
translation_vector = center_of_box - center_of_mass
u_final.atoms.translate(translation_vector)
with mda.Writer("glycam06-cellulose-Ibeta-para-fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True) as W:
    W.write(u_final.atoms)

# Cleanup any remaining temporary files
for temp_file in glob.glob("*temp*.pdb"):
    os.remove(temp_file)

