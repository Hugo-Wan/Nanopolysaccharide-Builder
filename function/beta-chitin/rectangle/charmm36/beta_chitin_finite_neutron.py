import MDAnalysis as mda
from MDAnalysis.coordinates import PDB
from MDAnalysis.core.universe import Merge
from MDAnalysis.core.topologyattrs import Elements
import warnings
import os
import glob
import math
import json
import sys
import re
warnings.filterwarnings("ignore", category=UserWarning)


a_trans     = 4.819 
b_trans     = 9.238976
c_trans     = 10.384001
gamma_angle = 97.156586

try:
    a_iterations = int(sys.argv[1])
    b_iterations = int(sys.argv[2])
    c_iterations = int(sys.argv[3])
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
unit_chain_input_file = os.path.join(main_folder_path, 'structure', 'beta_chitin', 'charmm36', 'unit-neutron.pdb')

unit_u = mda.Universe(unit_chain_input_file)


 #chain assembly    
 #------------------------------------------------------------------------------
unit_chain = []
for unit_i in range(1, c_iterations + 1):
    unit_u.atoms.positions -= [0, 0, c_trans]
    resid_1 = unit_i  * 2 - 1
    resid_2 = unit_i  * 2
    for j, atom in enumerate(unit_u.atoms):
        if j < 27:
            atom.residue.resid = resid_1
        else:
            atom.residue.resid = resid_2
    unit_chain_output = f"unit-{unit_i}.pdb"
    unit_chain.append(unit_chain_output)
    with mda.Writer(unit_chain_output, n_atoms=unit_u.atoms.n_atoms) as W:
        W.write(unit_u.atoms)
with open("unit-chain_temp.pdb", "w") as chain:
    for unit_chain_output in unit_chain:
         with open(unit_chain_output, "r") as unit_chain_pdb_file:
            for unit_chain_line in unit_chain_pdb_file:
                if unit_chain_line.startswith("ATOM"):
                    chain.write(unit_chain_line)
            os.remove(unit_chain_output)


 ##terminal side            
unit_finite_build_chain = "unit-chain_temp.pdb"
unit_finite_u = mda.Universe(unit_finite_build_chain)
unit_first_atom = unit_finite_u.atoms[0]
o4_atoms = unit_finite_u.select_atoms("name O4")
if not o4_atoms:
    raise ValueError("No atoms named 'O4' found in the structure.")
last_o4 = o4_atoms[-1]

new_positions = {
    'O1': unit_first_atom.position + [0.715, 0.167, 1.227],
    'HO1': unit_first_atom.position + [0.163, -0.114, 1.960],
    'HO4': last_o4.position + [0.915, -0.066, 0.281]
}

unit_new_atoms = {}
for name, pos in new_positions.items():
    unit_base_atom = unit_first_atom if name in ['O1', 'HO1'] else last_o4
    unit_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
    unit_new_uni.add_TopologyAttr('name', [name])
    unit_new_uni.add_TopologyAttr('type', [unit_base_atom.type])
    unit_new_uni.add_TopologyAttr('resname', [unit_base_atom.resname])
    unit_new_uni.add_TopologyAttr('resid', [unit_base_atom.resid])
    unit_new_uni.add_TopologyAttr('chainIDs', ['X'])
    unit_new_uni.add_TopologyAttr('segid', ['0']) 
    unit_new_uni.atoms.positions = [pos]
    unit_new_atoms[name] = unit_new_uni

unit_combined = Merge(unit_finite_u.atoms[:last_o4.index + 1], unit_new_atoms['HO4'].atoms, unit_finite_u.atoms[last_o4.index + 1:])
unit_combined = Merge(unit_combined.atoms[:2], unit_new_atoms['O1'].atoms, unit_combined.atoms[2:])
unit_combined = Merge(unit_combined.atoms[:3], unit_new_atoms['HO1'].atoms, unit_combined.atoms[3:])



for atom in unit_combined.atoms:
    atom.segment.segid = '0'  
    elements = [atom.name[0] for atom in unit_combined.atoms]
    unit_combined.add_TopologyAttr(Elements(elements))
unit_combined.atoms.write("unit-chain_temp.2.pdb")
#------------------------------------------------------------------------------


 #layer assembly 
unit_layer_input_file = "unit-chain_temp.2.pdb"
layer_unit_u = mda.Universe(unit_layer_input_file)
unit_layer = []  
for i in range(1, a_iterations + 1):
    layer_unit_u.atoms.positions += [a_trans, 0, 0]

    segid_increment_value = 1 
    for atom in layer_unit_u.atoms:
        numeric_part = atom.segid.replace('','')
        new_numeric_part = int(numeric_part) + segid_increment_value
    new_segid = f"{new_numeric_part}"
    atom.segment.segid = new_segid

    unit_layer_output = f"unit_layer_{i}.pdb"
    unit_layer.append(unit_layer_output)
    with mda.Writer(unit_layer_output, n_atoms=layer_unit_u.atoms.n_atoms) as W:
        W.write(layer_unit_u.atoms)
with open("unit-layer_temp.pdb", "w") as layer_file:
    for unit_layer_output in unit_layer:
        with open(unit_layer_output, "r") as unit_layer_pdb_file:
            for unit_layer_line in unit_layer_pdb_file:
                if unit_layer_line.startswith("ATOM"):
                    layer_file.write(unit_layer_line)
            os.remove(unit_layer_output) 
 


unit_output_file = "unit.pdb" 
with open(unit_output_file, "w") as unit_file:
    with open("unit-layer_temp.pdb", "r") as assembled_file:
        for line in assembled_file:
            if line.startswith("ATOM"):
                unit_file.write(line)


unit_input_file="unit.pdb"
crystal_structure = mda.Universe(unit_input_file)
crystal_structure_files = [] 


angle_rad = (gamma_angle - 90) * math.pi / 180
cosA = math.cos(angle_rad)
sinA = math.sin(angle_rad)

for crystal_i in range(1, b_iterations + 1):
    crystal_structure.atoms.positions += [-b_trans * sinA, b_trans * cosA , 0]
    crystal_structure_output = f"crystal_{crystal_i}.pdb"
    crystal_structure_files.append(crystal_structure_output)
    with mda.Writer(crystal_structure_output, n_atoms=crystal_structure.atoms.n_atoms) as writer:
        writer.write(crystal_structure.atoms)


for file_index, file_name in enumerate(crystal_structure_files):
    temp_structure = mda.Universe(file_name)
    segid_increment_value = file_index * a_iterations
    for segment in temp_structure.segments:
        matches = re.findall(r'\d+', segment.segid)
        if matches:
            numeric_part = int(matches[0])
            new_numeric_part = numeric_part + segid_increment_value
            new_segid = f"{new_numeric_part}"
            for atom in segment.atoms:
                atom.segment.segid = new_segid
    with mda.Writer(file_name, n_atoms=temp_structure.atoms.n_atoms) as writer:
        writer.write(temp_structure.atoms)
with open("unit_temp.pdb", "w") as combined_file:
    for file_name in crystal_structure_files:
        with open(file_name, "r") as file:
            for line in file:
                if line.startswith("ATOM"):
                    combined_file.write(line)
        os.remove(file_name) 

os.remove(unit_input_file) 

u_final = mda.Universe("unit_temp.pdb")
box_x = a_iterations * a_trans + 20
box_y = b_iterations * b_trans + 20
box_z = c_iterations * c_trans + 20
center_of_mass = u_final.atoms.center_of_mass()
center_of_box = [box_x / 2, box_y / 2, box_z / 2]
u_final.dimensions = [box_x, box_y, box_z, 90, 90, 90]
translation_vector = center_of_box - center_of_mass
u_final.atoms.translate(translation_vector)
with mda.Writer("charmm36-beta-chitin-fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True) as W:
    W.write(u_final.atoms)

 #Cleanup any remaining temporary files
for temp_file in glob.glob("*temp*.pdb"):
    os.remove(temp_file)

