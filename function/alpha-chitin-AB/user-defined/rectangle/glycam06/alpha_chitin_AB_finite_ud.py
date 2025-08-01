import MDAnalysis as mda
from MDAnalysis.coordinates import PDB
from MDAnalysis.core.universe import Merge
from MDAnalysis.core.topologyattrs import Elements
import warnings
import os
import glob
import json
import sys
import re
warnings.filterwarnings("ignore", category=UserWarning)
from psfgen import PsfGen


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
neutron_pdb_1 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'glycam06', 'AB_configuration', 'left-unit-finite.pdb')
user_defined_pdb_1 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'glycam06', 'AB_configuration', 'left-unit-finite_ud.pdb')

if os.path.exists(user_defined_pdb_1):
    unit_chain_input_file_1 = user_defined_pdb_1
else:
    unit_chain_input_file_1 = neutron_pdb_1

neutron_pdb_2 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'glycam06','AB_configuration', 'right-unit-finite.pdb')
user_defined_pdb_2 = os.path.join(main_folder_path, 'structure', 'alpha_chitin', 'glycam06', 'AB_configuration', 'right-unit-finite_ud.pdb')

if os.path.exists(user_defined_pdb_2):
    unit_chain_input_file_2 = user_defined_pdb_2
else:
    unit_chain_input_file_2 = neutron_pdb_2

l_u = mda.Universe(unit_chain_input_file_1)
r_u = mda.Universe(unit_chain_input_file_2)


 #left_chain assembly    
 #------------------------------------------------------------------------------
left_chain = []
for l_i in range(1, c_iterations + 1):
    l_u.atoms.positions += [0, 0, c_trans]
    resid_1 = (c_iterations - l_i + 1) * 2 
    resid_2 = (c_iterations - l_i + 1) * 2 + 1
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
        with open(l_chain_output, "r") as left_chain_pdb_file:
            for left_chain_line in left_chain_pdb_file:
                if left_chain_line.startswith("ATOM"):
                    l_chain.write(left_chain_line)
            os.remove(l_chain_output) 

 ##terminal side            
left_finite_build_chain = "left-chain_temp.pdb"
left_finite_u = mda.Universe(left_finite_build_chain)
left_first_atom = left_finite_u.atoms[0]
left_o4_atoms = left_finite_u.select_atoms("name O4")
left_H4_atoms = left_finite_u.select_atoms("name H4")

if not left_o4_atoms:
    raise ValueError("No atoms named 'O4' found in the structure.")
left_last_o4 = left_o4_atoms[-1]
left_last_h4 = left_H4_atoms[-1]


left_new_positions_ROH = {
    'O1': left_first_atom.position + [0.538, -0.449, 1.347],
    'HO1': left_first_atom.position + [0.129, -0.365, 1.932],
}

left_new_positions_0YB= {
    'O4':left_last_o4.position,
    'H4O': left_last_o4.position + [0.377, -0.192, -0.892]
}



left_new_atoms_ROH = {}
for name, pos in left_new_positions_ROH.items():
    if left_first_atom.name in ['O1', 'HO1']:
        left_base_atom = left_first_atom
    else:
        left_base_atom = left_finite_u.atoms[0]  
    left_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
    left_new_uni.add_TopologyAttr('name', [name])
    left_new_uni.add_TopologyAttr('type', [left_base_atom.type])
    left_new_uni.add_TopologyAttr('resname', ['ROH'])
    left_new_uni.add_TopologyAttr('resid', ['1'])
    left_new_uni.add_TopologyAttr('segid', ['0'])
    left_new_uni.atoms.positions = [pos]
    left_new_atoms_ROH[name] = left_new_uni



left_new_atoms_0YB = {}
for name, pos in left_new_positions_0YB.items():
    left_base_atom = left_last_o4  # Using the last O4 for properties
    left_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
    left_new_uni.add_TopologyAttr('name', [name])
    left_new_uni.add_TopologyAttr('type', [left_base_atom.type])
    left_new_uni.add_TopologyAttr('resname', ['4YB'])  # Resname set to '0GB' as specified
    left_new_uni.add_TopologyAttr('resid', [left_base_atom.resid])
    left_new_uni.add_TopologyAttr('segid', ['0'])
    left_new_uni.atoms.positions = [pos]
    left_new_atoms_0YB[name] = left_new_uni

left_combined = Merge(left_finite_u.atoms[:left_last_h4.index + 1], left_new_atoms_0YB['O4'].atoms,  left_finite_u.atoms[left_last_h4.index + 1:])
left_combined = Merge(left_combined.atoms[:left_last_h4.index + 2], left_new_atoms_0YB['H4O'].atoms, left_combined.atoms[left_last_h4.index + 2:])

left_combined = Merge(left_new_atoms_ROH['HO1'].atoms, left_combined.atoms)
left_combined = Merge(left_combined.atoms[:1], left_new_atoms_ROH['O1'].atoms, left_combined.atoms[1:])
left_all_o4_atoms = left_combined.select_atoms("name O4")
left_last_o4_to_remove = left_all_o4_atoms[-1]
left_atoms_to_keep = left_combined.select_atoms("not (index %d)" % left_last_o4_to_remove.index)
left_atoms_to_keep.write("left_chain_temp.1.pdb")


infile  = "left_chain_temp.1.pdb"
outfile = infile  # overwrite in place

with open(infile, 'r') as f:
    lines = f.readlines()

h1m_inds = [i for i, L in enumerate(lines)
             if L.startswith(("ATOM  ", "HETATM")) and L[12:16].strip() == "H1M"]
h2m_inds = [i for i, L in enumerate(lines)
             if L.startswith(("ATOM  ", "HETATM")) and L[12:16].strip() == "H2M"]

if h1m_inds and h2m_inds:
    last_h1 = h1m_inds[-1]
    last_h2 = h2m_inds[-1]
    h1m_line = lines.pop(last_h1)
    if last_h1 < last_h2:
        last_h2 -= 1
    lines.insert(last_h2 + 1, h1m_line)
    with open(outfile, 'w') as f:
        f.writelines(lines)



left_final=mda.Universe("left_chain_temp.1.pdb")
for atom in left_final.residues[-3:].atoms:  
    atom.residue.resname = '0YB'

for atom in left_final.atoms:
    atom.segment.segid = '0'  
    elements = [atom.name[0] for atom in left_final.atoms]
    left_final.add_TopologyAttr(Elements(elements))
left_final.atoms.write("left-chain_temp.2.pdb")

#------------------------------------------------------------------------------




##right_chain assembly      
##------------------------------------------------------------------------------ 
right_chain = []
right_side_start_segid=a_iterations
for r_i in range(1, c_iterations + 1):
    r_u.atoms.positions += [0, 0, c_trans]
    resid_1 = r_i  * 2  
    resid_2 = r_i  * 2 + 1
    for j, atom in enumerate(r_u.atoms):
        if j < 27:
            atom.residue.resid = resid_1
        else:
            atom.residue.resid = resid_2
    right_chain_output = f"r-{r_i}.pdb"
    right_chain.append(right_chain_output)
    with mda.Writer(right_chain_output, n_atoms=r_u.atoms.n_atoms) as W:
        W.write(r_u.atoms)
with open("right-chain_temp.pdb", "w") as chain:
    for r_chain_output in right_chain:
         with open(r_chain_output, "r") as right_chain_pdb_file:
            for right_chain_line in right_chain_pdb_file:
                if right_chain_line.startswith("ATOM"):
                    chain.write(right_chain_line)
            os.remove(r_chain_output)


 ##terminal side            
right_finite_build_chain = "right-chain_temp.pdb"
right_finite_u = mda.Universe(right_finite_build_chain)
right_first_atom = right_finite_u.atoms[0]
right_o4_atoms = right_finite_u.select_atoms("name O4")
right_H4_atoms = right_finite_u.select_atoms("name H4")

if not right_o4_atoms:
    raise ValueError("No atoms named 'O4' found in the structure.")
right_last_o4 = right_o4_atoms[-1]
right_last_h4 = right_H4_atoms[-1]



right_new_positions_ROH = {
    'O1': right_first_atom.position - [-0.538, -0.449, 1.247],
    'HO1': right_first_atom.position - [0.129, -0.365, 1.932],
}

right_new_positions_0YB= {
    'O4':right_last_o4.position,
    'H4O': right_last_o4.position - [0.377, -0.192, -0.892]
}




right_new_atoms_ROH = {}
for name, pos in right_new_positions_ROH.items():
    if right_first_atom.name in ['O1', 'HO1']:
        right_base_atom = right_first_atom
    else:
        right_base_atom = right_finite_u.atoms[0]  
    right_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
    right_new_uni.add_TopologyAttr('name', [name])
    right_new_uni.add_TopologyAttr('type', [right_base_atom.type])
    right_new_uni.add_TopologyAttr('resname', ['ROH'])
    right_new_uni.add_TopologyAttr('resid', ['1'])
    right_new_uni.add_TopologyAttr('segid', ['0'])
    right_new_uni.atoms.positions = [pos]
    right_new_atoms_ROH[name] = right_new_uni



right_new_atoms_0YB = {}
for name, pos in right_new_positions_0YB.items():
    right_base_atom = right_last_o4  # Using the last O4 for properties
    right_new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
    right_new_uni.add_TopologyAttr('name', [name])
    right_new_uni.add_TopologyAttr('type', [right_base_atom.type])
    right_new_uni.add_TopologyAttr('resname', ['4YB'])  # Resname set to '0GB' as specified
    right_new_uni.add_TopologyAttr('resid', [right_base_atom.resid])
    right_new_uni.add_TopologyAttr('segid', ['0'])
    right_new_uni.atoms.positions = [pos]
    right_new_atoms_0YB[name] = right_new_uni

right_combined = Merge(right_finite_u.atoms[:right_last_h4.index + 1], right_new_atoms_0YB['O4'].atoms,  right_finite_u.atoms[right_last_h4.index + 1:])
right_combined = Merge(right_combined.atoms[:right_last_h4.index + 2], right_new_atoms_0YB['H4O'].atoms, right_combined.atoms[right_last_h4.index + 2:])

right_combined = Merge(right_new_atoms_ROH['HO1'].atoms, right_combined.atoms)
right_combined = Merge(right_combined.atoms[:1], right_new_atoms_ROH['O1'].atoms, right_combined.atoms[1:])
right_all_o4_atoms = right_combined.select_atoms("name O4")
right_last_o4_to_remove = right_all_o4_atoms[-1]
right_atoms_to_keep = right_combined.select_atoms("not (index %d)" % right_last_o4_to_remove.index)
right_atoms_to_keep.write("right_chain_temp.1.pdb")


infile  = "right_chain_temp.1.pdb"
outfile = infile  # overwrite in place

with open(infile, 'r') as f:
    lines = f.readlines()

h1m_inds = [i for i, L in enumerate(lines)
             if L.startswith(("ATOM  ", "HETATM")) and L[12:16].strip() == "H1M"]
h2m_inds = [i for i, L in enumerate(lines)
             if L.startswith(("ATOM  ", "HETATM")) and L[12:16].strip() == "H2M"]

if h1m_inds and h2m_inds:
    last_h1 = h1m_inds[-1]
    last_h2 = h2m_inds[-1]
    h1m_line = lines.pop(last_h1)
    if last_h1 < last_h2:
        last_h2 -= 1
    lines.insert(last_h2 + 1, h1m_line)
    with open(outfile, 'w') as f:
        f.writelines(lines)



right_final=mda.Universe("right_chain_temp.1.pdb")
for atom in right_final.residues[-3:].atoms:  
    atom.residue.resname = '0YB'
right_side_start_segid=a_iterations
for atom in right_final.atoms:
    atom.segment.segid = f'{right_side_start_segid}'  
    elements = [atom.name[0] for atom in right_final.atoms]
    right_final.add_TopologyAttr(Elements(elements))
right_final.atoms.write("right-chain_temp.2.pdb")
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


 #left_layer assembly  
left_layer_input_file = "left-chain_temp.2.pdb"
layer_l_u = mda.Universe(left_layer_input_file)
left_layer = []  
for i in range(1, a_iterations + 1):
    layer_l_u.atoms.positions += [a_trans, 0, 0]

    segid_increment_value = 1 
    for atom in layer_l_u.atoms:
        numeric_part = atom.segid.replace('','')
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
    #j=i+a_iterations
    segid_increment_value = 1
    for atom in layer_r_u.atoms:
        numeric_part = atom.segid.replace('','')
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


 #unit_structure assembly  
left_layer_input_file = "left-layer_temp.pdb"
right_layer_input_file = "right-layer_temp.pdb"
unit_output_file = "unit.pdb" 
with open(unit_output_file, "w") as unit_file:
    with open("left-layer_temp.pdb", "r") as left_file:
        for line in left_file:
            if line.startswith("ATOM"):
                unit_file.write(line)
    with open("right-layer_temp.pdb", "r") as right_file:
        for line in right_file:
            if line.startswith("ATOM"):
                unit_file.write(line)


unit_input_file="unit.pdb"
crystal_structure = mda.Universe(unit_input_file)
crystal_structure_files = [] 

for crystal_i in range(1, b_iterations + 1):
    crystal_structure.atoms.positions += [0, -b_trans, 0]
    crystal_structure_output = f"crystal_{crystal_i}.pdb"
    crystal_structure_files.append(crystal_structure_output)
    with mda.Writer(crystal_structure_output, n_atoms=crystal_structure.atoms.n_atoms) as writer:
        writer.write(crystal_structure.atoms)


for file_index, file_name in enumerate(crystal_structure_files):
    temp_structure = mda.Universe(file_name)
    segid_increment_value = file_index * 2 * a_iterations
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
with mda.Writer("glycam06-alpha-chitin-AB-fcm.pdb", n_atoms=u_final.atoms.n_atoms, reindex=True) as W:
    W.write(u_final.atoms)

# Cleanup any remaining temporary files
for temp_file in glob.glob("*temp*.pdb"):
    os.remove(temp_file)

