import MDAnalysis as mda
from MDAnalysis.coordinates import PDB
from MDAnalysis.core.universe import Merge
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.coordinates.PDB import PDBWriter
import sys
import numpy as np
import math
import warnings
import os
import glob
import json


warnings.filterwarnings("ignore", category=UserWarning)


if len(sys.argv) != 5:
    print("Usage: python infinite-chain.py DDA pH unit_distance DP")
    sys.exit(1)

# Assign variables from command line arguments
DDA_target = float(sys.argv[3])                 # Deacetylation degree
pH = float(sys.argv[4])                         # pH condition
c_trans  = float(sys.argv[2])                   # Unit distance
c_iterations  = int(sys.argv[1])                # Degree of Polymerization



#DDA_target=0.55                #deacetylation degree ref:10.1021/acs.jchemed.7b00902J
#pH=3                           #pH condition              
#c_trans = 10.33
#c_iterations = 4



pKa=6.3    #pka of chitosan  ref:doi.org/10.1016/j.foodhyd.2022.108383
m_GlcNAc = 204.20  # GlcNAc unit mass
m_Glc = 179.17     # Glc unit mass


with open('config.json', 'r') as f:
    config = json.load(f)

# Get the main folder path from the configuration
main_folder_path = config['main_folder_path']
# Construct the path to the unit file
native_unit_input_file = os.path.join(main_folder_path, 'structure', 'chitosan', 'unit.pdb')
u = mda.Universe(native_unit_input_file)


chain_filenames = [] 

for i in range(1, c_iterations + 1):
    u.atoms.positions += [0, 0, c_trans]
    resid_1 = (c_iterations - i + 1) * 2 - 1
    resid_2 = (c_iterations - i + 1) * 2
    for j, atom in enumerate(u.atoms):
        if j < 27:
            atom.residue.resid = resid_1
        else:
            atom.residue.resid = resid_2
    chain_output = f"{i}.pdb"
    chain_filenames.append(chain_output)  # Append to the list of filenames
    with mda.Writer(chain_output, n_atoms=u.atoms.n_atoms) as W:
        W.write(u.atoms)

# Now use the list of filenames to write to the combined PDB file
with open("chitin-temp.pdb", "w") as chain_file:
    for chain_output in reversed(chain_filenames):
        with open(chain_output, "r") as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    chain_file.write(line)
        os.remove(chain_output)



chain_u = mda.Universe("chitin-temp.pdb")
residue_ids = [residue.resid for residue in chain_u.residues]
total_residues=len(residue_ids)
#print(total_residues)


def calculate_dda(x):
    numerator = (total_residues - x) * m_Glc
    denominator = numerator + x * m_GlcNAc
    return numerator / denominator if denominator != 0 else 0
closest_x = None
closest_dda = float('inf')  
for x in range(total_residues + 1): 
    dda = calculate_dda(x)
    if abs(dda - DDA_target) < abs(closest_dda - DDA_target):
        closest_x = x
        closest_dda = dda   
#print(f"aminiation_number:{total_residues-closest_x}")
#print(f"Actual_DDA={closest_dda:.4f}")


num_residues_to_modify = total_residues-closest_x
atoms_to_keep = chain_u.atoms.copy()  # Start with all atoms
if num_residues_to_modify > len(chain_u.residues):
    print("Number of modifications exceeds available residues. Operation cancelled.")
    sys.exit()

modified_residues = set()
selected_residue_details = []  
nh3_atoms_to_keep = chain_u.atoms[:] 
for _ in range(num_residues_to_modify):
    while True:
        nh3_random_residue = np.random.choice(chain_u.residues)
        if nh3_random_residue.resid not in modified_residues:
            modified_residues.add(nh3_random_residue.resid)
            selected_residue_details.append(nh3_random_residue.resid)
            #print(f"Randomly selected residue: {random_residue.resname} with ID {random_residue.resid}")

            if nh3_random_residue.resid == 1:
                nh3_atoms_to_remove = nh3_random_residue.atoms[8:15]  # Specific handling for residue ID 1
            else:
                nh3_atoms_to_remove = nh3_random_residue.atoms[8:15]  # General handling for other residues

            # Remove selected atoms from the atom group to keep
            nh3_atoms_to_keep = nh3_atoms_to_keep - nh3_atoms_to_remove
            break  # Exit the while loop once a unique residue is processed

# Create a new Universe with the remaining atoms
new_u = mda.Merge(nh3_atoms_to_keep)
new_u.atoms.write("modified-chain-temp.pdb")




###NH3_step
modified_chain = "modified-chain-temp.pdb"
modified_chain_u = mda.Universe(modified_chain)
all_nh3_atoms = list(modified_chain_u.atoms)

for sel_resid in selected_residue_details:
    sel_residue = modified_chain_u.residues[sel_resid - 1]  
    base_atom_index =  7
    base_atom = sel_residue.atoms[base_atom_index]

    # Conditional new positions based on whether the resid is odd or even
    if sel_resid % 2 == 1:  # Check if the residue ID is odd
        new_positions = {
            'HN1': base_atom.position + np.array([ 0.924, -0.082,  0.374]),
            'HN2': base_atom.position + np.array([-0.667,  0.256,  0.700]),
            'HN3': base_atom.position + np.array([-0.095, -0.607, -0.789])
        }
    else:  # Residue ID is even
        new_positions = {
            'HN1': base_atom.position + np.array([-0.529,  0.492, -0.691]),
            'HN2': base_atom.position + np.array([ 0.890,  0.441,  0.115]),
            'HN3': base_atom.position + np.array([-0.494, -0.013,  0.869])
        }

    # Starting position for the first new atom
    insert_pos = all_nh3_atoms.index(base_atom) + 1

    # Create and insert each new atom right after the previous one
    for name, pos in new_positions.items():
        new_uni = mda.Universe.empty(n_atoms=1, trajectory=True)
        new_uni.add_TopologyAttr('name', [name])
        new_uni.add_TopologyAttr('type', [base_atom.type])
        new_uni.add_TopologyAttr('resname', [base_atom.resname])
        new_uni.add_TopologyAttr('resid', [sel_residue.resid])
        new_uni.add_TopologyAttr('segid', ['CARB'])
        new_uni.add_TopologyAttr('chainIDs', ['N'])
        new_uni.atoms.positions = [pos]
        new_atom = new_uni.atoms[0]

        # Insert the new atom at the calculated position
        all_nh3_atoms.insert(insert_pos, new_atom)
        insert_pos += 1  # Update position for the next atom

chitosan_nh3 = mda.Merge(*[mda.AtomGroup([atom]) for atom in all_nh3_atoms])
chitosan_nh3.atoms.write("nh3-temp.pdb")





#pH calculation |  NH2_step
amination_unit = total_residues - closest_x  # where closest_x is from your previous script
ratio = 10 ** (pH - pKa)
y = (amination_unit * ratio) / (1 + ratio)
y_rounded = round(y)

if  y_rounded > 0:
    ratio = y_rounded / amination_unit
    if ratio > 0:
        actual_pH = math.log10(ratio) + 6.3
        actual_pH_rounded = round(actual_pH, 2)
        #print("Actual pH:", actual_pH_rounded)
else:
    actual_pH_rounded = pH
#print(f"NH2 number is {y_rounded}")

nh2_modified_chain = mda.Universe("nh3-temp.pdb")
np.random.seed(42)  # For reproducibility
nh2_random_residues = np.random.choice(selected_residue_details, size=y_rounded, replace=False)

nh2_atom_groups = []
nh2_modified_resids = set()  

# Iterate over each residue
for nh2_residue in nh2_modified_chain.residues:
    if nh2_residue.resid in nh2_random_residues and nh2_residue.resid not in nh2_modified_resids:
        nh2_modified_resids.add(nh2_residue.resid)  # Mark this residue ID as modified
        if nh2_residue.resid == 1:
            # Check if there are enough atoms to perform operations safely
            if len(nh2_residue.atoms) > 10:
                nh2_atoms_to_keep = nh2_residue.atoms[:10] + nh2_residue.atoms[11:]
                nh2_atoms_to_keep[7].name = 'N2'
                nh2_atoms_to_keep[8].name = 'HN21'
                nh2_atoms_to_keep[9].name = 'HN22'
        else:
            if len(nh2_residue.atoms) > 10:
                nh2_atoms_to_keep = nh2_residue.atoms[:10] + nh2_residue.atoms[11:]
                nh2_atoms_to_keep[7].name = 'N2'
                nh2_atoms_to_keep[8].name = 'HN21'
                nh2_atoms_to_keep[9].name = 'HN22'
        nh2_atom_groups.append(nh2_atoms_to_keep)
    else:
        if nh2_residue.resid not in nh2_modified_resids:
           nh2_atom_groups.append(nh2_residue.atoms)


def adjust_box_and_write(universe, output_file, c_iterations, c_trans):
    # Calculate the minimum and maximum dimensions of the atoms' positions
    min_dim = universe.atoms.positions.min(axis=0)
    max_dim = universe.atoms.positions.max(axis=0)
    center_of_molecule = (max_dim + min_dim) / 2
    
    # Buffers for x and y dimensions
    buffer_xy = 10.0  # 1 nm buffer
    # Calculation for z dimension based on input parameters
    box_z = c_iterations * c_trans
    new_box_dimensions = max_dim - min_dim + np.array([2 * buffer_xy, 2 * buffer_xy, 0])
    new_box_center = new_box_dimensions / 2
    new_box_center[2] = box_z / 2  
    shift_vector = new_box_center - center_of_molecule
    universe.atoms.positions += shift_vector
    universe.dimensions = np.array([new_box_dimensions[0], new_box_dimensions[1], box_z, 90, 90, 90])

    with PDBWriter(output_file, universe.atoms.n_atoms) as pdb:
        pdb.write(universe.atoms)

if nh2_atom_groups: 
    new_universe = mda.Merge(*nh2_atom_groups)  
    new_universe.atoms.write("nh2-temp.pdb")
    chitosan_chain = mda.Universe("nh2-temp.pdb")  
    adjust_box_and_write(chitosan_chain, "chitosan-ifcm.pdb", c_iterations, c_trans)
else:
    chitosan_chain = mda.Universe("nh3-temp.pdb")  
    adjust_box_and_write(chitosan_chain, "chitosan-ifcm.pdb",c_iterations, c_trans)

# Remove temporary files
for temp_file in glob.glob("*temp*.pdb"):
    os.remove(temp_file)

print(f"DDA: {closest_dda:.2f}, pH: {actual_pH_rounded}, Units: {y_rounded}")