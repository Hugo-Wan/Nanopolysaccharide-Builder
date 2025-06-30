import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis import Universe
import math
import warnings
import os
import glob
import json
from vmd import atomsel, molecule
import subprocess
warnings.filterwarnings("ignore", category=UserWarning)


 
structure_input_file = 'glycam06-cellulose-Ibeta-fcm.pdb'
u = mda.Universe(structure_input_file)


##-------------------Box dimension------------------
dims = u.dimensions
x=dims[0] + 50
y=dims[1] + 50
z=dims[2] + 20
coords = f"{x} {y} {z}"

#print(f"Box X length: {x:.3f} Å")
#print(f"Box Y length: {y:.3f} Å")
#print(f"Box Z length: {z:.3f} Å")
com = u.atoms.center_of_mass()
box_center = np.array([x/2, y/2, z/2])
translation = box_center - com
u.atoms.translate(translation)
u.dimensions = [x, y, z, dims[3], dims[4], dims[5]]
u.atoms.write("cellulose.temp.pdb")
##-------------------Box dimension------------------

##-------------------adding ions------------------
cation_name   = "Na+"
cation_number = 110

anion_name    = "Cl-"
anion_number  = 13
##-------------------adding ions------------------


structure_for_amber="cellulose.temp.pdb"
structure = molecule.load("pdb", structure_for_amber)
all_atoms = atomsel("all", molid=structure)
all_fragment = set(all_atoms.fragment)
all_resid=set(all_atoms.resid)
all_resid_last = sorted(all_resid)[-1]
num_fragment = len(all_fragment)
#print(num_fragment)
#print(all_resid_last)
end=num_fragment

def get_fragment_indices(mol_id, start_fragment):
    structure_index = []
    fragment_selection = atomsel(f"fragment {start_fragment}", molid=mol_id)
    indices = fragment_selection.index
    structure_index.extend(indices)
    return structure_index

def count_residues(pdb_filename):
    universe = mda.Universe(pdb_filename)
    num_residues = len(universe.residues)
    return num_residues


for i in range(0, num_fragment):
    selected_indices = get_fragment_indices(structure, int(i))
    ndx_i = selected_indices[0]
    ndx_j = selected_indices[-1]
    j=i+1
    #print(f"process chain {j}")
    selected_fragment = u.select_atoms(f"index {ndx_i} : {ndx_j}")
    selected_fragment.write(f"cellulose_temp_{j}.pdb") 



def generate_tleap_combine_command(num_glycans):
    glycans = [f"glycan_{i+1}" for i in range(num_glycans)]
    output_lines = ["source leaprc.GLYCAM_06j-1\n",
                    "source leaprc.water.tip3p\n"]
    for i, glycan in enumerate(glycans):
        pdb_file = f"cellulose_temp_{i+1}.pdb"
        load_and_bond = (
            f"{glycan} = loadpdb {pdb_file}\n"
        )
        output_lines.append(load_and_bond)
    
    combine_command = f"\n\ncombined_glycan = combine {{ {' '.join(glycans)} }}"
    output_lines.append(combine_command)
    

    output_lines.extend([
        "\nsolvatebox combined_glycan TIP3PBOX 20.0",
        "\ncheck combined_glycan", 
        f"\naddions combined_glycan {anion_name} {anion_number}",
        f"\naddions combined_glycan {cation_name} {cation_number}",
        "\nsaveAmberParm combined_glycan cell_solv.prmtop cell_solv.inpcrd",
        "\nsavepdb combined_glycan cell_solv.pdb",
        "quit"
    ])
    
    # Write the commands to a file
    with open("seq.in", "w") as file:
        file.write("\n".join(output_lines))

num_glycans = num_fragment  
generate_tleap_combine_command(num_glycans)

command = ['tleap', '-f', 'seq.in']
try:
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    print("tleap executed successfully:")
    print(result.stdout)  
except subprocess.CalledProcessError as e:
    print("An error occurred while running tleap:")
    print(e.stderr)  

for temp_file in glob.glob("*temp*.pdb"):
    os.remove(temp_file)
if os.path.exists("seq.in"):
    os.remove("seq.in")