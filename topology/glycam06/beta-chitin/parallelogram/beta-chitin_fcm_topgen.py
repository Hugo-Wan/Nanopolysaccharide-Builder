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





with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
structure_input_file = os.path.join(main_folder_path, 'glycam06-beta-chitin-para-fcm.pdb')
u = mda.Universe(structure_input_file)
structure = molecule.load("pdb", structure_input_file)
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
    selected_fragment.write(f"chitin_temp_{j}.pdb") 


def generate_tleap_combine_command(num_glycans, block_size=50):
    glycans = [f"glycan_{i+1}" for i in range(num_glycans)]
    output_lines = ["source leaprc.GLYCAM_06j-1", "source leaprc.water.tip3p", ""]
    combined_glycan_vars = []

    for glycan in glycans:
        pdb_file = f"{glycan.replace('glycan_', 'chitin_temp_')}.pdb"
        output_lines.append(f"{glycan} = loadpdb {pdb_file}")

    output_lines.append("") 

    for block_start in range(0, num_glycans, block_size):
        block_end = min(block_start + block_size, num_glycans)
        block_glycans = glycans[block_start:block_end]
        m = (block_start // block_size) + 1
        combined_name = f"combined_glycan_{m}"
        output_lines.append(f"{combined_name} = combine {{ {' '.join(block_glycans)} }}")
        combined_glycan_vars.append(combined_name)

    output_lines.append("")  

    output_lines.append(f"combined_final = combine {{ {' '.join(combined_glycan_vars)} }}")

    output_lines.extend([
        "",
        "saveAmberParm combined_final glycam06-beta-chitin-para-fcm.prmtop glycam06-beta-chitin-para-fcm.inpcrd",
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

#if j <= a: 
#    print(f"Processing forward chain polarity for segment {num_segid}")
#    for l in range(1, num_resid):
#        gen.patch(patchname="14BB", targets=[(num_segid, str(l)), (num_segid, str(l + 1))])
#    gen.patch(patchname="14BB", targets=[(num_segid, str(num_resid)), (num_segid, "1")])
#else: 
#    print(f"Processing reverse chain polarity for segment {num_segid}")
#    for m in range(num_resid, 1, -1):
#        gen.patch(patchname="14BB", targets=[(num_segid, str(m)), (num_segid, str(m - 1))])
#    gen.patch(patchname="14BB", targets=[(num_segid, "1"), (num_segid, str(num_resid))])
#count += 1  # Increment count for the next chain