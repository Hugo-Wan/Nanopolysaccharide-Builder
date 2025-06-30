from psfgen import PsfGen
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
warnings.filterwarnings("ignore", category=UserWarning)


with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
structure_input_file = os.path.join(main_folder_path, 'charmm36-c1b-4-bundles-icm.pdb')
topology_input_file = os.path.join(main_folder_path, 'topology', 'charmm36', 'top_all36_carb.rtf')
gen = PsfGen(output="/dev/null")  

u = mda.Universe(structure_input_file)
structure = molecule.load("pdb", structure_input_file)
all_atoms = atomsel("all", molid=structure)
all_fragment = set(all_atoms.fragment)
num_fragment = len(all_fragment)
print(num_fragment)

end=num_fragment


def get_fragment_indices(mol_id, start_fragment):
    structure_index = []
    fragment_selection = atomsel(f"fragment {start_fragment}", molid=mol_id)
    indices = fragment_selection.index
    structure_index.extend(indices)
    return structure_index



for i in range(0, num_fragment):
    selected_indices = get_fragment_indices(structure, int(i))
    ndx_i = selected_indices[0]
    ndx_j = selected_indices[-1]
    j=i+1
    print(f"process chain {j}")
    selected_fragment = u.select_atoms(f"index {ndx_i} : {ndx_j}")
    selected_fragment.write(f"cellulose_temp_{j}.pdb") 
    input_pdb = os.path.join(main_folder_path, f'cellulose_temp_{j}.pdb')
    gen.read_topology(topology_input_file)
    gen.add_segment(segid=str(j), pdbfile=input_pdb)
    num_segid = str(j)
    num_resid = len(gen.get_resids(segid=num_segid))
    print(f"resid number is {num_resid}")
    for l in range(1, num_resid):
        gen.patch(patchname="14BB", targets=[(num_segid, str(l)), (num_segid, str(l + 1))])
    gen.patch(patchname="14BB", targets=[(num_segid, str(num_resid)), (num_segid, str(1))])


gen.guess_coords()
gen.regenerate_angles()
gen.regenerate_dihedrals()
gen.write_psf(filename=os.path.join(main_folder_path, "charmm36-c1b-4-bundles-icm.psf"))


for temp_file in glob.glob("*temp*.pdb"):
    os.remove(temp_file)

