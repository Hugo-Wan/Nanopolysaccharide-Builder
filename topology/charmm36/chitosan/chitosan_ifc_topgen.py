from psfgen import PsfGen
import sys
import numpy as np
import math
import warnings
import os
import glob
import json

#int(sys.argv[1])     
with open('config.json', 'r') as f:
    config = json.load(f)
main_folder_path = config['main_folder_path']
structure_input_file = os.path.join(main_folder_path, 'charmm36-chitosan-ifcm.pdb')
topology_input_file = os.path.join(main_folder_path, 'topology', 'charmm36', 'top_all36_carb.rtf')
gen = PsfGen(output="/dev/null")  
gen.read_topology(topology_input_file)
gen.add_segment(segid="CARB", pdbfile=structure_input_file)
gen.read_coords(segid="CARB", filename=structure_input_file)
#print(len(gen.get_resids(segid="CARB")))
num_resid  = len(gen.get_resids(segid="CARB"))
###from reduced end to non-reduced end
num_residues = num_resid   
for i in range(1, num_residues):
    gen.patch(patchname="14BB", targets=[("CARB", str(i)), ("CARB", str(i + 1))])
gen.patch(patchname="14BB", targets=[("CARB", str(num_residues)), ("CARB", "1")])


gen.guess_coords()
gen.regenerate_angles()
gen.regenerate_dihedrals()

gen.write_psf(filename=os.path.join(main_folder_path, "charmm36-chitosan-ifcm.psf"))
#gen.write_pdb(filename=os.path.join(main_folder_path, "system.pdb")) 