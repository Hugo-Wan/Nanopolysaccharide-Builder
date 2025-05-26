# Nanopolysaccharide-Builder
Nanopolysaccharide Builder (NPB) is a Python-based package that enables the automatic generation of nanopolysaccharide structures with customizable biochemical topologies. NPB currently supports three model types: cellulose, chitin, and polysaccharide bundles. For each generated structure, corresponding topology files are provided, with support for both GLYCAM06j and CHARMM36 force fields.

# Features
1. Automated File Generation:
Generates geometry coordinate files (PDB format) and topology files (PSF for CHARMM36; PRMTOP for GLYCAM06).
2. Versatile Cross-Section Shapes:
Supports cellulose and chitin models with various cross-sectional geometries, including hexagonal, square, parallelogram, and rectangle. These shapes represent the diverse elementary structures found in biological load-bearing systems, such as plant and fungal cell walls and arthropod cuticles.
3. Customizable Surface Modification:
Enables surface modification at selected crystallographic planes, allowing the construction of user-defined nanopolysaccharide structures.
4. Microfibril Bundle Construction:
Facilitates the building of microfibril bundle models. Currently, only cellulose is supported, with options for parallel or antiparallel arrangements of unit structures.
(Extensible Material Support: Support for hemicellulose, and its linked with lignin and tannins, is under development.)

# Requirements
Python 3.10, Pyside, numpy PyQT5, MDAnalysis, VMD-python, psfgen, Ambertools.
NAMD 3.0/ GROMACS/ AMBER/ LAMMPS

# Installation
We recommend installing NPB on a Linux system, such as Rocky Linux, Ubuntu, or similar distributions. For Windows users, we suggest running NPB inside a virtual machine (using VMware) or through the Windows Subsystem for Linux (WSL).
For Mac users, we have encountered some issues with installing dependencies via conda. We are working on solutions, but for now, we also recommend using a Linux virtual machine to ensure compatibility and a smoother installation experience.
For Linux users, installation can be completed using conda. We recommend installing Miniconda  and creating a new environment for NPB to ensure a clean and safe setup.
````
conda create --name npb python=3.10 (optional)
conda activate npb (optional)
cd /Path/to/npb
python setup.py install
````
