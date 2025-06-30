resetpsf
psfcontext reset
package require autoionize
package require solvate
solvate ./cellulose.psf ./cellulose.pdb -s WT -b 1.8 -minmax {{0 0 0} {100.000   100.000  100.000}}
mol addfile solvate.psf
mol addfile solvate.pdb
mol delete top

autoionize -psf solvate.psf -pdb solvate.pdb -sc 0.8 -cation POT -anion CLA  -o system
mol addfile system.psf
mol addfile system.pdb
exit