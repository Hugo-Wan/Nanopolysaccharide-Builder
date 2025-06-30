#!/bin/bash
/home/hugowan/software/namd3/psfgen << ENDMOL
topology top_all36_carb.rtf

segment C1 {
  pdb ./1.pdb
  }
patch 14bb C1:1 C1:2
patch 14bb C1:2 C1:3
patch 14bb C1:3 C1:4





regenerate angles dihedrals
guesscoord

writepsf system.psf
ENDMOL
