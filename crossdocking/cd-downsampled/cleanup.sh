#!/bin/bash

# --------------------------------------------------
# Remove all files created after docking
# --------------------------------------------------

for pocket in $(ls carlos_cd)
#for pocket in AOFB
do
    # Cleanup 03_crustal.slurm
    #rm carlos_cd/${pocket}/PDB_Structures/*_default_ensemble_none_flexdist3.5*_p*.*
    #rm carlos_cd/${pocket}/PDB_Structures/*_default_ensemble_none_flexdist3.5*_crys*.*
    
    # Cleanup 05_spyrmsd.slurm
    rm carlos_cd/${pocket}/PDB_Structures/*.rmsds
    find carlos_cd/${pocket}/PDB_Structures/ -type d -name '*.d' -exec rm -r {} +
done