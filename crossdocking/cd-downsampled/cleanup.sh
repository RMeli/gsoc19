#!/bin/bash

# --------------------------------------------------
# Remove all files created after docking
# --------------------------------------------------

for pocket in $(ls carlos_cd)
do
    #rm carlos_cd/${pocket}/PDB_Structures/*_default_ensemble_none_flexdist3.5*_p*.*
    #rm carlos_cd/${pocket}/PDB_Structures/*_default_ensemble_none_flexdist3.5*_crys*.*
    rm carlos_cd/${pocket}/PDB_Structures/*_default_ensemble_none_flexdist3.5*.rmsds
done