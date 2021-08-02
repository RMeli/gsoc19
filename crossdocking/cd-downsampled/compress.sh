#!/bin/bash

# --------------------------------------------------
# Compress all files created after docking.
# All files should already be compressed, but
# failing test/cancelled jobs might leave spme files
# uncompressed (which need to be recompressed before
# re-starting the job).
# --------------------------------------------------

#for pocket in $(ls carlos_cd)
for pocket in AOFB
do
    gzip carlos_cd/${pocket}/PDB_Structures/*_default_ensemble_none_flexdist3.5*.pdb
    gzip carlos_cd/${pocket}/PDB_Structures/*_default_ensemble_none_flexdist3.5*.sdf
done