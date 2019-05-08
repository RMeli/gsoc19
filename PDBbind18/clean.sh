#!/bin/bash

# Removes useless files from PDBbind:
# - Binding pocket
# - Ligand .sdf file
# - PDF readme
# - Irrelevand INDEX files

for dataset in "refined" "other"
do
    cd ${dataset}

    for dir in $(ls -d ????)
    do
        # Remove ligand .sdf file
        rm ${dir}/*ligand.sdf

        # Remove binding pocket file
        rm ${dir}/*pocket.pdb

        # Remove PDF readme
        rm readme/*.pdf

        # Remove irrelevant INDEX files (NL, PN, PP)
        rm index/INDEX_*_NL.????
        rm index/INDEX_*_PN.????
        rm index/INDEX_*_PP.????
    done

    cd ..
done
