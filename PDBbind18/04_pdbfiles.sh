#!/bin/bash

# Create PDBFILES file used by gnina clustering.py script
# Used to create cross-validation folds

PDBFILES=pdbfiles.lst

rm -f ${PDBFILES}

for dataset in "test"
do
    for dir in $(ls -d ${dataset}/????) 
    do
        system=$(basename ${dir})

        # Absolute path
        path=$PWD
        
        # Input and output absolute paths
        recname=${path}/${dir}/${system}_protein.pdb
        sminame=${path}/${dir}/${system}_ligand.smi

        echo ${system} ${recname} ${sminame} >> ${PDBFILES}
    done
done