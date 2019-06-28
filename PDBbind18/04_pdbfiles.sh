#!/bin/bash

# Create PDBFILES file used by gnina clustering.py script
# Used to create cross-validation folds

PDBFILES=pdbfiles.lst

rm -f ${PDBFILES}

for dataset in "refined" "other"
do
    for dir in $(ls -d ${dataset}/????) 
    do
        system=$(basename ${dir})

        # Path relative to home directory
        path=$(echo $PWD | sed "s#${HOME}#~#g")
        
        # Input and output absolute paths
        recname=${path}/${dir}/${system}_protein.pdb
        sminame=${path}/${dir}/${system}_ligand.smi

        echo ${system} ${recname} ${sminame} >> ${PDBFILES}
    done
done