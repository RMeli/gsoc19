#!/bin/bash

# Compute SMILES strings for all the ligands

for dataset in "other"
do
    for dir in $(ls -d ${dataset}/????) 
    do
        system=$(basename ${dir})
        
        # Input and output names
        ligname=${dir}/${system}_ligand.mol2
        sminame=${dir}/${system}_ligand.smi

        # Convert mol2 to smi and extract SMILES string only
        smiles=$(obabel -imol2 -osmi ${ligname} | awk '{print $1}')

        # Save SMILES to file
        echo ${smiles} > ${sminame}
    done
done