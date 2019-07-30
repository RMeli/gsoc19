#!/bin/bash

# Strip water molecules from PDB files

for dataset in "test"
do
    for dir in $(ls -d ${dataset}/????) 
    do
        system=$(basename ${dir})
        
        # Absolute path
        path=$PWD

        # Input
        recname=${path}/${dir}/${system}_protein.pdb
        recwatname=${path}/${dir}/${system}_protein-wat.pdb

        mv ${recname} ${recwatname}

        grep -v "HOH" ${recwatname} > ${recname}
    done
done