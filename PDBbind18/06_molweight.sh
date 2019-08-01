#!/bin/bash

# List systems where ligand has a molecular mass lower than MAX_MAXX

max_mass=1000 # Dalton
list=pdbbind18.lst

rm -f ${list}

for dataset in "refined" "other"
do
    for dir in $(ls -d ${dataset}/????) 
    do
        system=$(basename ${dir})
        
        # Absolute path
        path=$PWD

        # Input
        sminame=${path}/${dir}/${system}_ligand.smi

        # Heavy atoms molecular weight
        hamw=$(python ../scripts/python/molweight.py ${sminame} | awk '{print $2}')

        if (( $(echo "${hamw} < ${max_mass}" | bc -l) )) 
        then
            echo ${dataset}/${system} >> ${list}
        else
            echo "Discarding ${dataset}/${system} (mw = ${hamw})..."
        fi
    done
done