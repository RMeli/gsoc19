#!/bin/bash

# Strip water molecules from PDB files

# List of all PDBbind18 systems (dataset/id)
list=pdbbind18-all.lst

# Number of CPUs (for parallel job)
n_cpus=12

# Strip waters from PDB files in PDBbind
nowat(){
    dir=${1}

    system=$(basename ${dir})

    # Input
    recname=${dir}/${system}_protein.pdb
    recwatname=${dir}/${system}_protein-wat.pdb

    mv ${recname} ${recwatname}

    grep -v "HOH" ${recwatname} > ${recname}
}

export -f nowat

parallel -j ${n_cpus} nowat ::: $(cat ${list})
