#!/bin/bash

n_cpus=8

source variables/paths

list=${database}/analysis/valid.lst

mkdir -p ${typedir}

# Export variables for GNU parallel
export typedir=${typedir}
export database=${database}
export pdbbind=${pdbbind}
export gtyper=${gtyper}

type(){
    # Get system name from argument
    sys=$1

    systempath=${database}/${sys}

    system=$(basename ${systempath})

    wdir=${typedir}/${system}
    mkdir -p ${wdir}

    # Type docking poses
    for lig in $(ls ${systempath}/${system}_ligand-*.pdb )
    do
        ligname=$(basename ${lig} .pdb)

        ${gtyper} ${lig} ${wdir}/${ligname}.gninatypes
    done

    # Type receptor poses
    for rec in $(ls ${systempath}/${system}_protein-*.pdb )
    do
        recname=$(basename ${rec} .pdb)

        ${gtyper} ${rec} ${wdir}/${recname}.gninatypes
    done 
}

# Export makeflex for GNU parallel
export -f type

# List all directories
sys=$(cat ${list})

# Reconstruct the dataset
parallel -j ${n_cpus} type ::: ${sys}