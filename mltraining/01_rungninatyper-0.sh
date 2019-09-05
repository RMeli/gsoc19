#!/bin/bash

n_cpus=6

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

    echo $sys

    systempath=${database}/${sys}

    system=$(basename ${systempath})

    wdir=${typedir}/${system}
    mkdir -p ${wdir}

    # Type docking poses
    lig=${systempath}/${system}_ligand-0.pdb
    ligname=$(basename ${lig} .pdb)

    ${gtyper} ${lig} ${wdir}/${ligname}.gninatypes

    # Type receptor poses
    rec=${systempath}/${system}_protein-0.pdb
    recname=$(basename ${rec} .pdb)

    ${gtyper} ${rec} ${wdir}/${recname}.gninatypes
}

# Export makeflex for GNU parallel
export -f type

# List all directories
sys=$(cat ${list})

# Reconstruct the dataset
parallel -j ${n_cpus} type ::: ${sys}