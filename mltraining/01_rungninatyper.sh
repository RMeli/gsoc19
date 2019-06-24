#!/bin/bash

source variables/paths

mkdir -p ${typedir}

for dataset in "test"
do
    datadir=${database}/${dataset}
    echo $datadir
    for systempath in $(ls -d ${datadir}/????)
    do
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

    done    
done
