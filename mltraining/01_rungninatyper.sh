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

        echo $system

        wdir=${typedir}/${system}
        mkdir -p ${wdir}

        # Type docking poses
        for lig in $(ls ${systempath}/${system}_ligand-*.pdb )
        do
            ligname=$(basename ${lig} .pdb)

            echo $ligname

            ${gtyper} ${lig} ${wdir}/${ligname}
        done

        # Type receptor poses
        for rec in $(ls ${systempath}/${system}_protein-*.pdb )
        do
            recname=$(basename ${rec} .pdb)

            echo $recname

            #${gtyper} ${rec} ${wdir}/${recname}
        done

    done    
done
