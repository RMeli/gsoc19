#!/bin/bash

noflexfile=../../docking/analysis/noflex.dat

source variables/paths

for noflex in $(cat ${noflexfile})
do
    system=$(basename ${noflex})

    # Remove flexible residues file
    rm ${noflex}/${system}_flex.pdb

    # Add original receptor
    cp ${pdbbind}/${noflex}/${system}_protein.pdb ${noflex}
done