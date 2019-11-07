#!/bin/bash

root=${PWD}
datadir=flexcd

datapath=${root}/${datadir}

cd ${datapath}

for pocket in $(ls -d *)
do
    cd ${pocket}

    for pdbfile in $(ls *.pdb)
    do
        # Fix END/ENDMDL for OBABEL
        grep -v "^END$" ${pdbfile} > tmp.pdb
        echo "END" >> tmp.pdb
        mv tmp.pdb ${pdbfile}
    done

    cd ${datapath}
done