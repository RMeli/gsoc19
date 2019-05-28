#!/bin/bash

# Paths
ddir=../../docking     # Docking directory
pdbbind=../../PDBbind18 # PDBbind18 directory

# Working directory
wd=$PWD

for dataset in "test"
do
    for dir in $(ls -d ${ddir}/${dataset}/????)
    do

 	    # PDB name
	    system=$(basename ${dir})

        # Database working directory
        dir=${dataset}/${system}
        mkdir -p ${dir}
        mkdir -p ${dir}/logs

        # Define ligand input and output names
        pdbin=${ddir}/${dataset}/${system}/dock.pdb
        pdbout=${dir}/${system}_ligand-.pdb
        pdbfix=${dir}/lig_fix.pdb

        # Fix END/ENDMDL for OBABEL
        sed 's#^ENDMDL$#<END>#g;s#^END$#<ENDMDL>#g' ${pdbin} > ${pdbfix}
        sed -i 's#<ENDMDL>#ENDMDL#g;s#<END>#END#g' ${pdbfix}

        # Split ligand
        obabel -m -ipdb ${pdbfix} -opdb -O ${pdbout} \
            2>&1 | tee ${dir}/logs/obabel_ligand.log
        rm ${pdbfix}

        # Define flexible residues input and output names
        pdbin=${ddir}/${dataset}/${system}/flex.pdb
        pdbout=${dir}/${system}_flex-.pdb
        pdbfix=${dir}/flex_fix.pdb

        # Fix END/ENDMDL for OBABEL
        sed 's#^ENDMDL$#<END>#g;s#^END$#<ENDMDL>#g' ${pdbin} > ${pdbfix}
        sed -i 's#<ENDMDL>#ENDMDL#g;s#<END>#END#g' ${pdbfix}

        # Split flexible residues
        obabel -m -ipdb ${pdbfix} -opdb -O ${pdbout} \
            2>&1 | tee ${dir}/logs/obabel_flexres.log
        rm ${pdbfix}

        # Combine flexible and rigid part of the receptor
        rigid=${pdbbind}/${dataset}/${system}/${system}_protein.pdb
        for flex in $(ls -d ${dir}/????_flex-*.pdb)
        do
            idx=$(basename ${flex} .pdb | sed 's#...._flex-##') # Flex index

            #python makeflex.py ${rigid} ${flex} ${dir}/${system}_protein-${idx}.pdb
        done
    done
done
