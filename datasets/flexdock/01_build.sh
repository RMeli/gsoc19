#!/bin/bash

# Number of CPUs
n_cpus=10

# List
list=lists/refined.lst

# Paths
source variables/paths

# Working directory
wd=$PWD

# Export variables for GNU parallel
export pdbbind=${pdbbind}
export ddir=${ddir}
export pscripts=${pscripts}

# Build single system from docking
makeflex(){
    # Get directory name from argument
    dir=$1

    echo "makeflex ${dir}"

    # PDB name
    system=$(basename ${dir})

    # Database working directory
    mkdir -p ${dir}
    mkdir -p ${dir}/logs

    # Define ligand input and output names
    pdbin=${ddir}/${dir}/dock.pdb
    pdbout=${dir}/${system}_ligand-.pdb
    pdbfix=${dir}/lig_fix.pdb

    # Fix END/ENDMDL for OBABEL
    grep -v "^END$" ${pdbin} > ${pdbfix}
    echo "END" >> ${pdbfix}

    # Split ligand
    obabel -m -ipdb ${pdbfix} -opdb -O ${pdbout} \
        2>&1 | tee ${dir}/logs/obabel_ligand.log
    mv ${pdbfix} ${dir}/${system}_ligand.pdb

    # Define flexible residues input and output names
    pdbin=${ddir}/${dir}/flex.pdb
    pdbout=${dir}/${system}_flex-.pdb
    pdbfix=${dir}/flex_fix.pdb

    # Fix END/ENDMDL for OBABEL
    grep -v "^END$" ${pdbin} > ${pdbfix}
    echo "END" >> ${pdbfix}

    # Split flexible residues
    obabel -m -ipdb ${pdbfix} -opdb -O ${pdbout} \
        2>&1 | tee ${dir}/logs/obabel_flexres.log
    mv ${pdbfix} ${dir}/${system}_flex.pdb

    # Combine flexible and rigid part of the receptor
    rigid=${pdbbind}/${dir}/${system}_protein.pdb
    for flex in $(ls -d ${dir}/????_flex-*.pdb)
    do
        idx=$(basename ${flex} .pdb | sed 's#...._flex-##') # Flex index

        python ${pscripts}/makeflex.py ${rigid} ${flex} ${dir}/${system}_protein-${idx}.pdb
    done
}

# Export makeflex for GNU parallel
export -f makeflex

# List all directories
dirs=$(cat ${list})

# Reconstruct the dataset
parallel -j ${n_cpus} makeflex ::: ${dirs}