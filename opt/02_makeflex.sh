#!/bin/bash

# Number of CPUs
n_cpus=10

# Paths
source variables/paths

# List
list=${pdbbind}/test.lst

# Working directory
wd=$PWD

# Export variables for GNU parallel
export pdbbind=${pdbbind}
export pscripts=${pscripts}

# Build single system from docking
makeflex(){
    # Get directory name from argument
    dir=$1

    echo "makeflex ${dir}"

    # PDB name
    system=$(basename ${dir})

    # Define ligand input and output names
    pdbin=${dir}/dock.pdb
    pdbout=${dir}/${system}_ligand-0.pdb

    grep -v "MODEL" ${pdbin} | grep -v "ENDMDL" > ${pdbout}

    # Define flexible residues input and output names
    pdbin=${dir}/flex.pdb
    pdbout=${dir}/${system}_flex-0.pdb

    grep -v "MODEL" ${pdbin} | grep -v "ENDMDL" > ${pdbout}

    # Combine flexible and rigid part of the receptor
    flex=${dir}/${system}_flex-0.pdb
    rigid=${pdbbind}/${dir}/${system}_protein.pdb

    python ${pscripts}/makeflex.py ${rigid} ${flex} ${dir}/${system}_protein-0.pdb
}

# Export makeflex for GNU parallel
export -f makeflex

# List all directories
dirs=$(cat ${list})

# Reconstruct the dataset
parallel -j ${n_cpus} makeflex ::: ${dirs}