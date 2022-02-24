#!/bin/bash

# Number of CPUs
n_cpus=12

optdir=$1

if [[ $optdir == "" ]]
then
  echo "OPTDIR must be specified."
  exit
fi

# Paths
source variables/paths

# Datasets
datasets="refined other"

# Working directory
wd=$PWD

# Export variables for GNU parallel
export pdbbind=${pdbbind}
export ddir=${ddir}
export pscripts=${pscripts}
export datasets=${datasets}

# Build single system from docking
makeflex(){
    # Get directory name from argument
    dir=$1

    # PDB name
    system=$(basename ${dir})

    # Combine flexible and rigid part of the receptor
    for dataset in ${datasets}
    do
        rigid=${pdbbind}/${dataset}/${system}/${system}_protein.pdb
        if [ -f ${rigid} ]
        then
            break
        else
            rigid=""
        fi
    done

    echo ${rigid}

    for flex in $(ls -d ${dir}/????_flex-*.pdb)
    do
        idx=$(basename ${flex} .pdb | sed 's#...._flex-##') # Flex index
        python ${pscripts}/makeflex.py ${rigid} ${flex} ${dir}/${system}_protein-${idx}.pdb
    done
}

# Export makeflex for GNU parallel
export -f makeflex

cd ${optdir}

# List all directories
dirs=$(ls -d minimized/*)

# Reconstruct the dataset
parallel -j ${n_cpus} makeflex ::: ${dirs}