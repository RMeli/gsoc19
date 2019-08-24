#!/bin/bash

# Number of CPUs
n_cpus=10

# Paths
source variables/paths

# Working directory
wd=$PWD

# Export variables for GNU parallel
export flexdock=${flexdock}

# Build single system from docking
augment(){
    # Get directory name from argument
    dir=$1
    
    cp ${dir}/*.pdb ${flexdock}/${dir}
}

# Export makeflex for GNU parallel
export -f augment

# Reconstruct the dataset
for dataset in "other"
do
    dirs=$(ls -d ${dataset}/*)
    parallel -j ${n_cpus} augment ::: ${dirs}
done