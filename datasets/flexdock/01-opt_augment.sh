#!/bin/bash

# Number of CPUs
n_cpus=10

# Paths
source variables/paths

# List
list=${pdbbind}/lists/test.lst

# Export variables for GNU parallel
export optimised=${optimised}

# Build single system from docking
augment(){
    # Get directory name from argument
    dir=$1
    
    cp ${optimised}/${dir}/*.pdb ${dir}
}

# Export makeflex for GNU parallel
export -f augment

# List all directories
dirs=$(cat ${list})

parallel -j ${n_cpus} augment ::: ${dirs}