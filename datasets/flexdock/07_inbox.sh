#!/bin/bash

n_cpus=8
fout="analysis/outbox.lst"
datasets="refined other"

source variables/paths

# Export variables for GNU parallel
export fout=${fout}
export pscripts=${pscripts}

# Score one system
inbox(){
    # Input variables from GNU parallel
    dir=$1
    dataset=$2

    # PDB name
	system=$(basename ${dir})
    
    for ligand in $(ls ${dir}/${system}_ligand-*.pdb)
    do
        # Ligand name and rank
        name=$(basename $ligand .pdb)
        rank=$( echo $name | sed "s#${system}_ligand-##g" )

        # Flexible residues
        flex=${dir}/${system}_flex-${rank}.pdb

        # Extract flexible residues from current protein and original crystal structure
        out=$(python ${pscripts}/inbox.py ${ligand} ${flex})

        ligin=$(echo $out | awk '{print $1}')
        flexin=$(echo $out | awk '{print $2}')

        if [ "$ligin" == "False" ] || [ "$flexin" == "False" ]
        then
            echo "${ligand} ${flex}" >> ${fout}
        fi

    done
}

# Export score for GNU parallel
export -f inbox

# Remove stale output file
rm -f ${fout}

for dataset in ${datasets}
do
    dirs=$(ls -d ${dataset}/????)

    parallel -j ${n_cpus} inbox ::: ${dirs} ::: ${dataset}

done