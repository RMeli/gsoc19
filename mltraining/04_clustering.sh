#!/bin/bash

n_cpus=8 # Number of CPUs for parallel calculations

# Paths
source variables/paths

# Export variables for GNU parallel
export rows=$(wc ${clusterdir}/seqs.dat | awk '{print $1}')
export gscripts=${gscripts}
export clusterdir=${clusterdir}

# Function to compute a single row of the similarity matrix
compute_row(){
    row=$1

    echo ">>> Computing row $((${row}+1)) of ${rows}..."

    python ${gscripts}/compute_row.py \
        --pdbseqs ${clusterdir}/seqs.dat \
        --row ${row} \
        --out ${clusterdir}/row-${row}
}

# Export function for GNU parallel
export -f compute_row

# Compute similarity of one system with all other systems
# Parallelised over all systems
s=$(seq 0 $((${rows}-1))) # All rows
parallel -j ${n_cpus} compute_row ::: ${s}

# Combine similarity measures
python ${gscripts}/combine_rows.py ${clusterdir}/row-*
mv matrix.pickle ${clusterdir}/matrix.pickle
