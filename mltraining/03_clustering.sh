#!/bin/bash

source variables/paths

n_cpus=8

# Extract ligands SMILES and receptors sequences
# Format: PDBCODE | SMILES | SEQ
python ${gscripts}/compute_seqs.py \
    --pdbfiles ../PDBbind18/pdbfiles.lst \
    --out ${clusterdir}/seqs.dat

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
for row in $(seq 0 ${n_cpus} $((${rows}-1)))
do
    s=$(seq ${row} $((${row} + ${n_cpus} - 1)))
    echo $s
    parallel compute_row ::: ${s}
done

# Combine similarity measures
python ${gscripts}/combine_rows.py ${clusterdir}/row-*
mv matrix.pickle ${clusterdir}/matrix.pickle
