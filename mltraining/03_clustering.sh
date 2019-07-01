#!/bin/bash

source variables/paths

# Extract ligands SMILES and receptors sequences
# Format: PDBCODE | SMILES | SEQ
python ${gscripts}/compute_seqs.py \
    --pdbfiles ../PDBbind18/pdbfiles.lst \
    --out ${clusterdir}/seqs.dat

# Compute similarity of one system with all other systems
# TODO: Parallelize this loop
rows=$(wc ${typedir}/seqs.dat | awk '{print $1}')
for row in $(seq 0 $((${rows}-1)))
do
    echo ">>> Computing row $((${row}+1)) of ${rows}..."

    python ${gscripts}/compute_row.py \
        --pdbseqs ${typedir}/seqs.dat \
        --row ${row} \
        --out ${clusterdir}/row-${row}
done

# Combine similarity measures
python ${gscripts}/combine_rows.py ${clusterdir}/row-*
mv matrix.pickle ${clusterdir}/matrix.pickle

# Create folds
python ${gscripts}/clustering.py \
    --cpickle  ${clusterdir}/matrix.pickle \
    --input ${typedir}/all.types \
    --output ${typedir}/all

