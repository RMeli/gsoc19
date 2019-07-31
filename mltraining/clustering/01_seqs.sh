#!/bin/bash

source ../variables/paths

# Extract ligands SMILES and receptors sequences
# Format: PDBCODE | SMILES | SEQ
python ${gscripts}/compute_seqs.py \
    --pdbfiles ../../PDBbind18/pdbfiles.lst \
    --out seqs.dat