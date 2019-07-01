#!/bin/bash

source variables/paths

python ${gscripts}/compute_seqs.py \
    --pdbfiles ../PDBbind18/pdbfiles.lst \
    --out ${typedir}/seqs.dat