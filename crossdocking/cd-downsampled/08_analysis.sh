#!/bin/bash

# --------------------------------------------------
# Perform analysis:
#   - RMSD distributions
# --------------------------------------------------

outdir=analysis/plots
mkdir -p ${outdir}

pscripts="../../scripts/python"

python ${pscripts}/distrmsd.py \
    carlos_cd/rmsds.csv \
    -mr 3 -b 0.1 \
    -opath ${outdir} \
    --ligrmsd rmsd --flexrmsd flexrmsd --flexmaxrmsd fmaxrmsd