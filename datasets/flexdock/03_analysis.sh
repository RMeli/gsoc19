#!/bin/bash

source variables/paths

python ${pscripts}/distrmsd.py \
    analysis/allscores.csv \
    -mr 3 -b 100 \
    -opath analysis/plots
