#!/bin/bash

source variables/paths

python ${pscripts}/distrmsd.py \
    analysis/allscores.csv \
    -mr 3 -b 0.1 \
    -opath analysis/plots

#python ${pscripts}/rmsdthreshold.py \
#    analysis/allscores.csv \
#    -opath analysis/plots