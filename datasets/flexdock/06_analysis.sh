#!/bin/bash

source variables/paths

python ${pscripts}/distrmsd.py \
    analysis/allscores.csv \
    -mr 5 -b 100 \
    -opath analysis/plots

python ${pscripts}/rmsdthreshold.py \
    analysis/allscores.csv \
    -opath analysis/plots