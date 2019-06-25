#!/bin/bash

source variables/paths

python3.6 ${pscripts}/distrmsd.py \
    analysis/allscores.csv \
    -mr 5 -b 100 \
    -opath analysis/plots