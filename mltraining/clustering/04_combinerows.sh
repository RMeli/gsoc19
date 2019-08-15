#!/bin/bash

# Paths
spath=../../scripts/python

# Combine similarity measures
python ${gscripts}/combine_rows_lowmem.py rows/row-* --out dist-lsim.pickle
