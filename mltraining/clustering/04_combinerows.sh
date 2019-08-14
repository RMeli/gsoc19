#!/bin/bash

# Paths
spath=../../scripts/python

# Combine similarity measures
python ${spath}/combine_rows_lowmem.py rows.lst -out dist-lsim.pickle
