#!/bin/bash

# Paths
source ../variables/paths

# Combine similarity measures
python ${gscripts}/combine_rows_lowmem.py rows/row-* --out dist-lsim.pickle