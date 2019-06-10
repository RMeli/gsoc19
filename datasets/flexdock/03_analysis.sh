#!/bin/bash

source variables/paths

python3.6 ${pscripts}/distrmsd.py allscores.csv -mr 2 -b 100 -opath analysis