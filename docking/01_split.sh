#!/bin/bash

# Split PDBbind
tdir=splits
mkdir ${tdir}
split -l 100 -d -a 3 ${pdbbind}/pdbbind18.list ${tdir}/split