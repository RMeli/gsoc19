#!/bin/bash

source variables/paths

# Split PDBbind
tdir=splits
mkdir ${tdir}
split -l 100 -d -a 3 ${pdbbind}/lists/pdbbind18.lst ${tdir}/split
