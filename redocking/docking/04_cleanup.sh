#!/bin/bash

failed=analysis/failed.lst
noflex=analysis/noflex.lst
noconf=analysis/noconf.lst

# Remove systems failed
for dir in $(cat ${failed})
do
    rm -r ${dir}
done

# Remove systems without flexible residues
for dir in $(cat ${noflex})
do
    rm -r ${dir}
done

# Remove systems for which no conformation was found
for dir in $(cat ${noconf})
do
    rm -r ${dir}
done