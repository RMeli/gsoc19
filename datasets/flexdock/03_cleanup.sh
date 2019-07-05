#!/bin/bash

problematic=analysis/problematic.lst

# Remove systems failed
for dir in $(cat ${problematic})
do
    rm -r ${dir}
done