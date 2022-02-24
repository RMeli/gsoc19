#!/bin/bash

invalid=analysis/invalid.lst

# Remove systems failed
for dir in $(cat ${invalid})
do
    rm -r ${dir}
done