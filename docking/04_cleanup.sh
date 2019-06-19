#!/bin/bash

failed=analysis/failed.dat

for dir in $(cat ${failed})
do
    rm -r ${dir}
done