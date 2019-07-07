#!/bin/bash

n_nodes=12
n_cores=24

n=$(( ${n_nodes} * ${n_cores} ))

rows=$(cat seqs.dat | wc -l)

echo "Number of rows: $rows"

rm -f mapping.dat

mkdir -p maps

row=0
while [ $row -lt $rows ]
do
    for i in $(seq 0 $(($n - 1)))
    do
	echo $i $row >> mapping.dat

	row=$(( $row + 1 ))

	if [ $row -ge $rows ]
 	then
	    break
	fi
    done
done

split -d -l ${n} mapping.dat maps/map

rm mapping.dat