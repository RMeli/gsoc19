#!/bin/bash

sroot="${HOME}/git/GNINA-1.0/analysis_scripts"

mkdir -p plots

for cnn in "rescore" "none"
do
	Cs=""
	Ns=""

	if [ $cnn == "rescore" ]
	then
		# CNN data in CSV file
		usecols="0 2 7"
	elif [ $cnn == "none" ]
	then
		# No CNN data in CSV file
		usecols="0 2 4"
	fi

	# Generate filenames to compare and corresponding names
	for e in 8 16
	do
		name="cd-carlos-CNN${cnn}-e${e}-d3.5"

		Cs="${Cs} ${name}/${name}.csv"
		Ns="${Ns} CNN${cnn}-E${e}"
	done

	echo $Cs
	echo $Ns

	python3 ${sroot}/generate_RMSD_graphs.py \
                        -C ${Cs} -N ${Ns} \
                        -U 2850 --usecols ${usecols} \
                        --figname plots/RMSD-CNN${cnn} --line_graph --bar_graph
done