#!/bin/bash

# Generate file with gnina commands
# Conform to standard docking scripts for analysis

flexdist=3.5
seed=42

# Data root
dataroot=/net/pulsar/home/koes/paf46/Research/gnina1.0
ifile=${dataroot}/carlos_cd_ds_cd.txt

mkdir -p gcommands

for exhaustiveness in 8 16
do
	for cnnscore in rescore none
	do
		ofile="gcommands/gcommands-CNN${cnnscore}-e${exhaustiveness}.txt"
		rm -rf ${ofile}

		root="cd-carlos-CNN${cnnscore}-e${exhaustiveness}-d3.5"

		for idx in $(seq 1 8153)
		do
			isystems=$(sed -n "${idx}p" ${ifile})

			# Ligand and receptor names (with pocket)
			ligand=$(echo ${isystems} | cut -f2 -d " ")
			receptor=$(echo ${isystems} | cut -f1 -d " ")

			# Path to ligand and receptor files
			ligpath=${dataroot}/${ligand}
			recpath=${dataroot}/${receptor}

			# Pocket, ligand and receptor names
			pocket=$(dirname ${ligand})
			lig=$(basename ${ligand} .sdf)
			rec=$(basename ${receptor} .pdb)

			# Ligand and receptor output names
			ligout=${pocket}/flexlig-${lig}-${rec}.sdf
			recout=${pocket}/flexrec-${lig}-${rec}.pdb

			# Command for flexible docking
			# dataroot and .gz output extension added for compatibility with analysis scripts
			echo "gnina -r i${dataroot}/${receptor} -l ${dataroot}/${ligand} --flexdist_ligand ${dataroot}/${ligand} --flexdist ${flexdist} --autobox_ligand ${dataroot}/${ligand} --exhaustiveness ${exha
ustiveness} --num_modes 10 --seed ${seed} --cnn_scoring ${cnnscore} --no_gpu --cpu 1 --out ${root}/${ligout}.gz --out_flex ${root}/${recout}.gz" >> ${ofile}

		done
	done
done
