#!/bin/bash
#SBATCH --job coalescer
#SBATCH --partition=any_cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm/slurm-%A.out
#SBATCH --error=slurm/slurm-%A.err

# Working directory and username
wdir=$(pwd)
user=$(whoami)

date
hostname
echo ${SLURM_JOB_ID}

sroot="${HOME}/git/GNINA-1.0/analysis_scripts"
iroot="${wdir}/gcommands"

#python3 ${sroot}/obrms_calc.py -i ${iroot}/test.txt -d "" --getscores

# CNNrescore
for dir in "cd-carlos-CNNrescore-e8-d3.5" \
	"cd-carlos-CNNrescore-e16-d3.5"
do
	python3 ${sroot}/coalescer.py -r ${wdir}/${dir}/carlos_cd/ -d ${wdir}/dirs.txt -s "" -o ${wdir}/${dir}/${dir}.csv --getscores
done

# CNNnone
for dir in "cd-carlos-CNNnone-e8-d3.5" \
        "cd-carlos-CNNnone-e16-d3.5"
do
	python3 ${sroot}/coalescer.py -r ${wdir}/${dir}/carlos_cd/ -d ${wdir}/dirs.txt -s "" -o ${wdir}/${dir}/${dir}.csv
done
