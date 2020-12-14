#!/bin/bash
#SBATCH --job obrms
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
python3 ${sroot}/obrms_calc.py -i ${iroot}/gcommands-CNNrescore-e8.txt -d "" --getscores
python3 ${sroot}/obrms_calc.py -i ${iroot}/gcommands-CNNrescore-e16.txt -d "" --getscores

# CNNnone
python3 ${sroot}/obrms_calc.py -i ${iroot}/gcommands-CNNnone-e8.txt -d ""
python3 ${sroot}/obrms_calc.py -i ${iroot}/gcommands-CNNnone-e16.txt -d ""
