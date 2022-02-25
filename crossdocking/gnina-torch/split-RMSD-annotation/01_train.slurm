#!/bin/bash
#SBATCH --job-name=train
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --partition=devel
#SBATCH --clusters=htc
#SBATCH --array=0-2
#SBATCH --output=slurm/train/%x.%A_%a.out
#SBATCH --time=00:10:00
#SBATCH --account=stat-ecr

hostname
nvidia-smi

which python

prefix="nc" # cluster/nc
annotation="flex1"

epochs=10

if [ ${annotation} = "" ]
then
    trainfolder=training/ligand/${modelname}
else
    trainfolder=training/${annotation}/${modelname}
fi
mkdir -p ${trainfolder}

python -m gnina.training \
        ../../cd-downsampled/files/SPLIT${annotation}${prefix}train${SLURM_ARRAY_TASK_ID}.types \
        --testfile ../../cd-downsampled/files/SPLIT${annotation}${prefix}test${SLURM_ARRAY_TASK_ID}.types \
        --flexlabel_pos 1 \
        --data_root ../../cd-downsampled/ \
        --ligmolcache lig.molcache2 \
        --recmolcache rec.molcache2 \
        --iterations ${epochs} \
        --batch_size 64 \
        --balanced \
        --test_every 2 \
        --checkpoint_every 2 \
        --num_checkpoints 2 \
        --out_dir ${trainfolder} \
        --log_file training${SLURM_ARRAY_TASK_ID}.log \
        --checkpoint_prefix training${SLURM_ARRAY_TASK_ID}