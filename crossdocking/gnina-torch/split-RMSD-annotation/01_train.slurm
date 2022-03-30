#!/bin/bash
#SBATCH --job-name=train-gt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --partition=medium
##SBATCH --partition=devel
#SBATCH --clusters=htc
#SBATCH --array=0-2
#SBATCH --output=slurm/train/%x.%A_%a.out
#SBATCH --time=72:00:00
##SBATCH --time=00:10:00
#SBATCH --account=stat-ecr

hostname
nvidia-smi

which python

# Suffix to distinguish different runs
suffix="loss"

prefix="nc" # cluster/nc
annotation="flex1"
model="default2018"

epochs=500

if [ -n "${suffix}" ];
then
    suffix="-${suffix}"
fi

trainfolder=training/${annotation}/${model}-${prefix}${suffix}
mkdir -p ${trainfolder}

stratify_min=0
stratify_max=0
stratify_step=0
scale_flexpose_loss=1.0

if [ "${suffix}" == "stratified" ]
then
    stratify_min=0
    stratify_max=1
    stratify_step=0.5
fi

if [ "${suffix}" == "loss" ]
then
    scale_flexpose_loss=5.0
fi


python -m gnina.training \
        ../../cd-downsampled/files/SPLIT${annotation}${prefix}train${SLURM_ARRAY_TASK_ID}.types \
        --testfile ../../cd-downsampled/files/SPLIT${annotation}${prefix}test${SLURM_ARRAY_TASK_ID}.types \
        --flexlabel_pos 1 \
        --stratify_pos 1 \
        --stratify_min ${stratify_min} \
        --stratify_max ${stratify_max} \
        --stratify_step ${stratify_step} \
        --scale_flexpose_loss ${scale_flexpose_loss} \
        --data_root ../../cd-downsampled/ \
        --model ${model} \
        --iterations ${epochs} \
        --batch_size 64 \
        --balanced \
        --test_every 10 \
        --checkpoint_every 30 \
        --num_checkpoints 20 \
        --out_dir ${trainfolder} \
        --log_file training${SLURM_ARRAY_TASK_ID}.log \
        --checkpoint_prefix training${SLURM_ARRAY_TASK_ID}