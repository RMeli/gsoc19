#!/bin/bash
#SBATCH --job cd-test-3.0
#SBATCH --partition=any_gpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task 4
#SBATCH --gres=gpu:1
#SBATCH -C C5
#SBATCH --array=0-3
#SBATCH --output=slurm/slurm_%A_%a.out
#SBATCH --error=slurm/slurm_%A_%a.err

flexdist=3.0
autoboxadd=4

# Root output folder
EXPERIMENT="docking-d${flexdist}"

# Working directory and username
work_dir=$(pwd)
user=$(whoami)

# Create output folders
outdir="${work_dir}/${EXPERIMENT}"
mkdir -p ${outdir}/docking
mkdir -p ${outdir}/info

# Load CUDA
module load cuda/10.2

# Job directory on node
job_dir="${user}_${SLURM_JOB_ID}.dcb.private.net"
mkdir /scr/$job_dir
cd /scr/$job_dir

# Determine list of ligand/receptor pairs to run
# Based on slurm array job id
printf -v splitfile "cd-m2_%03d" ${SLURM_ARRAY_TASK_ID}
echo $splitfile

# Save list of ligand/receptor pairs name
echo $splitfile | tee $SLURM_ARRAY_TASK_ID.txt

# Copy list with ligand/receptor pairs
rsync -a ${work_dir}/lists/${splitfile} .

# Starting time
date | tee -a $SLURM_ARRAY_TASK_ID.txt

# Data root
dataroot=${work_dir}/data/CrossDocked
echo "Data root: ${dataroot}" | tee -a $SLURM_ARRAY_TASK_ID.txt

# Loop over ligand/receptor pairs
while read line
do
    # Ligand and receptor names (with pocket)
    ligand=$(echo ${line} | cut -f1 -d " ").gz
    receptor=$(echo ${line} | cut -f2 -d " ").gz

    echo "Ligand: ${ligand}" | tee -a $SLURM_ARRAY_TASK_ID.txt
    echo "Receptor: ${receptor}" | tee -a $SLURM_ARRAY_TASK_ID.txt

    # Path to ligand and receptor files
    ligpath=${dataroot}/${ligand}
    recpath=${dataroot}/${receptor}

    # Pocket, ligand and receptor names
    pocket=$(dirname ${ligand})
    lig=$(basename ${ligand} .sdf.gz)
    rec=$(basename ${receptor} .pdb.gz)

    echo "Pocket: ${pocket}" | tee -a $SLURM_ARRAY_TASK_ID.txt
    echo "Lig: ${lig}" | tee -a $SLURM_ARRAY_TASK_ID.txt
    echo "Rec: ${rec}" | tee -a $SLURM_ARRAY_TASK_ID.txt

    # Create dir for pocket
    mkdir -p ${pocket}

    # Copy ligand and receptor files on node
    cp ${ligpath} ${ligand}
    cp ${recpath} ${receptor}

    # Ligand and receptor output names
    ligout=${pocket}/flexlig-${lig}-${rec}.sdf
    recout=${pocket}/flexrec-${lig}-${rec}.pdb

    # Run flexible docking
    gnina -r ${receptor} -l ${ligand} \
        --flexdist_ligand ${ligand} --flexdist ${flexdist} \
	--autobox_ligand ${ligand} --autobox_add ${autoboxadd} \
	--exhaustiveness 8 --num_modes 20 --cpu ${SLURM_TASKS_PER_NODE} \
        --out ${ligout} --out_flex ${recout} \
	2>&1 | tee ${pocket}/$(basename ${lig} .sdf)-$(basename ${rec} .pdb)-gnina.log

    # Compress output
    gzip ${ligout}
    gzip ${recout}

    # Remove input
    rm ${ligand} ${receptor}

    # Copy docking results on head node
    rsync -ar ${pocket} ${outdir}/docking/

    #break

done < ${splitfile}

# Node and end time
hostname | tee -a $SLURM_ARRAY_TASK_ID.txt
date | tee -a $SLURM_ARRAY_TASK_ID.txt

# Copy info on head node
rsync -ra *.txt ${outdir}/info