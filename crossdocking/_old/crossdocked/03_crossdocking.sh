#!/bin/bash

#SBATCH --job name=crossdocking
#SBATCH --partition=any_cpu
#SBATCH --nodes=1
##SBATCH --array=0-5
#SBATCH --array=0-862
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4
#SBATCH -o slurm/slurm_%A_%a.out
#SBATCH -e slurm/slurm_%A_%a.err

# Load smina module
# Compiled by dkoes
module load smina/smina

# Working directory and username
work_dir=$(pwd)
user=$(whoami)

# Job directory on node
job_dir="${user}_${SLURM_JOB_ID}.dcb.private.net"

mkdir /scr/$job_dir
cd /scr/$job_dir

# Determine list of ligand/receptor pairs to run
# Based on slurm array job id
printf -v splitfile "cd-m2_%03d" ${SLURM_ARRAY_TASK_ID}
echo $splitfile

# Save list of ligand/receptor pairs name
echo $splitfile>$SLURM_ARRAY_TASK_ID.txt

# Copy list with ligand/receptor pairs
rsync -a ${work_dir}/lists/${splitfile} .

# Starting time
date>>$SLURM_ARRAY_TASK_ID.txt

# Data root
dataroot=${work_dir}/data/CrossDocked
echo "Data root: ${dataroot}"

# Loop over ligand/receptor pairs
while read line
do 
    # Ligand and receptor names (with pocket)
    ligand=$(echo ${line} | cut -f1 -d " ").gz
    receptor=$(echo ${line} | cut -f2 -d " ").gz

    echo "Ligand: ${ligand}"
    echo "Receptor: ${receptor}"

    # Path to ligand and receptor files
    ligpath=${dataroot}/${ligand}
    recpath=${dataroot}/${receptor}
	
    # Pocket, ligand and receptor names
    pocket=$(dirname ${ligand})
    lig=$(basename ${ligand} .sdf.gz)
    rec=$(basename ${receptor} .pdb.gz)

    echo "Pocket: ${pocket}"
    echo "Lig: ${lig}"
    echo "Rec: ${rec}"

    # Create dir for pocket
    mkdir -p ${pocket} 

    # Copy ligand and receptor files on node
    cp ${ligpath} ${ligand}
    cp ${recpath} ${receptor}

    # Ligand and receptor output names
    ligout=${pocket}/flexlig-${lig}.sdf
    recout=${pocket}/flexrec-${rec}.pdb

    # Run flexible docking
    smina -r ${receptor} -l ${ligand} \
        --flexdist_ligand ${ligand} --flexdist 3 \
	--autobox_ligand ${ligand} --autobox_add 8 \
	--exhaustiveness 8 --num_modes 20 --cpu ${SLURM_TASKS_PER_NODE} \
        --out ${ligout} --out_flex ${recout} \
	2>&1 | tee ${pocket}/$(basename ${lig} .sdf)-$(basename ${rec} .pdb)-smina.log

    # Compress output
    gzip ${ligout}
    gzip ${recout}

    # Remove input
    rm ${ligand} ${receptor}

    # Copy docking results on head node
    rsync -ar ${pocket} ${work_dir}/docking/

    #break
    
done < ${splitfile}

# Node and end time
hostname>>$SLURM_ARRAY_TASK_ID.txt
date>>$SLURM_ARRAY_TASK_ID.txt

# Copy info on head node
rsync -ra *.txt ${work_dir}/info
