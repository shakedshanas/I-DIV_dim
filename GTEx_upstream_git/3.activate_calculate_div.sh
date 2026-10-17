#!/bin/bash
#SBATCH --job-name=div
#SBATCH --array=0-10000
#SBATCH --partition=hive1d
#SBATCH --output=outputs/%A_%a.output

ngenes=$(zcat 1.outputs/gene_annotation.csv.gz | wc -l)
smallidx=$(($ngenes - 10000))

python 3.calculate_div.py --geneidx $SLURM_ARRAY_TASK_ID

if [ $SLURM_ARRAY_TASK_ID -le $smallidx ] ; then 
    ACTUAL_JOB_ID=$((SLURM_ARRAY_TASK_ID + 10000))
    python 3.calculate_div.py --geneidx $ACTUAL_JOB_ID 
fi
