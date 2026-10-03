#!/bin/bash
#SBATCH --job-name=div
#SBATCH --array=0-10000
#SBATCH --partition=hive1d
#SBATCH --output=outputs/%A_%a.output

ngenes=$(zcat outputdir.1/gene_metadata.csv.gz | wc -l)
smallidx=$(($ngenes - 10000))

python 2.calculate_div.py --geneidx $SLURM_ARRAY_TASK_ID --datatype tpm

if [ $SLURM_ARRAY_TASK_ID -le $smallidx ] ; then 
    ACTUAL_JOB_ID=$((SLURM_ARRAY_TASK_ID + 10000))
    python 2.calculate_div.py --geneidx $ACTUAL_JOB_ID --datatype tpm
fi
