#!/bin/bash
#SBATCH --job-name=silh
#SBATCH --array=0-1080
#SBATCH --partition=hive1d
#SBATCH --output=outputs/process.%A_%a.output

# PCA + silhouette score per tissue pair, from the limma results of step 2.
# --pair is 0-indexed here (the script adds 1 to find the limma file of that pair).
# Each task appends one row to <basedon>_limma/overview_cluster_results.csv.
python 3.limma_process_res.py --pair $SLURM_ARRAY_TASK_ID --basedon SMTSD
if [ $SLURM_ARRAY_TASK_ID -lt 300 ] ; then
    python 3.limma_process_res.py --pair $SLURM_ARRAY_TASK_ID --basedon SMTS
fi
