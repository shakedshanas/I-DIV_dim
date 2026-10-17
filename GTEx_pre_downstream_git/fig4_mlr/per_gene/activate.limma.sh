#!/bin/bash
#SBATCH --job-name=limma
#SBATCH --array=0-50
#SBATCH --partition=hive7d
#SBATCH --output=outputs/%A_%a.limma.SMTSD.output

# module load R/4.3.3   
# module load gcc/13.2.0 


export TMPDIR=/lustre1/home/martinmikl/sshanas/tmp/$SLURM_JOB_ID   # or your scratch path
mkdir -p "$TMPDIR"
# Rscript limma.R --tindex=$SLURM_ARRAY_TASK_ID --basedon=SMTSD
Rscript limma_age_per_tissue.R --tindex=$SLURM_ARRAY_TASK_ID --basedon=SMTSD
Rscript limma_age_per_tissue.R --tindex=$SLURM_ARRAY_TASK_ID --basedon=SMTS
# Rscript limma_SMTS.R --tindex=$SLURM_ARRAY_TASK_ID --basedon=SMTS
rm -rf "$TMPDIR"
