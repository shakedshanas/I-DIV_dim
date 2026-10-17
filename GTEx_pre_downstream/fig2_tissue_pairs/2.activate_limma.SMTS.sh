#!/bin/bash

#SBATCH --job-name=lovar
# DE_limma.R pairs are 1-indexed: 25 SMTS -> 300 pairs
#SBATCH --array=1-300
#SBATCH --partition=hive1d
#SBATCH --output=outputs/SMTS.%A_%a.output

# 1081 pairs of SMTSD
# 300 pairs of SMTS


# before running R, run on the shell:
# module load R/4.3.3 ; module load gcc/13.2.0 && Rscript -e 'install.packages(c("data.table","nlme","optparse"), repos="https://cloud.r-project.org", lib="<your R library>")'
mkdir -p ${SCRATCH_DIR:-/tmp}
export TMPDIR=${SCRATCH_DIR:-/tmp}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$TMPDIR"
Rscript DE_limma.R --pair=$SLURM_ARRAY_TASK_ID --basedon=SMTS