#!/bin/bash

#SBATCH --job-name=lovar
#SBATCH --array=0-1081
#SBATCH --partition=hive1d
#SBATCH --output=outputs/SMTSD.%A_%a.output

# 1080 pairs of STMSD
# 300 pairs of SMTS

# before running R, run on the shell:
# module load R/4.3.3 ; module load gcc/13.2.0 && Rscript -e 'install.packages(c("data.table","nlme","optparse"), repos="https://cloud.r-project.org", lib="/lustre1/home/martinmikl/sshanas/R_packages")'
mkdir -p /lustre1/home/martinmikl/sshanas/tmp
export TMPDIR=/lustre1/home/martinmikl/sshanas/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$TMPDIR"
Rscript DE_limma.R --pair=$SLURM_ARRAY_TASK_ID --basedon=SMTSD
