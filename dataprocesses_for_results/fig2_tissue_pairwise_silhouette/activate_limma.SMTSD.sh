#!/bin/bash

#SBATCH --job-name=lovar
#SBATCH --array=0-1081
#SBATCH --partition=hive1d
#SBATCH --output=outputs/%A_%a.output

#### RUN TEST: sbatch -J lv_tst -o 0.out -e o.err --wrap "python lowvar.py --pair 0"
# 1080 pairs of STMSD
# 301 pairs of SMTS



# # python lowvar.py --pair $SLURM_ARRAY_TASK_ID --basedon SMTS
# python lowvar.py --pair $SLURM_ARRAY_TASK_ID --basedon SMTSD
# # python lowvar.py --pair $SLURM_ARRAY_TASK_ID --basedon SMTSD


# if [ $SLURM_ARRAY_TASK_ID -lt 302 ]; then
#     python lowvar.py --pair $SLURM_ARRAY_TASK_ID --basedon SMTS
# fi

# before running R, run on the shell:
# module load R/4.3.3 ; module load gcc/13.2.0 && Rscript -e 'install.packages(c("data.table","nlme","optparse"), repos="https://cloud.r-project.org", lib="/lustre1/home/martinmikl/sshanas/R_packages")'
mkdir -p /lustre1/home/martinmikl/sshanas/tmp
export TMPDIR=/lustre1/home/martinmikl/sshanas/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$TMPDIR"
Rscript DE_limma.R --pair=$SLURM_ARRAY_TASK_ID --basedon=SMTSD

# if [ $SLURM_ARRAY_TASK_ID -lt 303 ]; then
    # RscriptR DE_limma.R --pair=$SLURM_ARRAY_TASK_ID --basedon=SMTS
# fi
