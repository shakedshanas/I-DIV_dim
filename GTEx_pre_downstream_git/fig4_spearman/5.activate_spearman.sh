#!/bin/bash
#SBATCH --job-name=spr
#SBATCH --array=0-50
#SBATCH --partition=hive7d
#SBATCH --output=outputs/%A_%a.output


python spearman.py --tindex $SLURM_ARRAY_TASK_ID --basedon SMTS
python spearman.py --tindex $SLURM_ARRAY_TASK_ID --basedon SMTSD