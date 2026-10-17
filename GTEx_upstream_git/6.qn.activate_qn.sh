#!/bin/bash
#SBATCH --job-name=qn
#SBATCH --array=0-47
#SBATCH --partition=hive1d
#SBATCH --output=outputs/q.%A_%a.output
#SBATCH --mem=10G

python 6.qn.py --data I --basedon SMTS --tidx $SLURM_ARRAY_TASK_ID
python 6.qn.py --data G --basedon SMTS --tidx $SLURM_ARRAY_TASK_ID
python 6.qn.py --data I --basedon SMTSD --tidx $SLURM_ARRAY_TASK_ID
python 6.qn.py --data G --basedon SMTSD --tidx $SLURM_ARRAY_TASK_ID