#!/bin/bash
#SBATCH --job-name=lg_correct
#SBATCH --array=0-1
#SBATCH --partition=hive1d
#SBATCH --output=outputs/%A_%a.output
##SBATCH --dependency=afterok:9265576


lg=$(($SLURM_ARRAY_TASK_ID+1))
sbatch -o outputs/tpm.i.$lg.out --wrap "python 4.correct_data.py --dataname readtpm_g --grouptype i --lgtype $lg"
sbatch -o outputs/div.i.$lg.out --wrap "python 4.correct_data.py --dataname div1 --grouptype i --lgtype $lg"

sbatch -o outputs/tpm.ii.$lg.SMTS.out --wrap "python 4.correct_data.py --dataname readtpm_g --grouptype ii --basedon SMTS --lgtype $lg"
sbatch -o outputs/div.ii.$lg.SMTS.out --wrap "python 4.correct_data.py --dataname div1 --grouptype ii --basedon SMTS --lgtype $lg"

sbatch -o outputs/tpm.ii.$lg.SMTSD.out --wrap "python 4.correct_data.py --dataname readtpm_g --grouptype ii --basedon SMTSD --lgtype $lg"
sbatch -o outputs/div.ii.$lg.SMTSD.out --wrap "python 4.correct_data.py --dataname div1 --grouptype ii --basedon SMTSD --lgtype $lg"
