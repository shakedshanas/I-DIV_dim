#!/bin/bash
#SBATCH --job-name=transcriptome_div
#SBATCH --array=0-6
#SBATCH --partition=hive1d
#SBATCH --output=outputs/%A_%a.output

if [ ! -d outputdir.5 ]; then mkdir outputdir.5 ; fi
python 5.calculate_transcriptome_div.py --gf $SLURM_ARRAY_TASK_ID
