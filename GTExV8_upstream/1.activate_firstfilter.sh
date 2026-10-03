#!/bin/bash
#SBATCH --job-name=filterdata
#SBATCH --array=0-3
#SBATCH --partition=hive1d
#SBATCH --output=outputs/%A_%a.output

if [ ! -d outputs ]; then mkdir outputs ; fi
if [ ! -d outputdir.1 ]; then mkdir outputdir.1 ; fi
if [ ! -d outputdir.2 ]; then mkdir outputdir.2 ; fi
if [ ! -d outputdir.3 ]; then mkdir outputdir.3 ; fi
if [ ! -d outputdir.4 ]; then mkdir outputdir.4 ; fi
if [ ! -d outputdir.5 ]; then mkdir outputdir.5 ; fi

python 1.firstfilter.py --data $SLURM_ARRAY_TASK_ID