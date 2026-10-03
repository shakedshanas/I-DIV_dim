#!/bin/bash

#SBATCH --job-name=corrs
#SBATCH --array=0-10000
#SBATCH --partition=queen
#SBATCH --output=%A_%a.output
##SBATCH --dependency=afterok:11298830

ngenes=$(zcat /lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/gene_metadata.csv.gz | wc -l)
smallidx=$(expr $ngenes - 10000)
echo $SLURM_ARRAY_TASK_ID
python T-DIV.py -i $SLURM_ARRAY_TASK_ID -s all   -f other_and_HK-RP -z rows
python T-DIV.py -i $SLURM_ARRAY_TASK_ID -s SMTSD   -f other_and_HK-RP -z rows
python T-DIV.py -i $SLURM_ARRAY_TASK_ID -s SMTS   -f other_and_HK-RP -z rows
python  mI-DIV.py -i $SLURM_ARRAY_TASK_ID -s all -z rows
python  mI-DIV.py -i $SLURM_ARRAY_TASK_ID -s SMTSD -z rows
python  mI-DIV.py -i $SLURM_ARRAY_TASK_ID -s SMTS  -z rows

if [ $SLURM_ARRAY_TASK_ID -le $smallidx ] ; then 
    ACTUAL_JOB_ID=$((SLURM_ARRAY_TASK_ID + 10000))
    echo $ACTUAL_JOB_ID
    python T-DIV.py -i $ACTUAL_JOB_ID -s all -f other_and_HK-RP -z rows
    python T-DIV.py -i $ACTUAL_JOB_ID -s SMTSD   -f other_and_HK-RP -z rows
    python T-DIV.py -i $ACTUAL_JOB_ID -s SMTS   -f other_and_HK-RP -z rows
    python  mI-DIV.py -i $ACTUAL_JOB_ID -s all -z rows
    python  mI-DIV.py -i $ACTUAL_JOB_ID -s SMTSD -z rows
    python  mI-DIV.py -i $ACTUAL_JOB_ID -s SMTS  -z rows

fi
