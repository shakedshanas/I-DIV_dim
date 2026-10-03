#!/bin/bash
#SBATCH --job-name=corrs
#SBATCH --array=0-10000
#SBATCH --partition=hive1d
#SBATCH --output=outputs/spct_corrs_%A_%a.output


ngenes=$(zcat /lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/gene_metadata.csv.gz | wc -l)
smallidx=$(expr $ngenes - 10000)

python SMTS_wo_sext/div_nis.py -i $SLURM_ARRAY_TASK_ID
python SMTSD_brain_regions/div_nis.py -i $SLURM_ARRAY_TASK_ID

python SMTS_wo_sext/div_tpm.py -i $SLURM_ARRAY_TASK_ID
python SMTSD_brain_regions/div_tpm.py -i $SLURM_ARRAY_TASK_ID

python SMTS_wo_sext/tpm_nis.py -i $SLURM_ARRAY_TASK_ID
python SMTSD_brain_regions/tpm_nis.py -i $SLURM_ARRAY_TASK_ID

if [ $SLURM_ARRAY_TASK_ID -le $smallidx ] ; then 
    ACTUAL_JOB_ID=$((SLURM_ARRAY_TASK_ID + 10000))

    python SMTS_wo_sext/div_nis.py -i $ACTUAL_JOB_ID
    python SMTSD_brain_regions/div_nis.py -i $ACTUAL_JOB_ID

    python SMTS_wo_sext/div_tpm.py -i $ACTUAL_JOB_ID
    python SMTSD_brain_regions/div_tpm.py -i $ACTUAL_JOB_ID

    python SMTS_wo_sext/tpm_nis.py -i $ACTUAL_JOB_ID
    python SMTSD_brain_regions/tpm_nis.py -i $ACTUAL_JOB_ID

fi