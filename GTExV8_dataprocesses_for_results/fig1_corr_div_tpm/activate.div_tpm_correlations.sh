#!/bin/bash

#SBATCH --job-name=corrs
#SBATCH --array=0-10000
#SBATCH --partition=hive1d
#SBATCH --output=outputs/%A_%a.output

ngenes=$(zcat /lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/gene_metadata.csv.gz | wc -l)
smallidx=$(expr $ngenes - 10000)

DIVPATH=/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/div1.csv.gz
TPMPATH=/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/readtpm_g.csv.gz
NISPATH=/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/nisoforms_expressed.csv.gz

# python div_tpm_correlations.py -i $SLURM_ARRAY_TASK_ID -b all -d1 $DIVPATH -d2 $TPMPATH -o div_tpm
# python div_tpm_correlations.py -i $SLURM_ARRAY_TASK_ID -b SMTS -d1 $DIVPATH -d2 $TPMPATH -o div_tpm
# python div_tpm_correlations.py -i $SLURM_ARRAY_TASK_ID -b SMTSD -d1 $DIVPATH -d2 $TPMPATH -o div_tpm

# python div_tpm_correlations.py -i $SLURM_ARRAY_TASK_ID -b all -d1 $DIVPATH -d2 $NISPATH -o div_nis
# python div_tpm_correlations.py -i $SLURM_ARRAY_TASK_ID -b SMTS -d1 $DIVPATH -d2 $NISPATH -o div_nis
# python div_tpm_correlations.py -i $SLURM_ARRAY_TASK_ID -b SMTSD -d1 $DIVPATH -d2 $NISPATH -o div_nis

python tpm_nis_correlations.py -i $SLURM_ARRAY_TASK_ID -b all -o tpm_nis
python tpm_nis_correlations.py -i $SLURM_ARRAY_TASK_ID -b SMTS -o tpm_nis
python tpm_nis_correlations.py -i $SLURM_ARRAY_TASK_ID -b SMTSD -o tpm_nis
 

if [ $SLURM_ARRAY_TASK_ID -le $smallidx ] ; then 
    ACTUAL_JOB_ID=$((SLURM_ARRAY_TASK_ID + 10000))

    # python div_tpm_correlations.py -i $ACTUAL_JOB_ID -b all -d1 $DIVPATH -d2 $TPMPATH -o div_tpm
    # python div_tpm_correlations.py -i $ACTUAL_JOB_ID -b SMTS -d1 $DIVPATH -d2 $TPMPATH -o div_tpm
    # python div_tpm_correlations.py -i $ACTUAL_JOB_ID -b SMTSD -d1 $DIVPATH -d2 $TPMPATH -o div_tpm

    # python div_tpm_correlations.py -i $ACTUAL_JOB_ID -b all -d1 $DIVPATH -d2 $NISPATH -o div_nis
    # python div_tpm_correlations.py -i $ACTUAL_JOB_ID -b SMTS -d1 $DIVPATH -d2 $NISPATH -o div_nis
    # python div_tpm_correlations.py -i $ACTUAL_JOB_ID -b SMTSD -d1 $DIVPATH -d2 $NISPATH -o div_nis

    python tpm_nis_correlations.py -i $ACTUAL_JOB_ID -b all -o tpm_nis
    python tpm_nis_correlations.py -i $ACTUAL_JOB_ID -b SMTS -o tpm_nis
    python tpm_nis_correlations.py -i $ACTUAL_JOB_ID -b SMTSD -o tpm_nis
fi
