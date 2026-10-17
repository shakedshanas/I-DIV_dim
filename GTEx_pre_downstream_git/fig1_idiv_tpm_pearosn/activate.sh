#!/bin/bash

#SBATCH --job-name=corrs
#SBATCH --array=0-10000
#SBATCH --partition=hive1d
#SBATCH --output=outputs/%A_%a.output

DIVPATH=/lustre1/home/martinmikl/sshanas/outsourcedata/GTEx_upstream_git/div1.csv.gz
TPMPATH=/lustre1/home/martinmikl/sshanas/outsourcedata/GTEx_upstream_git/gene_tpm.csv.gz
NISPATH=/lustre1/home/martinmikl/sshanas/outsourcedata/GTEx_upstream_git/nisoforms_expressed.csv.gz
GATPATH=/lustre1/home/martinmikl/sshanas/outsourcedata/GTEx_upstream_git/gene_annotation.csv.gz


ngenes=$(zcat $GATPATH | wc -l)
smallidx=$(expr $ngenes - 10000)


# python pearsonr.py -i $SLURM_ARRAY_TASK_ID -b all -d1 $DIVPATH -d2 $TPMPATH -o div_tpm
python pearsonr.py -i $SLURM_ARRAY_TASK_ID -b SMTS -d1 $DIVPATH -d2 $TPMPATH -o div_tpm
# python pearsonr.py -i $SLURM_ARRAY_TASK_ID -b SMTSD -d1 $DIVPATH -d2 $TPMPATH -o div_tpm

# python pearsonr.py -i $SLURM_ARRAY_TASK_ID -b all -d1 $DIVPATH -d2 $NISPATH -o div_nis
python pearsonr.py -i $SLURM_ARRAY_TASK_ID -b SMTS -d1 $DIVPATH -d2 $NISPATH -o div_nis
# python pearsonr.py -i $SLURM_ARRAY_TASK_ID -b SMTSD -d1 $DIVPATH -d2 $NISPATH -o div_nis


if [ $SLURM_ARRAY_TASK_ID -le $smallidx ] ; then 
    ACTUAL_JOB_ID=$((SLURM_ARRAY_TASK_ID + 10000))

    # python pearsonr.py -i $ACTUAL_JOB_ID -b all -d1 $DIVPATH -d2 $TPMPATH -o div_tpm
    python pearsonr.py -i $ACTUAL_JOB_ID -b SMTS -d1 $DIVPATH -d2 $TPMPATH -o div_tpm
    # python pearsonr.py -i $ACTUAL_JOB_ID -b SMTSD -d1 $DIVPATH -d2 $TPMPATH -o div_tpm

    # python pearsonr.py -i $ACTUAL_JOB_ID -b all -d1 $DIVPATH -d2 $NISPATH -o div_nis
    python pearsonr.py -i $ACTUAL_JOB_ID -b SMTS -d1 $DIVPATH -d2 $NISPATH -o div_nis
    # python pearsonr.py -i $ACTUAL_JOB_ID -b SMTSD -d1 $DIVPATH -d2 $NISPATH -o div_nis

fi
