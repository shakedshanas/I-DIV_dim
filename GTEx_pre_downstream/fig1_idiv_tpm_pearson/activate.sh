#!/bin/bash

#SBATCH --job-name=corrs
#SBATCH --array=0-10000
#SBATCH --partition=hive1d
#SBATCH --output=outputs/%A_%a.output

DIVPATH=../../GTEx_upstream/div1.csv.gz
TPMPATH=../../GTEx_upstream/gene_tpm.csv.gz
NISPATH=../../GTEx_upstream/nisoforms_expressed.csv.gz
GATPATH=../../GTEx_upstream/gene_annotation.csv.gz


ngenes=$(zcat $GATPATH | wc -l)
smallidx=$(expr $ngenes - 10000)

# every array task handles gene i, and gene i + 10000 when it exists
run_gene () {
    # per SMTS tissue: I-DIV ~ TPM, I-DIV ~ NIS
    python pearsonr.py -i $1 -b SMTS -d1 $DIVPATH -d2 $TPMPATH -o div_tpm
    python pearsonr.py -i $1 -b SMTS -d1 $DIVPATH -d2 $NISPATH -o div_nis

    # all samples pooled, sex-specific tissues excluded: I-DIV ~ NIS, NIS ~ TPM, I-DIV ~ TPM
    python pearsonr.py -i $1 -b all_wo_sext -d1 $DIVPATH -d2 $NISPATH -o div_nis
    python pearsonr.py -i $1 -b all_wo_sext -d1 $NISPATH -d2 $TPMPATH -o nis_tpm
    python pearsonr.py -i $1 -b all_wo_sext -d1 $DIVPATH -d2 $TPMPATH -o div_tpm
}

run_gene $SLURM_ARRAY_TASK_ID

if [ $SLURM_ARRAY_TASK_ID -le $smallidx ] ; then
    run_gene $((SLURM_ARRAY_TASK_ID + 10000))
fi
