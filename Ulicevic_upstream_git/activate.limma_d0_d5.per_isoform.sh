#!/bin/bash

# run the three following commands before running: sbatch 5.activate_lme.sh
module load R/4.3.3 
module load gcc/13.2.0 
# Rscript -e 'install.packages(c("data.table","nlme","optparse"), repos="https://cloud.r-project.org", lib="/lustre1/home/martinmikl/sshanas/R_packages")'

sbatch -o DLL3.out -e DLL3.err -J DLL3 --wrap "Rscript limma_d0_d5.per_isoform.R --gene=ENSG00000090932" # DLL3
sbatch -o DHX15.out -e DHX15.err -J DHX15 --wrap "Rscript limma_d0_d5.per_isoform.R --gene=ENSG00000109606" # DHX15
sbatch -o RNPS1.out -e RNPS1.err -J RNPS1 --wrap "Rscript limma_d0_d5.per_isoform.R --gene=ENSG00000205937" # RNPS1
