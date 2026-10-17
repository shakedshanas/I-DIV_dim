#!/bin/bash

# isoform-level limma, day 0 vs day 5, for the three example genes of Fig. 5 (run after 2.calculate_div.py)
# R 4.3.3 with the packages in ../R_requirements.txt; on our cluster R is loaded with:
module load R/4.3.3 
module load gcc/13.2.0 

sbatch -o DLL3.out -e DLL3.err -J DLL3 --wrap "Rscript limma_d0_d5.per_isoform.R --gene=ENSG00000090932" # DLL3
sbatch -o DHX15.out -e DHX15.err -J DHX15 --wrap "Rscript limma_d0_d5.per_isoform.R --gene=ENSG00000109606" # DHX15
sbatch -o RNPS1.out -e RNPS1.err -J RNPS1 --wrap "Rscript limma_d0_d5.per_isoform.R --gene=ENSG00000205937" # RNPS1
