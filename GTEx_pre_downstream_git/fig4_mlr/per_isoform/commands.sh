sbatch -o PPARD.out -e PPARD.err -J PPARD --wrap "Rscript isoform_age_mlr_limma.R --gene=ENSG00000112033" # PPARD
sbatch -o FAM161B.out -e FAM161B.err -J FAM161B --wrap "Rscript isoform_age_mlr_limma.R --gene=ENSG00000156050" # FAM161B
sbatch -o KNDC1.out -e KNDC1.err -J KNDC1 --wrap "Rscript isoform_age_mlr_limma.R --gene=ENSG00000171798" # KNDC1 
