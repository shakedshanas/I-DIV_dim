#!/bin/bash
# Run every figure notebook from start to end (from this directory), one SLURM job per notebook.
# The notebooks are independent of each other; they read the outputs of ../GTEx_upstream,
# ../GTEx_pre_downstream and ../Ulicevic_upstream.

mkdir -p source_data supplementary_data tables_main

sbatch -J metrics_comparisons -o metrics_comparisons.out -e metrics_comparisons.err --wrap "jupyter nbconvert --to notebook --execute metrics_comparisons.ipynb --stdout > /dev/null"  # Fig. S1
sbatch -J intro               -o intro.out               -e intro.err               --wrap "jupyter nbconvert --to notebook --execute intro.ipynb --stdout > /dev/null"                # Fig. 1, S2
sbatch -J limma_equivalence   -o limma_equivalence.out   -e limma_equivalence.err   --wrap "jupyter nbconvert --to notebook --execute limma_equivalence.ipynb --stdout > /dev/null"    # Fig. 2, S3
sbatch -J transcriptome_div   -o transcriptome_div.out   -e transcriptome_div.err   --wrap "jupyter nbconvert --to notebook --execute transcriptome_div.ipynb --stdout > /dev/null"    # Fig. 3, S4
sbatch -J lme                 -o lme.out                 -e lme.err                 --wrap "jupyter nbconvert --to notebook --execute lme.ipynb --stdout > /dev/null"                  # Fig. 4, S5
sbatch -J Ulicevic            -o Ulicevic.out            -e Ulicevic.err            --wrap "jupyter nbconvert --to notebook --execute Ulicevic.ipynb --stdout > /dev/null"             # Fig. 5, S6
