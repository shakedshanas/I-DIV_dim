#!/bin/bash

sbatch -J filt -o outputs/2.filtgenes.out -e outputs/2.filtgenes.err --wrap "python 2.firstfilter.py --data gene_tpm"
sbatch -J filt -o outputs/2.filttanscript.out -e outputs/2.filttanscript.err --wrap "python 2.firstfilter.py --data transcript_tpm"