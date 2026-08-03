#!/bin/bash

for d in ii.1.SMTSD.resid.div1.csv.gz ii.1.SMTSD.resid.readtpm_g.csv.gz ii.1.SMTS.resid.div1.csv.gz ii.1.SMTS.resid.readtpm_g.csv.gz
    do
        fn=$(echo $d | cut -d '.' -f1-3,5)
        sbatch -J $fn -o outputs/$fn.out -e outputs/$fn.err --mem=20G --wrap "python 6.quantile_norm.py --df $d"
    done
