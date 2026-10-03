#!/bin/bash

# for f in *.csv ; do sbatch -o outputs/$f.3genefilter.out --wrap "python 3.genefilter.py --dfname $f" ; done
# sbatch -o outputs/readcount_g.3genefilter.out --wrap "python 3.genefilter.py --dfname readcount_g.csv.gz" 
sbatch -o outputs/readtpm.3genefilter.out  -J readtpm --wrap "python 3.genefilter.py --dfname outputdir.1/readtpm_g.csv.gz" 
sbatch -o outputs/gene_metadata.3genefilter.out -J gene_metadata --wrap "python 3.genefilter.py --dfname outputdir.1/gene_metadata.csv.gz"
sbatch -o outputs/div1.3genefilter.out -J div1 --wrap "python 3.genefilter.py --dfname outputdir.2/div1.csv"
sbatch -o outputs/nisoforms_expressed.3genefilter.out -J nisoforms_expressed --wrap "python 3.genefilter.py --dfname outputdir.2/nisoforms_expressed.csv"
sbatch -o outputs/total_expression.3genefilter.out -J total_expression--wrap "python 3.genefilter.py --dfname outputdir.2/total_expression.csv"