#!/bin/bash

sbatch -o outputs/4.genefilter.out -e outputs/4.genefilter.err  -J 4.genefilter --wrap "python 4.genefilter.py" 
