#!/bin/bash

# congrats, you are done processing the GTEx data. 
# next step is to clean, after this move to the directory ../GTEx_pre_downstream_git/ to process the data with statistics
# as a preperation for the downstream results including: figures 1-4

mv 4.outputs/* . && rm -r 4.outputs
rm -r 1.outputs/ 2.outputs/ 3.outputs/ 
rm outputs/*
gzip *.csv
rm -r qn_data/
