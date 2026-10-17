#!/bin/bash
sbatch -J TDIV -o outputs/5.calculate_transcriptome_div.out -e outputs/5.calculate_transcriptome_div.err --wrap "python 5.calculate_transcriptome_div.py"
