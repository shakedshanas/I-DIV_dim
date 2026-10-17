sbatch -J 1.combine_data -o 1.combine_data.out -e 1.combine_data.err --wrap "python 1.combine_data.py"
