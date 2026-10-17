sbatch -J 1.make_metadata_files -o outputs/1.make_metadata_files.out -e outputs/1.make_metadata_files.err --wrap "python 1.make_metadata_files.py"
