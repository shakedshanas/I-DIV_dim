# %%
import pandas as pd
import numpy as np
import argparse
import warnings
import os
import sys
sys.path.append('/lustre1/home/martinmikl/sshanas/bin')
from div import *
warnings.filterwarnings('ignore')
import time
# %% 
start_time = time.perf_counter()
''' Step 1 - load matadata and define parameters in argspars '''
''' technical infromation about the sample '''
geneids = pd.read_csv(os.getcwd() + '/4.outputs/'  + 'gene_annotation.csv.gz', compression='gzip', usecols=['accession_number'], index_col='accession_number').index
# %%
# import files using argspars
parser = argparse.ArgumentParser(description='Create diversity files')
parser.add_argument('--geneidx', type=int, choices=range(len(geneids)), help='index of genes')
parser.add_argument('--iinput', type=str, help='name of input file')
args = parser.parse_args()
# %%
gid = geneids[args.geneidx]
# %%
''' Step 2 - load data: one gene at a time'''

path_for_isoform_expression_matrix = args.iinput
rows = pd.read_csv(path_for_isoform_expression_matrix, sep=',', compression='gzip', usecols=['accession_number'])
rows_pd_skip = list(rows.loc[rows['accession_number'] != gid].index + 1) # +1 because of the header row in the file

gtex = pd.read_csv(path_for_isoform_expression_matrix, sep=',', compression='gzip', skiprows=rows_pd_skip, index_col = 'accession_number').drop(columns='transcript_id')

# %%    
''' Step 3 - calculate I-DIV '''
df = pd.Series(gtex.apply(lambda x: DIV(arr=x, spreading_percentage=1, tpm_or_count='tpm'), axis='rows'), name=gid).to_frame().T.rename_axis("accession_number", axis="columns")

''' Step 4 - save dataframe to file '''
dfname = 'IDIV'
filepath = os.getcwd() + '/' + args.iinput.replace('.I.', '.IDIV.').replace('.gz','')
df.to_csv(filepath, index=True, mode='a', header = not os.path.exists(filepath))  
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Script executed in {elapsed_time:.4f} seconds")