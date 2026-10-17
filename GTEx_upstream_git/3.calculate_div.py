
# %%
import pandas as pd
import numpy as np
import argparse
import warnings
import os
import sys
sys.path.append('/lustre1/home/martinmikl/sshanas/bin')
from div import *
import time
warnings.filterwarnings('ignore')

# %% 
''' Step 1 - load matadata and define parameters in argspars '''
''' technical infromation about the sample '''
start_time = time.perf_counter()
inputpath = os.getcwd() + '/1.outputs/'
# outputpath = os.getcwd() + '/outputdir.2/'



geneids = pd.read_csv(os.getcwd() + '/1.outputs/' + 'gene_annotation.csv.gz', compression='gzip', usecols=['accession_number'], index_col='accession_number').index
# %%
# import files using argspars
parser = argparse.ArgumentParser(description='aggregate transcript files to I-DIV or number of isoforms per gene')
parser.add_argument('--geneidx', type=int, choices=range(len(geneids)), help='index of genes')
args = parser.parse_args()


# %%
gid = geneids[args.geneidx]
# %%
''' Step 2 - load data: one gene at a time'''

rows = pd.read_csv(os.getcwd() + '/2.outputs/' + 'transcript_tpm.csv.gz', sep=',', compression='gzip', usecols=['accession_number'])
rows_pd_skip = list(rows.loc[rows['accession_number'] != gid].index + 1) # +1 because of the header row in the file

gtex = pd.read_csv(os.getcwd() + '/2.outputs/' + 'transcript_tpm.csv.gz', sep=',', skiprows=rows_pd_skip, compression='gzip').drop(columns='transcript_id')
gtex.set_index('accession_number', inplace=True)
# %%    
''' Step 3 - process data '''

dataframes = {}
dataframes['div1'] = pd.Series(gtex.apply(lambda x: DIV(arr=x, spreading_percentage=1, tpm_or_count='tpm'), axis='rows'), name=gid).to_frame().T.rename_axis("accession_number", axis="columns")
dataframes['nisoforms_expressed'] = pd.Series(gtex.apply(lambda x: len(x.values[x.values>0]), axis='rows'), name=gid).to_frame().T.rename_axis("accession_number", axis="columns")

''' Step 4 - save dataframes to files '''
for dfname, df in dataframes.items():
    filepath = os.getcwd() + '/3.outputs/' + f"{dfname}.csv"
    df.to_csv(filepath, index=True, mode='a', header = not os.path.exists(filepath))        

end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Script executed in {elapsed_time:.4f} seconds")