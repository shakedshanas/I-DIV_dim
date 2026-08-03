
# %%
import pandas as pd
import numpy as np
import argparse
import warnings
import os
from div import *
warnings.filterwarnings('ignore')

# %% 
''' Step 1 - load matadata and define parameters in argspars '''
''' technical infromation about the sample '''

inputpath = os.getcwd() + '/outputdir.1/'
outputpath = os.getcwd() + '/outputdir.2/'


if not os.path.exists( outputpath ):
    os.mkdir(outputpath  )
geneids = pd.read_csv(inputpath + 'gene_metadata.csv.gz', compression='gzip', usecols=['accession_number'], index_col='accession_number').index
# %%
# import files using argspars
parser = argparse.ArgumentParser(description='Create diversity files')
parser.add_argument('--geneidx', type=int, choices=range(len(geneids)), help='index of genes')
parser.add_argument('--datatype', type=str, choices=['tpm', 'count'], help="Which file to work with: 'TPM' or 'TRC'")
args = parser.parse_args()


# %%
gid = geneids[args.geneidx]
# %%
''' Step 2 - load data: one gene at a time'''

rows = pd.read_csv(inputpath + f'read{args.datatype}_t.csv.gz', sep=',', compression='gzip', usecols=['accession_number'])
# rows_pd_skip = filterrows(rows_pd=rows, column_filter_by='accession_number', gid=gid)

# rows_pd_skip = list(rows.loc[~rows['accession_number'] == gid].index + 1) # +1 because of the header row in the file
rows_pd_skip = list(rows.loc[rows['accession_number'] != gid].index + 1) # +1 because of the header row in the file

gtex = pd.read_csv(inputpath + f'read{args.datatype}_t.csv.gz', sep=',', skiprows=rows_pd_skip, compression='gzip')
gtex.set_index('accession_number', inplace=True)
# %%    
''' Step 3 - process data '''

dataframes = {}
dataframes['div1'] = pd.Series(gtex.apply(lambda x: DIV(arr=x, spreading_percentage=1, tpm_or_count=args.datatype), axis='rows'), name=gid).to_frame().T
dataframes['total_expression'] = pd.Series(gtex.sum(axis='rows'), name=gid).to_frame().T
dataframes['nisoforms_expressed'] = pd.Series(gtex.apply(lambda x: len(x.values[x.values>0]), axis='rows'), name=gid).to_frame().T

# dataframes['div0'] = pd.Series(gtex.apply(lambda x: DIV(arr=x, spreading_percentage=0, tpm_or_count=args.datatype), axis='rows'), name=gid).to_frame().T
# evenness1_min=pd.Series(gtex.apply(lambda x: DIV_min(arr=x, spreading_percentage=1, tpm_or_count=args.datatype, ceil=True), axis='rows'), name=gid).to_frame().T
    
# normalized diversity
# dataframes['div1_norm'] = dataframes['div1']-evenness1_min
# dataframes['div1_norm2'] = dataframes['div1_norm']/(1-evenness1_min)
# del evenness1_min


''' Step 4 - save dataframes to files '''
for dfname, df in dataframes.items():
    filepath = outputpath + f"{dfname}.csv"
    df.to_csv(filepath, index=True, mode='a', header = not os.path.exists(filepath))        