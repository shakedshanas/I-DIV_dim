import pandas as pd
import numpy as np
import math
import argparse
import warnings
import os


parser = argparse.ArgumentParser(description='Create diversity files')
parser.add_argument('--dfname', type=str, help="path and name of the dataframe that should be filtered")
args = parser.parse_args()


outputpath = os.getcwd() + '/outputdir.3/'

if not os.path.exists( outputpath ):
    os.mkdir(outputpath  )

''' filter genes to have at least 1% of the samples with a value > 0 in the diversity file. 
This is to reduce the number of genes that we have to load in the next steps, which will speed up the process.'''
df = pd.read_csv(os.getcwd() + '/outputdir.2/' + 'div1.csv', index_col='Unnamed: 0')
columnname = df.columns[3]
df = df[df[columnname] != columnname]
df = df.astype(float)

genes_filtered_in_div = (df.apply(lambda x: len(x[x>0]), axis='columns') / df.apply(lambda x: len(x), axis='columns')) 
genes_filtered_in_div = genes_filtered_in_div[genes_filtered_in_div > 0.01].index


del df

''' intersect the two sets of genes to get the final set of genes that we want to keep.'''
# genes_filtered_in = np.intersect1d(genes_filtered_in_div, genes_filtered_in_tpm)

# def filterrows(rows_pd, column_filter_by):
#     ''' GTEx matrices are too big to load all at once, so we filter the rows by gid.
#     This function returns the rows that should be skipped when loading the data.
    
#     Args:
#         rows_pd: row numbers of dataframe.
#         column_filter_by: The column name to filter by (column name of gene IDs).
#         gid: The gene ID to filter for.
    
#     Returns:
#         A list of row indices to skip when loading the data.
#     '''
#     rows_pd_stay = rows_pd.loc[rows_pd[column_filter_by].isin(genes_filtered_in)].index
#     rows_pd_skip = [x + 1  for x in rows_pd.index if x not in rows_pd_stay] 
    
#     return rows_pd_skip

if args.dfname.endswith('.gz'): 
    rows = pd.read_csv(args.dfname, usecols=[0], compression='gzip')
else:
    rows = pd.read_csv(args.dfname, usecols=[0])

colname =  rows.columns[0]

skiprows = list(rows.loc[~rows[colname].isin(genes_filtered_in_div)].index + 1) # +1 because of the header row in the file
del genes_filtered_in_div

# skiprows = filterrows(rows, column_filter_by=colname)
df = pd.read_csv(args.dfname, skiprows=skiprows)
df.set_index(colname, inplace=True)
if args.dfname == "outputdir.1/gene_metadata.csv.gz":
    div = pd.read_csv(os.getcwd() + '/outputdir.2/' + 'div1.csv', index_col='Unnamed: 0')
    genes_filtered_in_div = (div.apply(lambda x: len(x[x>0]), axis='columns') / div.apply(lambda x: len(x), axis='columns'))
    genes_filtered_in_div.name = 'DIVpercentile_gt0'
    df = df.join(genes_filtered_in_div)

    ''' filter genes to have at least 1% of the samples with a value > 0 in the diversity file. 
    This is to reduce the number of genes that we have to load in the next steps, which will speed up the process.'''

    tpm = pd.read_csv(os.getcwd() + '/outputdir.1/' + 'readtpm_g.csv.gz', index_col='accession_number', compression='gzip') 
    columnname = tpm.columns[3]
    tpm = tpm[tpm[columnname] != columnname]
    tpm = tpm.astype(float)

    genes_filtered_in_tpm = (tpm.apply(lambda x: len(x[x>0]), axis='columns') / tpm.apply(lambda x: len(x), axis='columns')) 
    genes_filtered_in_tpm.name = 'TPMpercentile_gt0'
    df = df.join(genes_filtered_in_tpm)

filename = args.dfname.split('/')[-1].replace('.csv', '').replace('.gz', '') + ".csv.gz"
df.to_csv(outputpath + filename, index=True, compression='gzip')
