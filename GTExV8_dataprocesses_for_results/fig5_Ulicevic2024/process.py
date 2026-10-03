 # %%
import pandas as pd
import os
import math
from numpy import array
import warnings
warnings.filterwarnings("ignore")
import sys
sys.path.append('/Volumes/data/shaked/scripts/functions')
from div import *

# %%
gene_dict = {file.split('.')[0][-6:]: pd.read_csv(file, compression='gzip', sep='\t', usecols=['gene_id', 'TPM']) for file in os.listdir() if file.endswith('genes.results.gz')}
gene_df = gene_dict['day0_1'].set_index('gene_id').drop(columns='TPM')
for k in sorted(gene_dict.keys()):
    gene_df = gene_df.join(gene_dict[k].set_index('gene_id')).rename(columns={'TPM': k})
gene_df.to_csv('tpm.csv', index=True)   
del gene_df

isoform_dict = {file.split('.')[0][-6:]: pd.read_csv(file, compression='gzip', sep='\t', usecols=['transcript_id', 'gene_id', 'TPM']) for file in os.listdir() if file.endswith('isoforms.results.gz')}
isoform_df = isoform_dict['day0_1'][['transcript_id', 'gene_id']]
for k in sorted(isoform_dict.keys()):
    isoform_df = isoform_df.merge(isoform_dict[k], on=['transcript_id', 'gene_id']).rename(columns={'TPM': k})
isoform_df.drop(columns='transcript_id', inplace=True)
del isoform_dict
 
dataframes = {}
dataframes['evenness0'] = pd.DataFrame()
dataframes['evenness1'] = pd.DataFrame()
evenness1_min = pd.DataFrame()
dataframes['sum_pergene'] = pd.DataFrame()
dataframes['sum_nisoforms'] = pd.DataFrame()


# Iterate through groups and perform operations
for gene, group in isoform_df.groupby(by='gene_id'):
    # group.drop(columns='accession_number', inplace=True)
    
        #Those lines are for calculating our evenness with different spreading
    dataframes['evenness0'][gene]=group.set_index('gene_id').apply(lambda x: DIV(x.values, 0, 'tpm'), axis=0)                      #spreading 0
    dataframes['evenness1'][gene]=group.set_index('gene_id').apply(lambda x: DIV(x.values, 1, 'tpm', ceil=False), axis=0)           #spreading 1
    
    #This line compute shannon_eveness_min
    evenness1_min[gene]=group.set_index('gene_id').apply(lambda x: DIV_min(x.values, 1, 'tpm', ceil=True), axis=0)   #spreading 1, minimum diversity

    #This line computes # of isoforms per gene and read counts per gene
    dataframes['sum_pergene'][gene] = group.set_index('gene_id').sum()
    # dataframes['sum_nisoforms'][gene] = group.apply(lambda x: x.sum(), axis=0)
    dataframes['sum_nisoforms'][gene] = group.set_index('gene_id').apply(lambda x: len(x.values[x.values>0]), axis=0)
     
# Those lines create the normalized diversity
dataframes['evenness1_norm'] = dataframes['evenness1']-evenness1_min
dataframes['evenness1_norm2'] = dataframes['evenness1_norm']/(1-evenness1_min)
del evenness1_min
for k in dataframes.keys():
    dataframes[k].to_csv(k + '.csv', index=True)