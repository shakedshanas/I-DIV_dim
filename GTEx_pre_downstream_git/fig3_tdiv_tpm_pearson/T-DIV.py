import pandas as pd
import numpy as np
import scipy.stats as st
import os
import argparse
import math
import gzip
import sys
sys.path.append('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025')
from div import *
''' 
edited on August 15: added the function filterrows and incoporated it to the script. Goal: make the code faster & less computationally intesive. 
edited on January 11th 2026: change of input data, transcriptome diversity will be by gene group (named other_only or R_only)

'''

# %%
#### PATHS ####
metadata_path = '/lustre1/home/martinmikl/sshanas/outsourcedata/GTEx_upstream_git/metadata.csv.gz'
tpm_path = '/lustre1/home/martinmikl/sshanas/outsourcedata/GTEx_upstream_git/gene_tpm.zs.csv.gz'
###############
parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', type=int, help='index of genes')
parser.add_argument('-s', type=str, choices=['SMTS', 'SMTSD','all'], help='samples grouped or all samples')
args = parser.parse_args()

# %%
ngenes_percentile_threshold = 1

def filterrows(rows_pd):
    rows_pd_stay = rows_pd[rows_pd['accession_number'].str.contains(gid)].index
    rows_pd_skip = [x + 1  for x in rows_pd.index if x not in rows_pd_stay] 
    
    return rows_pd_skip

# %%

metadata = pd.read_csv(metadata_path, index_col='SAMPID', compression='gzip')
metadata = metadata[metadata[f'gt0'] >= metadata[f'percentile_ngenes_' + str(ngenes_percentile_threshold)]]

rows_pd = pd.read_csv(f'readtpm_g.zscore.{args.z}.csv.gz',compression='gzip', usecols=['accession_number'])
gid = rows_pd.loc[int(args.i), 'accession_number']

rows_pd_skip = filterrows(rows_pd)
tpm = pd.read_csv(tpm_path, skiprows=rows_pd_skip, compression='gzip', index_col='accession_number')

rawtpm = pd.read_csv(tpm_path, compression='gzip', index_col=0, skiprows=rows_pd_skip, index_col='accession_number')
del rows_pd, rows_pd_skip

# pearsonn_df_heg = pd.DataFrame()
pearsonn_df = pd.DataFrame()

if args.s == 'all':
    metadata = metadata.join( pd.Series(tpm.loc[gid], name = 'zscore_logtpm_' + gid)  )
    metadata = metadata.join( pd.Series(rawtpm.loc[gid], name = 'rawtpm_' + gid)  )

    metadata = metadata[metadata['rawtpm_' + gid] > 0]
    if len(metadata) >= 10:
        pearsonn_df.loc[gid, 'coeff'], pearsonn_df.loc[gid, 'pval'] = st.pearsonr(metadata[f'div'], metadata['zscore_logtpm_' + gid])
    else:
        # pearsonn_df_heg.loc[gid, 'coeff'], pearsonn_df_heg.loc[gid, 'pval'] = 0,1
        pearsonn_df.loc[gid, 'coeff'], pearsonn_df.loc[gid, 'pval'] = 0,1

elif args.s != 'all':
    tissues_list = metadata[args.s].value_counts()[metadata[args.s].value_counts() > 10].index
    # samplegroups = metadata[args.s].unique()
    for smpl in tissues_list:
        samplelist = metadata[metadata[args.s] == smpl ].index
        temp_metadata = metadata.loc[samplelist].join( tpm.loc[gid] ).join( pd.Series(rawtpm.loc[gid], name = 'rawtpm_' + gid)  )
        temp_metadata = temp_metadata[temp_metadata['rawtpm_' + gid] > 0]
        if len(temp_metadata) >= 10:
            pearsonn_df.loc[gid, smpl + '_coeff'], pearsonn_df.loc[gid, smpl + '_pval'] = st.pearsonr(temp_metadata[f'div'], temp_metadata[gid])
        else:
            # pearsonn_df_heg.loc[gid, smpl + '_coeff'], pearsonn_df_heg.loc[gid, smpl + '_pval'] = 0,1
            pearsonn_df.loc[gid, smpl + '_coeff'], pearsonn_df.loc[gid, smpl + '_pval'] = 0,1
   
fn = '.'.join(['pearson','T-DIV',args.s,'csv'])
pearsonn_df.to_csv(fn, index=True, mode='a', header = not os.path.exists(fn))  