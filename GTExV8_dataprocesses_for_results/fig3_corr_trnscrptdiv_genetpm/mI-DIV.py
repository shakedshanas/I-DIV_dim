# %%
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

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', type=int, help='index of genes') # choices=range(15705)
parser.add_argument('-s', type=str, choices=['SMTS', 'SMTSD','all'], help='samples grouped or all samples')
parser.add_argument('-z', type=str, choices=['columns', 'rows'], help='zscore on columns or rows')
# parser.add_argument('-f', type=str, choices=['RP','HK-RP','no_filter','other_and_HK', 'other_only', 'RP-HK', 'RP_U_HK'], default='', help='gene_filtering_for_metrics')

args = parser.parse_args()

# %%


def filterrows(rows_pd):
    rows_pd_stay = rows_pd[rows_pd['accession_number'].str.contains(gid)].index
    rows_pd_skip = [x + 1  for x in rows_pd.index if x not in rows_pd_stay] 
    
    return rows_pd_skip

# %%
ngenes_percentile_threshold = 1
metadata = pd.read_csv(f'/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.5/no_filter/metadata.csv.gz', index_col='SAMPID', compression='gzip')
# metadata = metadata[metadata[f'gt0'] >= metadata[f'percentile_ngenes_' + str(ngenes_percentile_threshold)]] # 16 July 2026 got commented
    
rows_pd = pd.read_csv(f'readtpm_g.zscore.{args.z}.csv.gz',compression='gzip', usecols=['accession_number'])
gid = rows_pd.loc[int(args.i), 'accession_number']

rows_pd_skip = filterrows(rows_pd)

rawtpm = pd.read_csv(f'/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/readtpm_g.csv.gz', compression='gzip', index_col=0, skiprows=rows_pd_skip)
tpm = pd.read_csv(f'readtpm_g.zscore.{args.z}.csv.gz', skiprows=rows_pd_skip, compression='gzip', index_col = ['accession_number'])    
del rows_pd, rows_pd_skip

# pearsonn_df_heg = pd.DataFrame()
pearsonn_df = pd.DataFrame()

if args.s == 'all':
    metadata = metadata.join( pd.Series(tpm.loc[gid], name = 'zscore_logtpm_' + gid)  )
    metadata = metadata.join( pd.Series(rawtpm.loc[gid], name = 'rawtpm_' + gid)  )
    metadata = metadata[metadata['rawtpm_' + gid] > 0]
    if len(metadata) >= 10:
        # pearsonn_df_heg.loc[gid, 'coeff'], pearsonn_df_heg.loc[gid, 'pval'] = st.pearsonr(metadata[f'ngenes_{ngenes_percentile_threshold}percentile_HEG_meanisofdiv'], metadata['zscore_logtpm_' + gid])
        pearsonn_df.loc[gid, 'coeff'], pearsonn_df.loc[gid, 'pval'] = st.pearsonr(metadata[f'meanisofdiv'], metadata['zscore_logtpm_' + gid])
    else:
        # pearsonn_df_heg.loc[gid, 'coeff'], pearsonn_df_heg.loc[gid, 'pval'] = 0,1
        pearsonn_df.loc[gid, 'coeff'], pearsonn_df.loc[gid, 'pval'] = 0,1

elif args.s != 'all':
    tissues_list = metadata[args.s].value_counts()[metadata[args.s].value_counts() > 10].index
    samplegroups = metadata[args.s].unique()
    for smpl in tissues_list:
        samplelist = metadata[metadata[args.s] == smpl ].index
        temp_metadata = metadata.loc[samplelist].join( tpm.loc[gid] ).join( pd.Series(rawtpm.loc[gid], name = 'rawtpm_' + gid)  )
        temp_metadata = temp_metadata[temp_metadata['rawtpm_' + gid] > 0]
        if len(temp_metadata) >= 10:
            # pearsonn_df_heg.loc[gid, smpl + '_coeff'], pearsonn_df_heg.loc[gid, smpl + '_pval'] = st.pearsonr(temp_metadata[f'ngenes_{ngenes_percentile_threshold}percentile_HEG_meanisofdiv'], temp_metadata[gid])
            pearsonn_df.loc[gid, smpl + '_coeff'], pearsonn_df.loc[gid, smpl + '_pval'] = st.pearsonr(temp_metadata[f'meanisofdiv'], temp_metadata[gid])
        else:
            # pearsonn_df_heg.loc[gid, smpl + '_coeff'], pearsonn_df_heg.loc[gid, smpl + '_pval'] = 0,1
            pearsonn_df.loc[gid, smpl + '_coeff'], pearsonn_df.loc[gid, smpl + '_pval'] = 0,1

# fn = '.'.join(['pearson',str(ngenes_percentile_threshold) + '_nreads_heg','mI-DIV','no_filter',args.s,args.z,'csv'])
# pearsonn_df_heg.to_csv(fn, index=True, mode='a', header = not os.path.exists(fn))              
fn = '.'.join(['pearson','entire','mI-DIV','no_filter',args.s,args.z,'csv'])
pearsonn_df.to_csv(fn, index=True, mode='a', header = not os.path.exists(fn))  
