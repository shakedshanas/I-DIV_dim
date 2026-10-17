# %%
import pandas as pd
import numpy as np
import scipy.stats as st
import sys
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import argparse
# %%

#### PATHS ####
inpt = '/lustre1/home/martinmikl/sshanas/outsourcedata/GTEx_upstream_git/'
###############

def corr_ages_per_tissue(matrix: str):
    
    corr_df_spearman = pd.DataFrame()
    tempdf = gtex_dict[matrix].replace(np.inf,np.nan).replace(-np.inf,np.nan).dropna(axis = 'columns').join(metadata['AGE_int'])
    tempdf.to_csv(baseon + '/' + '.'.join([t, 'spearmanr', matrix, 'raw', 'csv', 'gz']), index=True, compression='gzip')    
    for gene in tqdm(tempdf.drop(columns=['AGE_int']).columns):
        corr_df_spearman.loc[gene, 'coef'], corr_df_spearman.loc[gene, 'p_val'] = st.spearmanr(tempdf[gene], tempdf['AGE_int'])
    
    corr_df_spearman.assign(**{baseon: t}).to_csv(baseon + '/' + '.'.join([t, 'spearmanr', matrix, 'csv', 'gz']), index=True, compression='gzip')

metadata = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.1/sample_metadata.csv.gz', index_col='SAMPID', compression='gzip')
parser = argparse.ArgumentParser(description='run Spearmans correlation on every gene per tissue')
parser.add_argument('--tindex', type=int, choices=range(metadata.SMTSD.nunique()) , help='index of tissue in a list made in this script')
parser.add_argument('--basedon', choices=['SMTS', 'SMTSD'] , help="Should Spearman's test run per SMTS or SMTSD tissue")
args = parser.parse_args()

baseon = args.basedon
tissues_for_process=metadata[baseon].value_counts()[metadata[baseon].value_counts() > 10].index
if args.tindex > len(tissues_for_process):
    sys.exit()
    
t = tissues_for_process[int(args.tindex)]


''' filter samples '''
metadata = metadata[metadata[baseon] == t]
metadata = metadata.assign(age_int = metadata.age_grouped.str[0].astype('int'))
print(int(args.tindex), t)
print('nsamples: ' + str(len(metadata)))

''' load data '''
gtex_dict = {}
gtex_dict = {'div': pd.read_csv(inpt + t + '.IDIV.csv.gz', usecols = ['accession_number'] + list(metadata.index)).drop_duplicates('accession_number').set_index('accession_number'),
             'tpm':  pd.read_csv(inpt + t + '.G.csv.gz', usecols = ['accession_number'] + list(metadata.index) ).drop_duplicates('accession_number').set_index('accession_number')}


for metric in ['tpm', 'div']:
    gtex_dict[metric] = gtex_dict[metric].replace(np.inf,np.nan).replace(-np.inf,np.nan)
    gtex_dict[metric].dropna(axis='index', inplace=True)
    gtex_dict[metric] = gtex_dict[metric].T

common_genes = np.intersect1d(gtex_dict['div'].columns, gtex_dict['tpm'].columns)

for metric in ['tpm', 'div']:
    gtex_dict[metric] = gtex_dict[metric][common_genes]
    print(metric, str(gtex_dict[metric].isna().sum(axis = 'rows').sum()), 'axis=rows')
    print(metric, str(gtex_dict[metric].isna().sum(axis = 'columns').sum()), 'axis=columns')
    '''correlation of TPM/DIV with AGE for every gene in tissue or tissue.'''
    corr_ages_per_tissue(matrix=metric)