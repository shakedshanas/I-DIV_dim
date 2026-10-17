# %%
import pandas as pd
import numpy as np
import argparse
import os
import warnings
warnings.filterwarnings("ignore")
import time
# %%

def quantile_normalize(mat):
    """Quantile-normalizes columns (samples) of a 2D numpy array to a common
    reference distribution (the row-wise mean of the sorted columns)."""
    n_rows = len(mat)
    ranks = np.apply_along_axis(lambda col: pd.Series(col).rank(method="average").to_numpy(), axis=0, arr=mat ) # axis=0 means we are ranking each column
    ranks = pd.DataFrame(ranks, columns = mat.columns, index = mat.index)
    ref = np.sort(mat, axis=0).mean(axis=1) # compute the row-wise mean of the sorted values (the "reference distribution")
    xp = np.arange(1, n_rows + 1) # create an array of ranks from 1 to n_rows
    # normalized = np.empty_like(mat, dtype=float)
    normalized = pd.DataFrame(columns=mat.columns, index=mat.index)
    for c in mat.columns:
        normalized.loc[:, c] = np.interp(ranks.loc[:, c], xp, ref)
    return normalized
# %%

parser = argparse.ArgumentParser(description='Quantile normalization of dataframes')
parser.add_argument('--basedon', type=str, choices=['SMTS', 'SMTSD'], help="split the samples to SMTS to SMTSD tissues")
parser.add_argument('--data', type=str, choices=['I', 'G'], help="quantile normalize isoform or gene expression data")
parser.add_argument('--tidx', type=int ,help="tissue index")
args = parser.parse_args()

start_time = time.perf_counter()
gene_info = pd.read_csv(os.getcwd() + '/4.outputs/gene_annotation.csv.gz', index_col='accession_number')
metadata = pd.read_csv(os.getcwd() + '/4.outputs/metadata.csv.gz', compression='gzip', index_col='SAMPID', usecols=['SAMPID',args.basedon])

t = metadata[args.basedon].unique()[args.tidx]
metadata = metadata.loc[metadata[args.basedon] == t]
i_cols = list(metadata.index)
if args.data == 'I':
    df = pd.read_csv(os.getcwd() + '/4.outputs/transcript_tpm.csv.gz', usecols = ['accession_number', 'transcript_id'] + i_cols, index_col='transcript_id')
    transcriptid_to_geneid = df['accession_number']
    quantile_normalize(mat=df.drop(columns='accession_number')).join(transcriptid_to_geneid).to_csv(os.getcwd() + '/qn_data/' + args.basedon + '/' + '.'.join([t,args.data,'csv','gz']), compression='gzip', index = True)
    pd.DataFrame(columns=['accession_number'] + i_cols).to_csv(os.getcwd() + '/qn_data/' + args.basedon + '/' + '.'.join([t,'IDIV','csv']), index=True) # meake headers for calculatio of I-DIV on quantile normalized isoforms
    
elif args.data == 'G':
    df = pd.read_csv(os.getcwd() + '/4.outputs/gene_tpm.csv.gz', usecols = ['accession_number'] + i_cols, index_col='accession_number')
    quantile_normalize(mat=np.log2(df+1)).to_csv(os.getcwd() + '/qn_data/' + args.basedon + '/' + '.'.join([t,args.data,'csv','gz']), compression='gzip', index = True)
    
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Script executed in {elapsed_time:.4f} seconds")