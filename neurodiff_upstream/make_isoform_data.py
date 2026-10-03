# %%
import pandas as pd
import os


# m = 'TPM' # length
def isoform_based_data(m:str):
    gene_dict = {file.split('.')[0][-6:]: pd.read_csv(file, compression='gzip', sep='\t', usecols=['transcript_id', 'gene_id', m]) for file in os.listdir() if file.endswith('isoforms.results.gz')}
    i_df = gene_dict['day0_1'].drop(columns=m)
    for k in sorted(gene_dict.keys()):
        i_df = i_df.merge(gene_dict[k], on=['transcript_id', 'gene_id']).rename(columns={m: k})
    i_df.to_csv(m + '.isoform.csv', index=True)

isoform_based_data(m='TPM')
isoform_based_data(m='length')