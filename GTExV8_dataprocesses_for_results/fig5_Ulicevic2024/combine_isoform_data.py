# %%
import pandas as pd
import numpy as np
import os
# %%

isoform_file_name_list = [f for f in os.listdir('.') if f.endswith('isoforms.results.gz')]
genes = pd.read_csv('/Users/shakedshanas/Library/CloudStorage/GoogleDrive-martin_lab@sci.haifa.ac.il/My Drive/General/Datasets/Ulicevic2024/RNAseq/genes.csv', index_col='gene_id', usecols=['gene_id'])
df = pd.read_csv(isoform_file_name_list[0], compression='gzip', usecols=[0,1], sep='\t', index_col=[0])
df = df[df.gene_id.isin(genes.index)]
# %%
for f in isoform_file_name_list:
    samplename = f.split('.')[0].split('_')[-2:] # day_repeat
    samplename = '_'.join(samplename)
    smpldf = pd.read_csv(f, compression='gzip', sep='\t', index_col='transcript_id', usecols=['transcript_id', 'TPM']).rename(columns={'TPM': samplename  })
    df = df.join(smpldf, how='inner')

df.to_csv('isoforms_tpm_combined.csv', compression='gzip', index=True)
# %%
