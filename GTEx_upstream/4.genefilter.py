# %%
import pandas as pd
import numpy as np
import math
import argparse
import warnings
import os
import scipy.stats as st
# %%
if not os.path.exists( os.getcwd() + '/4.outputs/' ):
    os.mkdir(os.getcwd() + '/4.outputs/'  )

''' filter genes to have at least 1% of the samples with a value > 0 in the diversity file. 
This is to reduce the number of genes that we have to load in the next steps, which will speed up the process.'''
div = pd.read_csv(os.getcwd() + '/3.outputs/' + 'div1.csv', index_col='accession_number')
div = div.astype(float)

gene_info = pd.read_csv(os.getcwd() + '/1.outputs/' + 'gene_annotation.csv.gz', index_col='accession_number')
idiv_gt0_percentage = (div.apply(lambda x: len(x[x>0]), axis='columns') / div.apply(lambda x: len(x), axis='columns')) 
gene_info = gene_info.loc[idiv_gt0_percentage[idiv_gt0_percentage > 0.01].index].join(pd.Series (idiv_gt0_percentage,name = 'IDIV_gt0'))
gene_info.to_csv(os.getcwd() + '/4.outputs/'  + 'gene_annotation.csv.gz', index=True)

div.loc[gene_info.index].to_csv(os.getcwd() + '/4.outputs/' + 'div1.csv.gz', compression='gzip', index = True)
del div

nis = pd.read_csv(os.getcwd() + '/3.outputs/' + 'nisoforms_expressed.csv', index_col='accession_number')
nis.loc[gene_info.index].to_csv(os.getcwd() + '/4.outputs/' + 'nisoforms_expressed.csv.gz', compression='gzip', index = True)
del nis

tpm = pd.read_csv(os.getcwd() + '/2.outputs/' + 'gene_tpm.csv.gz', index_col='accession_number')
tpm.loc[gene_info.index].to_csv(os.getcwd() + '/4.outputs/' + 'gene_tpm.csv.gz', compression='gzip', index = True)
tpm.loc[gene_info.index].apply(st.zscore, axis = 'rows').to_csv(os.getcwd() + '/4.outputs/' + 'gene_tpm.zs.csv.gz', index=True, compression='gzip') # every sample: mean 0 st 1
del tpm 

i_tpm = pd.read_csv(os.getcwd() + '/2.outputs/' + 'transcript_tpm.csv.gz', index_col='transcript_id')
i_tpm.loc[i_tpm.accession_number.isin(gene_info.index)].to_csv(os.getcwd() + '/4.outputs/' + 'transcript_tpm.csv.gz', compression='gzip', index = True)
del i_tpm
