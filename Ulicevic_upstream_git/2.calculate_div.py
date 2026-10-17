 # %%
import pandas as pd
import os
import warnings
warnings.filterwarnings("ignore")
import sys
sys.path.append('/lustre1/home/martinmikl/sshanas/bin')
from div import *
import time
# %%

start_time = time.perf_counter()
isoform_df = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/Ulicevic_upstream_git/transcript_tpm.csv.gz').drop(columns = ['Unnamed: 0', 'transcript_id'])
dataframes = {}
dataframes['div1'] = pd.DataFrame(columns = isoform_df.columns)
dataframes['nisoforms_expressed'] = pd.DataFrame(columns = isoform_df.columns)
# Iterate through groups and perform operations
for gene, group in isoform_df.groupby(by='accession_number'):

    dataframes['div1'].loc[gene]=group.set_index('accession_number').apply(lambda x: DIV(arr=x, spreading_percentage=1, tpm_or_count='tpm'), axis = 'rows')
    dataframes['nisoforms_expressed'].loc[gene] = group.set_index('accession_number').apply(lambda x: len(x.values[x.values>0]), axis = 'rows')
     
for k in dataframes.keys():
    dataframes[k].rename_axis("accession_number", axis="columns").to_csv(k + '.csv.gz', index=True, compresion='gzip')
    
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Script executed in {elapsed_time:.4f} seconds")