import scipy.stats as st
import pandas as pd
import numpy as np
import os
import argparse
import warnings
warnings.filterwarnings("ignore")

def filterrows(rows_pd):
    # rows_pd_stay = rows_pd[rows_pd['Unnamed: 0'].str.contains(gid)].index
    # rows_pd_skip = [x + 1  for x in rows_pd.index if x not in rows_pd_stay] 
    rows_pd_skip = rows_pd[~rows_pd.iloc[:, 0].str.contains(gid)].index + 1
    return rows_pd_skip

def log2_trn(df):
    return np.log2(df + 1)

geneids = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/gene_metadata.csv.gz', index_col='accession_number', usecols=['accession_number'], compression='gzip').index 
metadata = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.5/no_filter/metadata.csv.gz', index_col='SAMPID', compression='gzip')
metadata = metadata[~metadata.SMTS.isin(['Uterus', 'Breast', 'Testis', 'Ovary', 'Vagina' ,'Prostate'])]
chosen_columns = metadata.index
parser = argparse.ArgumentParser(description='Pearson and spearman correlations of every gene by 2 of the metrics: gene expression, isoform diversity, number of isoforms')
parser.add_argument('-i', type=int, choices=range(len(geneids)), help="gene index")
args = parser.parse_args()

gid = geneids[args.i]


div_col0 = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/div1.csv.gz', compression='gzip', nrows=0).columns[0]
nis_col0 = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/nisoforms_expressed.csv.gz', compression='gzip', nrows=0).columns[0]

rows_pd_div = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/div1.csv.gz', compression='gzip', usecols= [div_col0])
rows_pd_skip_div = filterrows(rows_pd = rows_pd_div)

rows_pd_nis = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/nisoforms_expressed.csv.gz', compression='gzip', usecols= [nis_col0])
rows_pd_skip_nis = filterrows(rows_pd = rows_pd_nis)

del rows_pd_div, rows_pd_nis, geneids, metadata

spearman_df = pd.DataFrame()
pearson_df = pd.DataFrame()

div = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/div1.csv.gz', compression='gzip', skiprows = rows_pd_skip_div, index_col=div_col0, usecols= [div_col0] + list(chosen_columns))
nis = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/nisoforms_expressed.csv.gz', compression='gzip', skiprows = rows_pd_skip_nis, index_col=nis_col0, usecols= [nis_col0] + list(chosen_columns))


gn_df = div.loc[gid][div.loc[gid] > 0].to_frame().rename(columns = {gid: 'div'}).join( nis.loc[gid] ).rename(columns = {gid: 'nis'})

if len(gn_df) >= 10:
    spearman_df.loc[gid,'coeff'], spearman_df.loc[gid, 'pval'] = st.spearmanr(gn_df['div'], gn_df['nis'])
    pearson_df.loc[gid, 'coeff'], pearson_df.loc[gid, 'pval'] = st.pearsonr(gn_df['div'], gn_df['nis'])
else:
    spearman_df.loc[gid,'coeff'], spearman_df.loc[gid, 'pval'] = 0, 1
    pearson_df.loc[gid, 'coeff'], pearson_df.loc[gid, 'pval'] = 0, 1

output_path = '/lustre1/home/martinmikl/sshanas/outsourcedata/downstream_2025_6/fast_1_2/corr_div_tpm/SMTS_wo_sext/'
spearman_df_nf = output_path + 'div_nis.spearman.csv'
pearson_df_nf = output_path + 'div_nis.pearson.csv'
spearman_df.to_csv(spearman_df_nf, index=True, mode='a', header = not os.path.exists(spearman_df_nf))
pearson_df.to_csv(pearson_df_nf, index=True, mode='a', header = not os.path.exists(pearson_df_nf))  