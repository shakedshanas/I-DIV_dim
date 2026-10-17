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

def is_tpm(path):
    ''' gene expression is correlated as log2(TPM+1) '''
    return os.path.basename(path) == 'gene_tpm.csv.gz'

def keep_positive_d1(path):
    ''' I-DIV and number of expressed isoforms are only correlated in samples where the gene is expressed (value > 0) '''
    return os.path.basename(path) in ['div1.csv.gz', 'nisoforms_expressed.csv.gz']

sex_specific_SMTS = ['Uterus', 'Breast', 'Testis', 'Ovary', 'Vagina', 'Prostate']

geneids = pd.read_csv('../../GTEx_upstream/gene_annotation.csv.gz', index_col='accession_number', usecols=['accession_number'], compression='gzip').index 
metadata = pd.read_csv('../../GTEx_upstream/metadata.csv.gz', index_col='SAMPID', compression='gzip')

parser = argparse.ArgumentParser(description='Pearson and spearman correlations of every gene by 2 of the metrics: gene expression, isoform diversity, number of isoforms')
parser.add_argument('-i', type=int, choices=range(len(geneids)), help="gene index")
parser.add_argument('-d1', type=str, help="path to isoform diversity matrix - including name of file")
parser.add_argument('-d2', type=str, help="path to gene expression matrix or numbef of isoforms matrix - including name of file")
parser.add_argument('-b' , choices=['SMTS', 'SMTSD', 'all', 'all_wo_sext'], help="correlation per SMTS/SMTSD tissue, or across all samples ('all_wo_sext': all samples except the sex-specific tissues)")
parser.add_argument('-o', '--out', type=str, help="output file")
args = parser.parse_args()

gid = geneids[args.i]
rows_pd_d1 = pd.read_csv(args.d1, compression='gzip', usecols= ['accession_number'])
rows_pd_skip_d1 = filterrows(rows_pd = rows_pd_d1)

rows_pd_d2 = pd.read_csv(args.d2, compression='gzip', usecols= ['accession_number'])
rows_pd_skip_d2 = filterrows(rows_pd = rows_pd_d2)




del rows_pd_d1, rows_pd_d2, geneids


if args.b == 'all_wo_sext':
    metadata = metadata[~metadata.SMTS.isin(sex_specific_SMTS)]

pearson_df = pd.DataFrame()
if args.b in ['SMTS', 'SMTSD']:
    for t in metadata[args.b].unique(): 
        samples_list = metadata[metadata[args.b]==t].index
        
        d1 = pd.read_csv(args.d1, compression='gzip', skiprows = rows_pd_skip_d1 ,usecols = ['accession_number'] + list(samples_list), index_col='accession_number')
        d2 = pd.read_csv(args.d2, compression='gzip', skiprows = rows_pd_skip_d2 ,usecols = ['accession_number'] + list(samples_list), index_col='accession_number')
        
        if is_tpm(args.d2):
            d2 = log2_trn(df=d2)
        
        if keep_positive_d1(args.d1):
            gn_df = d1.loc[gid][d1.loc[gid] > 0].to_frame().rename(columns = {gid: 'd1'}).join( d2.loc[gid] ).rename(columns = {gid: 'd2'})
        else:
            gn_df = d1.loc[gid].to_frame().rename(columns = {gid: 'd1'}).join( d2.loc[gid] ).rename(columns = {gid: 'd2'})
        
        if len(gn_df) >= 10:
            pearson_df.loc[gid, t + '_coeff'], pearson_df.loc[gid, t + '_pval'] = st.pearsonr(gn_df['d1'], gn_df['d2'])
        else:
            pearson_df.loc[gid, t + '_coeff'], pearson_df.loc[gid, t + '_pval'] = 0, 1


else: 
    samples_list = list(metadata.index)
    d1 = pd.read_csv(args.d1, compression='gzip', skiprows = rows_pd_skip_d1, usecols = ['accession_number'] + samples_list, index_col='accession_number')
    d2 = pd.read_csv(args.d2, compression='gzip', skiprows = rows_pd_skip_d2, usecols = ['accession_number'] + samples_list, index_col='accession_number')

    if is_tpm(args.d2):
        d2 = log2_trn(df=d2)

    if keep_positive_d1(args.d1):
        gn_df = d1.loc[gid][d1.loc[gid] > 0].to_frame().rename(columns = {gid: 'd1'}).join( d2.loc[gid] ).rename(columns = {gid: 'd2'})
    else:
        gn_df = d1.loc[gid].to_frame().rename(columns = {gid: 'd1'}).join( d2.loc[gid] ).rename(columns = {gid: 'd2'})
    
    if len(gn_df) >= 10:
        pearson_df.loc[gid, 'coeff'], pearson_df.loc[gid, 'pval'] = st.pearsonr(gn_df['d1'], gn_df['d2'])
    else:
        pearson_df.loc[gid, 'coeff'], pearson_df.loc[gid, 'pval'] = 0, 1

if os.path.isdir(args.out) == False:
    os.mkdir(args.out)

pearson_df_nf = args.out + '/' + args.b + '.pearson.csv'
pearson_df.to_csv(pearson_df_nf, index=True, mode='a', header = not os.path.exists(pearson_df_nf))  