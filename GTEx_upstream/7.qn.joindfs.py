#%% 
import pandas as pd
import glob
import argparse
import os
import time
#%% 
parser = argparse.ArgumentParser(description='Join dataframaes into one')
parser.add_argument('--data', type=str, choices=['G', 'I'], help="quantile normalize isoform or gene expression data")
parser.add_argument('--basedon', type=str, choices=['SMTS', 'SMTSDbrain', 'SMTSD'], help="what to join into a single dataframe")
args = parser.parse_args()
start_time = time.perf_counter()
empt_df = pd.read_csv('4.outputs/gene_annotation.csv.gz', usecols=['accession_number'])

def joindataframes_genelevel(df, filelist):

    for f in filelist:
        if args.data=='G':
            df = df.merge(pd.read_csv(f).drop_duplicates(subset='accession_number'), on='accession_number', how='left')
        elif args.data=='IDIV':
            df = df.merge(pd.read_csv(f).drop_duplicates(subset='accession_number'), on='accession_number', how='left')

    return df


def joindataframes_isoformlevel(df, filelist):

    for f in filelist:
        df = df.merge(pd.read_csv(f).drop_duplicates(subset='transcript_id').drop(columns='accession_number'), on='transcript_id', how='left')
    return df


if args.data != 'I':
    if args.basedon=='SMTSDbrain':
        files = glob.glob(os.getcwd() + '/qn_data/' + f'SMTSD/Brain - *.{args.data}.csv.gz')
        df = joindataframes_genelevel(df = empt_df, filelist = files)
        columns_list = [c for c in df.columns if 'Unnamed' not in c]
        df[columns_list].set_index('accession_number').to_csv(os.getcwd() +'/'+ f'QN.Brain.{args.data}.csv.gz', compression='gzip', index=True)
    elif args.basedon=='SMTSD':
        files = glob.glob(os.getcwd() + '/qn_data/' + f'SMTSD/*.{args.data}.csv.gz') 
        df = joindataframes_genelevel(df = empt_df, filelist = files)
        columns_list = [c for c in df.columns if 'Unnamed' not in c]
        df[columns_list].set_index('accession_number').to_csv(os.getcwd() +'/'+ f'QN.SMTSD.{args.data}.csv.gz', compression='gzip', index=True)
    elif args.basedon=='SMTS':
        files = glob.glob(os.getcwd() + '/qn_data/' + f'SMTS/*.{args.data}.csv.gz') 
        df = joindataframes_genelevel(df = empt_df, filelist = files)
        columns_list = [c for c in df.columns if 'Unnamed' not in c]
        df[columns_list].set_index('accession_number').to_csv(os.getcwd() +'/'+ f'QN.SMTS.{args.data}.csv.gz', compression='gzip', index=True)

elif args.data == 'I':
    if args.basedon=='SMTSDbrain':
        files = glob.glob(os.getcwd() + '/qn_data/' + f'SMTSD/Brain - *.I.csv.gz')
        empt_df = pd.read_csv(files[0], usecols=['transcript_id','accession_number']).drop_duplicates(subset='transcript_id')
        df = joindataframes_isoformlevel(df = empt_df, filelist = files)
        df.set_index('transcript_id').to_csv(os.getcwd() +'/'+ f'QN.Brain.I.csv.gz', compression='gzip', index=True)
    elif args.basedon=='SMTSD':
        files = glob.glob(os.getcwd() + '/qn_data/' + f'SMTSD/*.I.csv.gz') 
        empt_df = pd.read_csv(files[0], usecols=['transcript_id','accession_number']).drop_duplicates(subset='transcript_id')
        df = joindataframes_isoformlevel(df = empt_df, filelist = files)
        df.set_index('transcript_id').to_csv(os.getcwd() +'/'+ f'QN.{args.basedon}.I.csv.gz', compression='gzip', index=True)
    elif args.basedon=='SMTS':
        files = glob.glob(os.getcwd() + '/qn_data/' + f'SMTS/*.I.csv.gz') 
        empt_df = pd.read_csv(files[0], usecols=['transcript_id','accession_number']).drop_duplicates(subset='transcript_id')
        df = joindataframes_isoformlevel(df = empt_df, filelist = files)
        df.set_index('transcript_id').to_csv(os.getcwd() +'/'+ f'QN.{args.basedon}.I.csv.gz', compression='gzip', index=True)

end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Script executed in {elapsed_time:.4f} seconds")