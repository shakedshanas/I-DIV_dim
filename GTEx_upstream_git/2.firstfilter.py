import pandas as pd
import numpy as np
import gffpandas.gffpandas as gffpd
import os
import argparse
import gzip
import time

parser = argparse.ArgumentParser(description='filter genes')
parser.add_argument('--data', type=str, choices=['gene_tpm','transcript_tpm'], help='data to filter')
args = parser.parse_args()

start_time = time.perf_counter()
raw_data_path = '/'.join(os.getcwd().split('/')[:-1]) + '/' + 'input_files_git/'
outputpath = os.getcwd() + '/2.outputs/'
if not os.path.exists( outputpath ):
    os.mkdir(outputpath  )
if not os.path.exists( os.getcwd() + '/3.outputs/' ):
    os.mkdir(os.getcwd() + '/3.outputs/'  )



def skip_rows(rows, columnname:str):
    rows[columnname] = rows[columnname].str.split('.', expand=True)[0]
    
    rows_pd_stay = rows[rows[columnname].isin(gene_annotation.index.to_list()  )].index
    rows_pd_skip = [x for x in rows.index if x not in rows_pd_stay]
    rows_pd_skip = [0,1] + [x + 3 for x in rows_pd_skip]    
    return rows_pd_skip

gene_annotation = pd.read_csv('1.outputs/' + 'gene_annotation.csv.gz', compression='gzip', index_col=[0])
metadata = pd.read_csv('1.outputs/' + 'sample_metadata.csv.gz', compression='gzip', index_col=[0])
    
if args.data == 'gene_tpm':
    rows_pd = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', sep='\t', skiprows=2, compression='gzip', usecols= ['Name'])
    rows_pd_skip = skip_rows(rows = rows_pd, columnname = 'Name')
    
    gtex = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', sep='\t', skiprows=rows_pd_skip, usecols= ['Name'] + metadata.index.to_list(), compression='gzip')
    gtex['accession_number'] = gtex.Name.str.split('.', expand=True)[0]
    # gtex = gtex[gtex.Name.isin(gene_annotation.index.to_list())]
    gtex.drop(columns='Name', inplace=True)
    gtex.set_index('accession_number', inplace=True)

if args.data == 'transcript_tpm':
    rows_pd = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz', sep='\t', skiprows=2, compression='gzip', usecols= ['gene_id','transcript_id'])
    rows_pd_skip = skip_rows(rows = rows_pd, columnname = 'gene_id')
    
    gtex = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz', sep='\t', skiprows=rows_pd_skip, usecols=['gene_id','transcript_id'] + metadata.index.to_list())
    gtex = gtex.join(gtex.gene_id.str.split('.', expand=True)[0]).rename(columns={0: 'accession_number'})
    # gtex = gtex[gtex.accession_number.isin(gene_annotation.index.to_list())]
    gtex.set_index('accession_number', inplace=True)
    gtex.drop(columns = 'gene_id', inplace=True)
    
    pd.DataFrame(columns = ['accession_number'] + metadata.index.to_list()).to_csv(os.getcwd() + '/3.outputs/' + 'div1.csv', index=False)
    pd.DataFrame(columns = ['accession_number'] + metadata.index.to_list()).to_csv(os.getcwd() + '/3.outputs/' + 'nisoforms_expressed.csv', index=False)

gtex.to_csv(outputpath  + args.data + '.csv.gz', compression='gzip', index=True)

end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Script executed in {elapsed_time:.4f} seconds")