# %%
import pandas as pd
import numpy as np
import gffpandas.gffpandas as gffpd
import gseapy
import os
import argparse
import gzip
import time

''' make gene annotation data '''
start_time = time.perf_counter()
bm = gseapy.Biomart()

raw_data_path = '/'.join(os.getcwd().split('/')[:-1]) + '/' + 'input_files_git/'
gene_annotation = gffpd.read_gff3(raw_data_path + 'gencode.v26.annotation.gff3').filter_feature_of_type(['gene']).attributes_to_columns()[['seq_id', 'gene_id', 'gene_name', 'gene_type', 'start', 'end']]
gene_annotation['accession_number'] = gene_annotation.gene_id.str.split('.', expand=True)[0]
gene_annotation = gene_annotation[(gene_annotation.gene_type=='protein_coding') & (~gene_annotation.seq_id.isin(['chrX', 'chrY', 'chrM']))]
gene_annotation.set_index('accession_number', inplace=True)    

''' count exons and gene length '''
exon_annotation = gffpd.read_gff3(raw_data_path + 'gencode.v26.annotation.gff3').filter_feature_of_type(['exon']).attributes_to_columns()[['seq_id', 'gene_id', 'gene_name', 'gene_type', 'start', 'end']]
exon_annotation = exon_annotation.loc[exon_annotation.accession_number.isin(gene_annotation.index)]
exon_annotation['gene_name_start_end'] = exon_annotation.apply(lambda x: '_'.join(sorted([x['gene_name'], str(x['start']), str(x['end'])])), axis='columns')
exon_annotation.drop_duplicates('gene_name_start_end', inplace=True) # now we get unique exons
exon_annotation['accession_number'] = exon_annotation.gene_id.str.split('.', expand=True)[0]
gene_annotation = gene_annotation.join(pd.Series(exon_annotation.accession_number.value_counts(), name = 'nexons'))
gene_annotation = gene_annotation.assign(length = gene_annotation['end'] - gene_annotation['start'])

''' get ribosomal genes '''

go_bp = gseapy.parser.download_library('GO_Cellular_Component_2025', organism='Human')

# Filter for ribosomal-related gene sets
ribosomal_sets = {k: v for k, v in go_bp.items() if 'ribosom' in k.lower()}
ribosomal_genes = list(set([gene for genes_list in ribosomal_sets.values() for gene in genes_list]))
del ribosomal_sets

RP_results = bm.query(dataset='hsapiens_gene_ensembl',
                attributes=['ensembl_gene_id', 'external_gene_name', 'external_synonym'], # , 'go_id'
                filters={'external_gene_name': ribosomal_genes})

RP_ids = list(RP_results.ensembl_gene_id.unique())

gene_annotation['RP_genes'] = gene_annotation.index.isin(RP_ids)

if not os.path.isdir('1.outputs'):
    os.mkdir('1.outputs')

gene_annotation.to_csv('1.outputs/gene_annotation.csv.gz', index=True, compression='gzip')
# %%



''' make sample metadata '''
tissue_systems = {
    'Skin - Sun Exposed (Lower leg)': 'Integumentary System',
    'Adipose - Subcutaneous': 'Integumentary System',
    'Muscle - Skeletal': 'Musculoskeletal System',
    'Artery - Tibial': 'Cardiovascular System',
    'Nerve - Tibial': 'Nervous System',
    'Heart - Atrial Appendage': 'Cardiovascular System',
    'Heart - Left Ventricle': 'Cardiovascular System',
    'Lung': 'Respiratory System',
    'Skin - Not Sun Exposed (Suprapubic)': 'Integumentary System',
    'Minor Salivary Gland': 'Digestive System',
    'Brain - Cortex': 'Nervous System',
    'Brain - Cerebellum': 'Nervous System',
    'Whole Blood': 'Cardiovascular System',
    'Cells - Cultured fibroblasts': 'Connective Tissue',
    'Thyroid': 'Endocrine System',
    'Artery - Aorta': 'Cardiovascular System',
    'Pancreas': 'Endocrine System',
    'Adipose - Visceral (Omentum)': 'Digestive System',
    'Pituitary': 'Endocrine System',
    'Esophagus - Gastroesophageal Junction': 'Digestive System',
    'Uterus': 'Reproductive System',
    'Breast - Mammary Tissue': 'Reproductive System',
    'Testis': 'Reproductive System',
    'Colon - Sigmoid': 'Digestive System',
    'Brain - Hippocampus': 'Nervous System',
    'Brain - Anterior cingulate cortex (BA24)': 'Nervous System',
    'Brain - Caudate (basal ganglia)': 'Nervous System',
    'Brain - Nucleus accumbens (basal ganglia)': 'Nervous System',
    'Brain - Putamen (basal ganglia)': 'Nervous System',
    'Brain - Frontal Cortex (BA9)': 'Nervous System',
    'Brain - Cerebellar Hemisphere': 'Nervous System',
    'Brain - Hypothalamus': 'Nervous System',
    'Liver': 'Digestive System',
    'Artery - Coronary': 'Cardiovascular System',
    'Esophagus - Mucosa': 'Digestive System',
    'Esophagus - Muscularis': 'Digestive System',
    'Adrenal Gland': 'Endocrine System',
    'Cells - EBV-transformed lymphocytes': 'Immune System',
    'Brain - Substantia nigra': 'Nervous System',
    'Brain - Spinal cord (cervical c-1)': 'Nervous System',
    'Ovary': 'Reproductive System',
    'Vagina': 'Reproductive System',
    'Stomach': 'Digestive System',
    'Colon - Transverse': 'Digestive System',
    'Prostate': 'Reproductive System',
    'Kidney - Cortex': 'Urinary System',
    'Brain - Amygdala': 'Nervous System'
}


''' Step 1 - load sample attributes: technical infromation about the samples '''
annot          = pd.read_csv(raw_data_path + 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t') # explains the samples 
subject_annot  = pd.read_csv(raw_data_path + 'GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t') # explains the sex & age of the donors in each sample

subject_annot = subject_annot.join(subject_annot.SUBJID.str.split('-', expand=True)) # extract source and donor from subjectID
# subject_annot.rename(columns={'SEX': 'sex', 'AGE': 'age', 'DTHHRDY': 'death'}, inplace=True)

annot = annot[['SAMPID', 'SMTS', 'SMTSD', 'SMNABTCHT', 'SMTSISCH', 'SMRIN', 'SMGEBTCH']].join(annot.SAMPID.str.split('-', expand=True)).drop(columns=[2,3,4,5,6])
metadata = subject_annot.merge(annot, right_on=[0, 1], left_on=[0, 1], how='left').set_index('SAMPID').drop(columns='SUBJID')
metadata.rename(columns={0: 'source', 1: 'donor'}, inplace=True) # ,'SMTS': 'tissue', 'SMTSD': 'subtissue', 'SMNABTCHT': 'isolation', 'SMTSISCH': 'ischemic', 'SMRIN': 'smpl_quality', 'SMGEBTCH': 'batch'}, inplace=True)

del annot, subject_annot
''' Step 2 - filter samples 1: death types, RIN, number of samples in SMTSD '''
# categorized only 2 death types
metadata = metadata[metadata.DTHHRDY.isin([1,2])] # death 1 and 2 samples will be from healthy donors.

''' impute the SMTSISCH time and RIN'''
# impute missing smaple quality from the mean of the SMTSD. 

smplq_imput_dict_smplq = metadata.groupby('SMTSD')['SMRIN'].mean().to_dict()
smplq_imput_dict_SMTSISCH = metadata.groupby('SMTSD')['SMTSISCH'].mean().to_dict()
for i in metadata.index:
    if pd.isna(metadata.loc[i, 'SMRIN']) == True:
        metadata.loc[i, 'SMRIN'] = smplq_imput_dict_smplq[metadata.loc[i,'SMTSD']]
    if pd.isna(metadata.loc[i,'SMTSISCH']) == True:
        metadata.loc[i,'SMTSISCH'] = smplq_imput_dict_SMTSISCH[metadata.loc[i,'SMTSD']]
    elif metadata.loc[i, 'SMTSISCH'] < 0:
        metadata.loc[i, 'SMTSISCH'] = 0

# add "singleton" to batches that appear once
batch_counts = metadata['SMGEBTCH'].value_counts()
singleton_batches = batch_counts[batch_counts == 1].index
metadata.loc[metadata['SMGEBTCH'].isin(singleton_batches), 'SMGEBTCH'] = 'singleton'

# filter of RIN is or greater than 6
metadata = metadata[metadata.SMRIN >= 6] 
smpl_count = metadata.SMTSD.value_counts()
many_samples_SMTSDs = smpl_count[smpl_count >= 15].index
metadata = metadata[metadata.SMTSD.isin(many_samples_SMTSDs)] # use SMTSDs with at least 15 samples
metadata['AGE_grouped'] = metadata.AGE.map({'70-79': '60-79', '60-69': '60-79', '50-59': '50-59', '40-49': '20-49',  '30-39': '20-49', '20-29': '20-49'})
metadata = metadata.assign(AGE_int = metadata.AGE_grouped.str[0].astype(int))

del smplq_imput_dict_smplq, smpl_count, many_samples_SMTSDs

''' Step 3 - first filter of genes: non-coding, X, Y, mitochondria  '''
# gene_annotation = pd.read_csv('/lustre1/home/martinmikl/sshanas/GENCODE/gene_names_to_id.26.csv.gz', compression='gzip', index_col='accession_number')

''' fit the GTEx columns and the metadata samples to each other'''

gtex_columns = {}
gtex_columns['gene_tpm'] = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', sep='\t', skiprows=2, nrows=0, compression='gzip').drop(columns=['Name', 'Description']).columns
gtex_columns['transcript_tpm'] = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz', sep='\t', nrows=0, skiprows=2, compression='gzip').drop(columns=['transcript_id', 'gene_id']).columns
gtex_columns['metadta'] = metadata.index.to_list()

common_columns = set.intersection(*map(set, gtex_columns.values()))
common_columns = list(common_columns)

metadata = metadata[metadata.index.isin(common_columns)]
metadata['system'] = metadata.SMTSD.map(tissue_systems)

gtex_index = {}
gtex_index['gene_tpm'] = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', sep='\t', skiprows=2, usecols=['Name']).Name.str.split('.', expand=True)[0].values
gtex_index['transcript_tpm'] = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz', sep='\t', skiprows=2, usecols=['gene_id']).gene_id.str.split('.', expand=True)[0].values
gtex_index['metadta'] = gene_annotation.index.to_list()

common_index = set.intersection(*map(set, gtex_index.values())) # create a single set from all the genes across the dictionary 
common_index = list(common_index) # make it as a list

values, counts = np.unique(gtex_index['transcript_tpm'], return_counts=True)
value_counts = dict(zip(values, counts))

gene_annotation = gene_annotation[gene_annotation.index.isin(common_index)].join(pd.Series(value_counts, name = 'npt_isoforms') ) # number of potential isoforms
metadata.to_csv(os.getcwd() + '/1.outputs/' + 'sample_metadata.csv.gz', compression='gzip', index=True)
gene_annotation.to_csv(os.getcwd() + '/1.outputs/' + 'gene_annotation.csv.gz', compression='gzip', index=True)


end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Script executed in {elapsed_time:.4f} seconds")