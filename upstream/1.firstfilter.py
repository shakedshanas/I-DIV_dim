import pandas as pd
# import gffpandas.gffpandas as gffpd
import os
import argparse
import gzip
from time import sleep

parser = argparse.ArgumentParser(description='filter genes')
parser.add_argument('--data', type=int, choices=[0,1,2,3], help='types of data to choose')
args = parser.parse_args()

raw_data_path = '/'.join(os.getcwd().split('/')[:-1]) + '/'
outputpath = os.getcwd() + '/outputdir.1/'
if not os.path.exists( outputpath ):
    os.mkdir(outputpath  )


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

def create_metadata_smpl_genes():
    ''' Step 1 - load sample attributes: technical infromation about the samples '''
    annot          = pd.read_csv(raw_data_path + 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.csv', sep='\t') # explains the samples 
    subject_annot  = pd.read_csv(raw_data_path + 'GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.csv', sep='\t') # explains the sex & age of the donors in each sample

    subject_annot = subject_annot.join(subject_annot.SUBJID.str.split('-', expand=True)) # extract source and donor from subjectID
    # subject_annot.rename(columns={'SEX': 'sex', 'AGE': 'age', 'DTHHRDY': 'death'}, inplace=True)

    annot = annot[['SAMPID', 'SMTS', 'SMTSD', 'SMNABTCHT', 'SMTSISCH', 'SMRIN', 'SMGEBTCH']].join(annot.SAMPID.str.split('-', expand=True)).drop(columns=[2,3,4,5,6])
    metadata = subject_annot.merge(annot, right_on=[0, 1], left_on=[0, 1], how='left').set_index('SAMPID').drop(columns='SUBJID')
    metadata.rename(columns={0: 'source', 1: 'donor'}, inplace=True) # ,'SMTS': 'tissue', 'SMTSD': 'subtissue', 'SMNABTCHT': 'isolation', 'SMTSISCH': 'ischemic', 'SMRIN': 'smpl_quality', 'SMGEBTCH': 'SMGEBTCH'}, inplace=True)

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
    metadata = metadata[metadata.SMRIN >= 6] # sample quality above 6
    smpl_count = metadata.SMTSD.value_counts()
    many_samples_SMTSDs = smpl_count[smpl_count >= 15].index
    metadata = metadata[metadata.SMTSD.isin(many_samples_SMTSDs)] # use SMTSDs with at least 15 samples
    metadata = metadata.assign(age_decade = metadata.AGE.str[0].astype(int))
    metadata['age_grouped'] = metadata.AGE.map({'70-79': '60-79', '60-69': '60-79', '50-59': '50-59', '40-49': '20-49',  '30-39': '20-49', '20-29': '20-49'})
    del smplq_imput_dict_smplq, smpl_count, many_samples_SMTSDs

    ''' Step 3 - first filter of genes: non-coding, X, Y, mitochondria  '''
    ''' EDIT April 6th 2026: gene metadata and first filtering is done elsewhere as a preparation for this process''' 
    # annotation = gffpd.read_gff3('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/gencode.v26.annotation.gff3')

    # gene_annotation = annotation.filter_feature_of_type(['gene']).attributes_to_columns()[['seq_id', 'gene_id', 'gene_name', 'gene_type', 'start', 'end']]
    # gene_annotation['accession_number'] = gene_annotation.gene_id.str.split('.', expand=True)[0]
    # gene_annotation = gene_annotation[gene_annotation.gene_type=='protein_coding']
    # gene_annotation = gene_annotation[~gene_annotation.seq_id.isin(['chrX', 'chrY', 'chrM'])]
    # gene_annotation['length'] = gene_annotation.end.astype(int) - gene_annotation.start.astype(int)
    # # gene_annotation.set_index('accession_number', inplace=True)


    # # gene_annotation = gene_annotation[gene_annotation.transcript_type.isin(['protein_coding',None])]

    # exon_annotation = annotation.filter_feature_of_type(['exon']).attributes_to_columns()[['seq_id', 'type', 'transcript_type', 'gene_name', 'gene_type', 'start', 'end']]
    # exon_annotation = exon_annotation[exon_annotation.gene_name.isin(gene_annotation.gene_name.values)] # use only protein coding genes in the exon df
    # exon_annotation['gene_name_start_end'] = exon_annotation.apply(lambda x: '_'.join([x.gene_name, str(x.start), str(x.end)]), axis='columns')
    # exon_annotation.drop_duplicates('gene_name_start_end', inplace=True)
    # exon_annotation = exon_annotation.merge(pd.Series(exon_annotation.gene_name.value_counts(), name = 'nexons'), right_index=True, left_on='gene_name').drop_duplicates('gene_name')
    # gene_annotation = gene_annotation.merge(exon_annotation[['gene_name', 'nexons']], on = 'gene_name')
    # gene_annotation.set_index('accession_number', inplace=True)
    # del exon_annotation, annotation
    gene_annotation = pd.read_csv('/lustre1/home/martinmikl/sshanas/GENCODE/gene_names_to_id.26.csv.gz', compression='gzip', index_col='accession_number')

    ''' fit the GTEx columns and the metadata samples to each other'''

    gtex_columns = {}
    gtex_columns['gene_reads'] = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz', sep='\t', skiprows=2, nrows=0, compression='gzip').drop(columns=['Name', 'Description']).columns
    gtex_columns['gene_tpm'] = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', sep='\t', skiprows=2, nrows=0, compression='gzip').drop(columns=['Name', 'Description']).columns
    gtex_columns['transcript_reads'] = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz', sep='\t', nrows=0, skiprows=2, compression='gzip').drop(columns=['transcript_id', 'gene_id']).columns
    gtex_columns['transcript_tpm'] = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz', sep='\t', nrows=0, skiprows=2, compression='gzip').drop(columns=['transcript_id', 'gene_id']).columns
    gtex_columns['metadta'] = metadata.index.to_list()

    common_columns = set.intersection(*map(set, gtex_columns.values()))
    common_columns = list(common_columns)
    
    metadata = metadata[metadata.index.isin(common_columns)]
    metadata['system'] = metadata.SMTSD.map(tissue_systems)

    gtex_index = {}
    gtex_index['gene_reads'] = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz', sep='\t', skiprows=2, usecols=['Name']).Name.str.split('.', expand=True)[0].values
    gtex_index['gene_tpm'] = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', sep='\t', skiprows=2, usecols=['Name']).Name.str.split('.', expand=True)[0].values
    gtex_index['transcript_reads'] = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz', sep='\t', skiprows=2, usecols=['gene_id']).gene_id.str.split('.', expand=True)[0].values
    gtex_index['transcript_tpm'] = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz', sep='\t', skiprows=2, usecols=['gene_id']).gene_id.str.split('.', expand=True)[0].values
    gtex_index['metadta'] = gene_annotation.index.to_list()
   
    common_index = set.intersection(*map(set, gtex_index.values())) # create a single set from all the genes across the dictionary 
    common_index = list(common_index) # make it as a list
    
    gene_annotation = gene_annotation[gene_annotation.index.isin(common_index)].join(pd.Series(gtex_index['transcript_reads'].value_counts() , name = 'npt_isoforms') ) # number of potential isoforms
    metadata.to_csv(outputpath + 'sample_metadata.csv.gz', compression='gzip', index=True)
    gene_annotation.to_csv(outputpath + 'gene_metadata.csv.gz', compression='gzip', index=True)
    
    return metadata, gene_annotation

def skip_rows(rows, columnname:str):
    rows[columnname] = rows[columnname].str.split('.', expand=True)[0]
    
    rows_pd_stay = rows[rows[columnname].isin(gene_annotation.index.to_list()  )].index
    rows_pd_skip = [x for x in rows.index if x not in rows_pd_stay]
    rows_pd_skip = [0,1] + [x + 3 for x in rows_pd_skip]    
    return rows_pd_skip

if os.path.isfile(outputpath + 'sample_metadata.csv.gz') and os.path.isfile(outputpath + 'gene_metadata.csv.gz'):
    gene_annotation = pd.read_csv(outputpath + 'gene_metadata.csv.gz', compression='gzip', index_col=[0])
    metadata = pd.read_csv(outputpath + 'sample_metadata.csv.gz', compression='gzip', index_col=[0])
else:
    metadata, gene_annotation = create_metadata_smpl_genes()
    
d = {0:'readcount_g', 
     1: 'readtpm_g',
     2: 'readcount_t', 
     3: 'readtpm_t'}

if args.data == 0:    
    rows_pd = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz', sep='\t', skiprows=2, usecols= ['Name'], compression='gzip')
    rows_pd_skip = skip_rows(rows = rows_pd, columnname = 'Name')
    
    gtex = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz', sep='\t', skiprows=rows_pd_skip, usecols= ['Name'] + metadata.index.to_list(), compression='gzip')
    gtex['accession_number'] = gtex.Name.str.split('.', expand=True)[0]
    # gtex = gtex[gtex.Name.isin(gene_annotation.index.to_list())]
    gtex.drop(columns='Name', inplace=True)
    gtex.set_index('accession_number', inplace=True)
    
if args.data == 1:
    rows_pd = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', sep='\t', skiprows=2, compression='gzip', usecols= ['Name'])
    rows_pd_skip = skip_rows(rows = rows_pd, columnname = 'Name')
    
    # gtex = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', sep='\t', skiprows=2, usecols= ['Name'] + metadata.index.to_list(), compression='gzip')
    gtex = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', sep='\t', skiprows=rows_pd_skip, usecols= ['Name'] + metadata.index.to_list(), compression='gzip')
    gtex['accession_number'] = gtex.Name.str.split('.', expand=True)[0]
    # gtex = gtex[gtex.Name.isin(gene_annotation.index.to_list())]
    gtex.drop(columns='Name', inplace=True)
    gtex.set_index('accession_number', inplace=True)

if args.data == 2:
    rows_pd = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz', sep='\t', skiprows=2, compression='gzip', usecols= ['gene_id'])
    rows_pd_skip = skip_rows(rows = rows_pd, columnname = 'gene_id')
    
    # gtex = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz', sep='\t', skiprows=2, usecols= ['gene_id'] + metadata.index.to_list(), compression='gzip')
    gtex = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz', sep='\t', skiprows=rows_pd_skip, usecols= ['gene_id'] + metadata.index.to_list(), compression='gzip')
    gtex = gtex.join(gtex.gene_id.str.split('.', expand=True)[0]).rename(columns={0: 'accession_number'})
    # gtex = gtex[gtex.accession_number.isin(gene_annotation.index.to_list())]
    gtex.set_index('accession_number', inplace=True)
    gtex.drop(columns = 'gene_id', inplace=True)

if args.data == 3:
    rows_pd = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz', sep='\t', skiprows=2, compression='gzip', usecols= ['gene_id'])
    rows_pd_skip = skip_rows(rows = rows_pd, columnname = 'gene_id')
    
    gtex = pd.read_csv(raw_data_path + 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz', sep='\t', skiprows=rows_pd_skip, usecols=['gene_id'] + metadata.index.to_list())
    gtex = gtex.join(gtex.gene_id.str.split('.', expand=True)[0]).rename(columns={0: 'accession_number'})
    # gtex = gtex[gtex.accession_number.isin(gene_annotation.index.to_list())]
    gtex.set_index('accession_number', inplace=True)
    gtex.drop(columns = 'gene_id', inplace=True)



gtex.to_csv(outputpath  + d[args.data] + '.csv.gz', compression='gzip', index=True)
