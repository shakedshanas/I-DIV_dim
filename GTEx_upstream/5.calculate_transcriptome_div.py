# %%
import pandas as pd
import numpy as np
import sys
import os
sys.path.append('..')  # div.py lives in the repository root
from div import *
import scipy.stats as st
import time
# %%
'''
On the sample level, I add tothe metadata: 
1. calculate transcriptome diversity (T-DIV), mean isoform diversity (µ I-DIV) per sample.
2. add tissue abbreviations columns we use in the paper 
3. add supplementary information of cell heterogneity estimation 

I make more directories for the quantile normalization steps
'''
start_time = time.perf_counter()
gene_annot = pd.read_csv(os.getcwd() + '/4.outputs/gene_annotation.csv.gz', index_col='accession_number')
metadata = pd.read_csv(os.getcwd() + '/1.outputs/sample_metadata.csv.gz', compression='gzip', index_col=0)
deconvolution_supplementary_tables_path = '../input_files/Donovan2020/'




''' three gene lists ''' 
all_genes_list = gene_annot.index
RP_excluded = gene_annot.loc[gene_annot.RP_genes == False].index
RP_only = gene_annot.loc[gene_annot.RP_genes == True].index
del gene_annot

tpm = pd.read_csv(os.getcwd() + '/4.outputs/gene_tpm.csv.gz', index_col='accession_number')
metadata['div'] = tpm.loc[RP_excluded].apply(lambda x: DIV(x, spreading_percentage = 1, tpm_or_count = 'count'), axis = 'rows') # T-DIV
metadata['RP_div'] = tpm.loc[RP_only].apply(lambda x: DIV(x, spreading_percentage = 1, tpm_or_count = 'count'), axis = 'rows') # T-DIV
metadata['15834_div'] = tpm.loc[all_genes_list].apply(lambda x: DIV(x, spreading_percentage = 1, tpm_or_count = 'count'), axis = 'rows') # T-DIV
metadata['ngenes_expressed'] = tpm.apply(lambda x: x.gt(0).sum(), axis = 'rows')
del tpm


div = pd.read_csv(os.getcwd() + '/4.outputs/div1.csv.gz', index_col='accession_number')
metadata['meanisofdiv'] = div.apply(lambda x: x.mean(), axis = 'rows')

del div
nis = pd.read_csv(os.getcwd() + '/4.outputs/nisoforms_expressed.csv.gz', index_col='accession_number')
metadata['meannis'] = nis.apply(lambda x: x.mean(), axis = 'rows')
del nis

''' 
in some cases, I will want to use T-DIV, I will want to use samples which express more than 1 percentile of the distribution of number of genes expressed. 
this could be used as a kind of threshold for number of genes expressed, comapring all the other samples I have. 
'''
bottom_1_percentile = int(np.percentile(metadata.ngenes_expressed, 1)) 
metadata['percentile_ngenes_1'] = bottom_1_percentile

''' cellular composition diversity '''
cellular_compositions = {}
smpls_idcs = [5,6,7,8,9,11,12,13,14,15,18]
for i in smpls_idcs:
    
    ''' 
    most of the supplementary data I am using here is estimated by mouse gene expression signatures. 
    but some are estimated by human gene expression signature.
    I had made a decision with myself to use only mouse, for persistency 
    
    
    https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-14561-0/MediaObjects/41467_2020_14561_MOESM28_ESM.pdf
    
    It goes as such: 
    Supplementary Data 4 legend. Cellular composition estimates of bulk RNA-seq from GTEx Artery (Aorta) using mouse aorta signature genes ;; i=5
    Supplementary Data 5 legend. Cellular composition estimates of bulk RNA-seq from GTEx Heart (Atrial appendage) using mouse heart (subsetted to atrium cells) signature genes ;; i=6
    Supplementary Data 6 legend. Cellular composition estimates of bulk RNA-seq from GTEx Brain (Amygdala, Anterior cingulate cortex (BA24), Caudate (basal ganglia), Cerebellar Hemisphere, Cerebellum, Cortex, Frontal Cortex (BA9), Hippocampus, Hypothalamus, 
        Nucleus accumbens (basal ganglia), Putamen (basal ganglia), Spinal cord (cervical c-1), Substantia nigra) using mouse brain non-microglia signature genes ;; i=7
    Supplementary Data 7 legend. Cellular composition estimates of bulk RNA-seq from GTEx Colon (Sigmoid and transverse) using mouse colon signature genes ;; i=8
    Supplementary Data 8 legend. Cellular composition estimates of bulk RNA-seq from GTEx Adipose (Subcutaneous, Visceral (Omentum)) using mouse fat signature genes ;; i=9
    Supplementary Data 9 legend. Cellular composition estimates of bulk RNA-seq from GTEx Liver using human liver signature genes ;; i=10 ## HUMAN
    Supplementary Data 10 legend. Cellular composition estimates of bulk RNA-seq from GTEx Liver using mouse liver signature genes ;; i=11
    Supplementary Data 11 legend. Cellular composition estimates of bulk RNA-seq from GTEx Breast (mammary) using mouse mammary signature genes ;; i=12
    Supplementary Data 12 legend. Cellular composition estimates of bulk RNA-seq from GTEx Muscle (skeletal muscle) using mouse muscle signature genes ;; i=13
    Supplementary Data 13 legend. Cellular composition estimates of bulk RNA-seq from GTEx Pancreas using mouse pancreas signature genes ;; i=14
    Supplementary Data 14 legend. Cellular composition estimates of bulk RNA-seq from GTEx Skin (Not sun exposed (Suprapubic), Sun exposed (Lower leg)) using mouse skin signature genes ;; i=15
    Supplementary Data 15 legend. Cellular composition estimates of bulk RNA-seq from GTEx Skin (Not sun exposed (Suprapubic), Sun exposed (Lower leg)) using human epidermis signature genes ;; i=16 ## HUMAN
    Supplementary Data 16 legend. Cellular composition estimates of bulk RNA-seq from GTEx Spleen using mouse spleen signature genes ;; i=17
    Supplementary Data 17 legend. Cellular composition estimates of bulk RNA-seq from GTEx Heart (Left ventricle) using mouse heart (subsetted to ventricle cells) signature genes ;; i=18    
    '''
    
    cellular_compositions[i] = pd.read_csv(deconvolution_supplementary_tables_path + f'41467_2020_14561_MOESM{i}_ESM.csv', index_col='Input.Sample').join(metadata.SMTSD).assign(n = i).dropna().drop_duplicates()
    cell_columns = cellular_compositions[i].drop(columns = ['P.value', 'SMTSD','n']).columns
    cellular_compositions[i]['Donovan_cell_div'] = cellular_compositions[i].apply(lambda row: -1 * sum([row[col] * np.log(row[col])  for col in cell_columns if row[col] > 0 ])   , axis = 'columns') # shannon entropy calculation

# because there are two heart samples that are focusing on different genes (?) - I calculated the mean diversity to the samples from different cell composition estimates.
s = cellular_compositions[6]['Donovan_cell_div'].to_frame(name = '6').join(pd.Series(cellular_compositions[18]['Donovan_cell_div'], name = '18')).mean(axis='columns')
s.name = 'Donovan_cell_div' 
cellular_compositions[6] = cellular_compositions[6].drop(columns = 'Donovan_cell_div').join(s)
del cellular_compositions[18]
cell_est_div = pd.concat([cellular_compositions[i] for i in cellular_compositions.keys() ])

metadata = metadata.join(pd.concat([cellular_compositions[i]['Donovan_cell_div'] for i in cellular_compositions.keys() ]))

cellular_div_corr_transcriptome_div = pd.DataFrame()
for subt in metadata.dropna().SMTSD.unique():
    temp_df  = metadata.loc[metadata['SMTSD'] == subt, ['div', 'Donovan_cell_div']].dropna()
    cellular_div_corr_transcriptome_div.loc[subt, 'coefficient'], cellular_div_corr_transcriptome_div.loc[subt, 'p-value'] = st.pearsonr(temp_df['div'], temp_df['Donovan_cell_div'])
    cellular_div_corr_transcriptome_div.loc[subt, 'nsamples'] = len(temp_df)
cellular_div_corr_transcriptome_div = cellular_div_corr_transcriptome_div.assign(fdr = st.false_discovery_control(cellular_div_corr_transcriptome_div['p-value'], method = 'bh'))
cellular_div_corr_transcriptome_div = cellular_div_corr_transcriptome_div.assign(mlogpfdr = -np.log10(cellular_div_corr_transcriptome_div['fdr']))

# cellular_div_corr_transcriptome_div.to_csv(source_data_path + '3d.csv.gz', compression='gzip', index=True)
cellular_div_corr_transcriptome_div.to_csv(os.getcwd() + '/4.outputs/CDIV_TDIV_pearosnr.csv.gz', compression='gzip', index=True)

''' add abbreviation columns '''
abb_SMTS = {t: t.replace(' Tissue', '').replace('Gland', 'G.').replace('Vessel', 'V.') for t in metadata.SMTS.unique()}
metadata['SMTS_abb'] = metadata.SMTS.map(abb_SMTS)
brain_regions_abb = {
    'Brain - Cortex': 'CTX',
    'Brain - Cerebellum': 'CBL',
    'Brain - Hippocampus': 'HIP',
    'Brain - Anterior cingulate cortex (BA24)': 'ACC',
    'Brain - Caudate (basal ganglia)': 'CAU',
    'Brain - Nucleus accumbens (basal ganglia)': 'NAC',
    'Brain - Putamen (basal ganglia)': 'PUT',
    'Brain - Frontal Cortex (BA9)': 'FCX',
    'Brain - Cerebellar Hemisphere': 'CBH',
    'Brain - Hypothalamus': 'HYP',
    'Brain - Substantia nigra': 'SNI',
    'Brain - Spinal cord (cervical c-1)': 'SPC',
    'Brain - Amygdala': 'AMY',
}
metadata['SMTSDbrain_abb'] = metadata.SMTSD.map(brain_regions_abb)
metadata.to_csv(os.getcwd() + '/4.outputs/metadata.csv.gz', compression='gzip', index=True)
del metadata
''' directories for the quantile normalization steps '''

if not os.path.isdir(os.getcwd() + '/qn_data/'):
    os.mkdir(os.getcwd() + '/qn_data/')
if not os.path.isdir(os.getcwd() + '/qn_data/SMTS/'):
    os.mkdir(os.getcwd() + '/qn_data/SMTS/')
if not os.path.isdir(os.getcwd() + '/qn_data/SMTSD/'):
    os.mkdir(os.getcwd() + '/qn_data/SMTSD/')
    
    
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Script executed in {elapsed_time:.4f} seconds")
