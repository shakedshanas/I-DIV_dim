import pandas as pd
import numpy as np
import sys
import os
sys.path.append('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025')
from div import *
import argparse

parser = argparse.ArgumentParser(description='filter genes')
parser.add_argument('--gf', type=int, choices=range(6))
args = parser.parse_args()

outputpath = os.getcwd() + '/outputdir.5/'

# if not os.path.exists( outputpath ):
#     os.mkdir(outputpath  )

gene_annot = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/gene_metadata.csv.gz', index_col='accession_number')
metadata = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.1/sample_metadata.csv.gz', compression='gzip', index_col=0)

dir_names = {0: 'no_filter',
             1: 'RP',  # RP
             2: 'other_only',   # RP and HK excluded
             3: 'HK-RP',    # only HK
             4: 'other_and_HK-RP', # RP excluded
             5: 'RP_U_HK'}

outputpath += dir_names[args.gf] + '/'
if not os.path.exists( outputpath ):
    os.mkdir(outputpath  )

if args.gf == 0:
    gene_ids = gene_annot.index
elif args.gf == 1: 
    gene_ids = gene_annot[(gene_annot.RP_only ==True) | (gene_annot.HK_and_RP == True)].index
elif args.gf == 2:
    gene_ids = gene_annot[(gene_annot.RP_only ==False) & (gene_annot.HK_and_RP == False) & (gene_annot.HK_only == False)].index
elif args.gf == 3:
    gene_ids = gene_annot[gene_annot.HK_only == True].index
elif args.gf == 4:
    gene_ids = gene_annot[(gene_annot.RP_only == False) & (gene_annot.HK_and_RP == False)].index
elif args.gf == 5:
    gene_ids = gene_annot[(gene_annot.HK_only == True) | (gene_annot.RP_only == True) | (gene_annot.HK_and_RP == True)].index
del gene_annot
print('number of genes in group from GENCODE: ' +  str(len(gene_ids)))
with open (outputpath + 'genes_gencode.txt', 'w') as f:
    for g in gene_ids:
        f.write(g + '\n')

# filter genes
tpm_i = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/readtpm_g.csv.gz', usecols=['accession_number'])
tpm_i = tpm_i[~tpm_i['accession_number'].isin(gene_ids)].index.tolist() # index list to skip
tpm_i = [i + 1 for i in tpm_i] # add +1 to every index, in order to not skip the header
tpm = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/readtpm_g.csv.gz', index_col='accession_number', skiprows=tpm_i)


metadata = metadata.join(pd.Series(tpm.gt(0).sum(axis='rows'), name = 'ngenes_expressed'))
metadata['ngenes_expressed'] = tpm.gt(0).sum(axis='rows')
prct_dict = {1: int(np.percentile(metadata.ngenes_expressed, 1)),
             5: int(np.percentile(metadata.ngenes_expressed, 5)), 
             10: int(np.percentile(metadata.ngenes_expressed, 10))}

for n in prct_dict.keys():
    k = 'ngenes_' + str(n) + 'percentile_HEG'
    metadata[k + '_div'] = tpm.apply(lambda x: DIV(x.nlargest(prct_dict[n]), spreading_percentage = 1, tpm_or_count = 'count'), axis = 'rows')
    metadata[k + '_gt0'] = tpm.apply(lambda x: x.nlargest(prct_dict[n]).gt(0).sum(), axis = 'rows')
    metadata[k + '_meantpm'] = tpm.apply(lambda x: x.nlargest(prct_dict[n]).mean(), axis = 'rows')

metadata['div'] = tpm.apply(lambda x: DIV(x, spreading_percentage = 1, tpm_or_count = 'count'), axis = 'rows')
metadata['gt0'] = tpm.apply(lambda x: x.gt(0).sum(), axis = 'rows')
metadata['meantpm'] = tpm.apply(lambda x: x.mean(), axis = 'rows')
tpm.to_csv(outputpath + 'tpm.csv.gz', compression='gzip', index=True)

del tpm, tpm_i

div_i = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/div1.csv.gz', usecols=['Unnamed: 0'])
div_i = div_i[~div_i['Unnamed: 0'].isin(gene_ids)].index.tolist() # index list to skip
div_i = [i + 1 for i in div_i] # add +1 to every index, in order to not skip the header
div = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/div1.csv.gz', index_col='Unnamed: 0', skiprows=div_i)

for n in prct_dict.keys():
    k = 'ngenes_' + str(n) + 'percentile_HEG'
    metadata[k + '_meanisofdiv'] = div.apply(lambda x: x.nlargest(prct_dict[n]).mean(), axis = 'rows')
metadata['meanisofdiv'] = div.apply(lambda x: x.mean(), axis = 'rows')
div.to_csv(outputpath + 'div1.csv.gz', compression='gzip', index=True)
del div, div_i

nis_cols = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/nisoforms_expressed.csv.gz', nrows=0).columns
nis_idxcol = nis_cols[0]
nis_i = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/nisoforms_expressed.csv.gz', usecols=[nis_idxcol])
nis_i = nis_i[~nis_i[nis_idxcol].isin(gene_ids)].index.tolist() # index list to skip
nis_i = [i + 1 for i in nis_i] # add +1 to every index, in order to not skip the header

nis = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/nisoforms_expressed.csv.gz', index_col=nis_idxcol, skiprows=nis_i)

for n in prct_dict.keys():
    k = 'ngenes_' + str(n) + 'percentile_HEG'
    metadata[k + '_meannis'] = nis.apply(lambda x: x.nlargest(prct_dict[n]).mean(), axis = 'rows')
metadata['meannis'] = nis.apply(lambda x: x.mean(), axis = 'rows')
nis.to_csv(outputpath + 'nis.csv.gz', compression='gzip', index=True)
del nis, nis_i

metadata = metadata.assign(percentile_ngenes_1 = prct_dict[1], percentile_ngenes_5 = prct_dict[5], percentile_ngenes_10 = prct_dict[10])
metadata.to_csv(outputpath + 'metadata.csv.gz', compression='gzip', index=True)