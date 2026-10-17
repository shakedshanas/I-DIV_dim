# %%
import pandas as pd
import numpy as np
import argparse 
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import euclidean_distances
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as st
from itertools import combinations
import os
import warnings
import pickle
import time
warnings.filterwarnings('ignore')
# %%

parser = argparse.ArgumentParser(description='process every pair of tissues to find genes with low variability')
parser.add_argument('--pair', type=int, help='pair number out of tissue_pairs: 1080 pairs of subtissues in (SMTSD) GTEx v8, 300 pairs of tissues in (SMTS) GTEx v8')
parser.add_argument('--basedon', type=str, choices= ['SMTSD', 'SMTS'], help='pairwise on the subtissue level (SMTSD) or on the organ level (SMTS)')
args = parser.parse_args()
start_time = time.perf_counter()
idx_tissuepair = args.pair + 1
basedon = args.basedon

### PATHS ###

input_path = os.getcwd() + f'/{basedon}_limma/'

#############
# %%

filenames_dict = {}
filenames_dict['div'] = [f for f in os.listdir(basedon + '_limma') if (f.startswith(str(idx_tissuepair) + '.')) & (f.endswith('.div.limma_DE_equiv.csv.gz'))][0]
filenames_dict['tpm'] = [f for f in os.listdir(basedon + '_limma') if (f.startswith(str(idx_tissuepair) + '.')) & (f.endswith('.tpm.limma_DE_equiv.csv.gz'))][0]

div_limma = pd.read_csv(input_path + filenames_dict['div'], compression='gzip', index_col='gene')
tpm_limma = pd.read_csv(input_path + filenames_dict['tpm'], compression='gzip', index_col='gene')

t1 = filenames_dict['div'].split('.')[1]
t2 = filenames_dict['div'].split('.')[2]
pair_list = [t1, t2]
ppair = (t1, t2)
metadata = pd.read_csv('../../GTEx_upstream/metadata.csv.gz', index_col='SAMPID', usecols=['SAMPID', basedon], compression='gzip')
metadata = metadata[metadata[basedon].isin(list(pair_list))]
metadata = metadata.sort_values(basedon)
metadata = metadata[basedon]
paired_tissue_samples = metadata.index.to_list()

gtex_tpm = pd.read_csv('../../GTEx_upstream/gene_tpm.csv.gz', compression='gzip', index_col='accession_number', usecols= ['accession_number'] + paired_tissue_samples)
gtex_tpm = np.log2(gtex_tpm + 1)
gtex_div = pd.read_csv('../../GTEx_upstream/div1.csv.gz', compression='gzip', index_col='accession_number', usecols= ['accession_number'] + paired_tissue_samples)
    
genes_lists = {}
genes_lists['all'] = gtex_tpm.index.to_list()
genes_lists['tpm_eq'] = tpm_limma.loc[(tpm_limma.AveExpr > 0) & (tpm_limma['classification'] == 'equivalent')].index.to_list()
genes_lists['div_eq'] = div_limma.loc[(div_limma.AveExpr > 0) & (div_limma['classification'] == 'equivalent')].index.to_list()

# sample_labels = metadata.map({ppair[lable]: lable for lable in [0,1]}).values
d = {}
d['pair_index'] = idx_tissuepair
d['t1'] = t1
d['t2'] = t2

# %% 
def make_pca(metric, gene_list, sample_list):
    n_components = 2
    pca = PCA(n_components=n_components)
    scaler = StandardScaler()
    
    scaled_features = scaler.fit_transform(metric.loc[gene_list, sample_list].T)
    pca_result = pd.DataFrame(pca.fit_transform(scaled_features), index = sample_list).join(metadata)
    var1 = pca.explained_variance_ratio_[0]
    var2 = pca.explained_variance_ratio_[1]

    label_dict = {l: i for i,l in enumerate(pca_result[basedon].unique())}
    silhouette_sc = silhouette_score(pca_result[[0,1]].values, pca_result[basedon].map(label_dict))


    # rank genes that drive separation
    loadings = pd.DataFrame(
        pca.components_.T, 
        index=gene_list, 
        columns=[f'PC{i+1}' for i in range(pca.n_components_)])


    return pca_result, var1, var2, silhouette_sc, loadings


# %%
# fig, axes = plt.subplots(nrows=2, ncols=len(genes_lists.keys()), figsize=(4 * len(genes_lists.keys()), 8))
for i,k in enumerate(genes_lists.keys()):
    pca_result, var1, var2, silhouette_sc, loadings = make_pca(metric=gtex_tpm, gene_list=genes_lists[k], sample_list=metadata.index)

    d['silhouette_' + k +'_tpm'] = silhouette_sc
    d['eclidean_pca_' + k +'_tpm'] = pd.DataFrame(euclidean_distances(pca_result.drop(columns=basedon),pca_result.drop(columns=basedon))).melt().drop_duplicates('value').value.mean()

for i,k in enumerate(genes_lists.keys()):
    pca_result, var1, var2, silhouette_sc, loadings = make_pca(metric=gtex_div, gene_list=genes_lists[k], sample_list=metadata.index)
    d['silhouette_' + k +'_div'] = silhouette_sc
    d['eclidean_pca_' + k +'_div'] = pd.DataFrame(euclidean_distances(pca_result.drop(columns=basedon),pca_result.drop(columns=basedon))).melt().drop_duplicates('value').value.mean()



# %%

fn= input_path + 'overview_cluster_results.csv'
pd.DataFrame(d.items()).set_index(0).T.to_csv(fn, mode='a', header=not os.path.exists(fn))
# %%
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Script executed in {elapsed_time:.4f} seconds")