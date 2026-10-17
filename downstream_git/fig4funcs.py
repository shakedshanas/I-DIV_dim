import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import scipy.stats as st
from statsmodels.stats.multitest import multipletests

import math
import os
import datetime
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from scipy.stats import gaussian_kde
from sklearn.cluster import KMeans
import gffpandas.gffpandas as gffpd
sns.set_context('paper')
import gseapy
from scipy.spatial.distance import cdist
from sklearn.metrics import silhouette_score
# gseapy.get_library_name()
from sklearn.model_selection import cross_val_predict, cross_val_score
import pprint
from sklearn.metrics import PredictionErrorDisplay
from sklearn.preprocessing import MinMaxScaler
import math
import warnings
from itertools import combinations
warnings.filterwarnings("ignore")
import sys
sys.path.append('/lustre1/home/martinmikl/sshanas/bin')
from div import *
# from plot_functions import *

sns.set_context('paper')
plt.rcParams.update({
    'font.size': 10 * 0.6,
    'axes.titlesize': 10 * 0.6,
    'axes.labelsize': 10 * 0.6,
    'xtick.labelsize': 10 * 0.6,
    'ytick.labelsize': 10 * 0.6,
    'legend.fontsize': 10 * 0.6,
    'figure.titlesize': 12 * 0.6,
    'legend.title_fontsize': 10 * 0.6,
    'legend.fontsize': 10 * 0.6,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans']
})
data_path = '/lustre1/home/martinmikl/sshanas/outsourcedata/GTEx_upstream_git/'
metadata = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/GTEx_upstream_git/metadata.csv.gz', index_col=0)
gene_info = pd.read_csv('/lustre1/home/martinmikl/sshanas/outsourcedata/GTEx_upstream_git/gene_annotation.csv.gz', index_col=0)
id_to_name = gene_info['gene_name'].to_dict()
name_to_id = {id_to_name[k] : k for k in id_to_name.keys()}
age_color_palette = dict(zip(metadata.AGE_grouped.unique(),sns.color_palette('Greens_r', n_colors = metadata.AGE_grouped.nunique() + 1)))
multiple_linear_regression_path = '/lustre1/home/martinmikl/sshanas/outsourcedata/GTEx_pre_downstream_git/fig4_mlr/per_isoform/'

figsoutput_path = 'figs_lme/'
spearman_path = '/lustre1/home/martinmikl/sshanas/outsourcedata/GTEx_pre_downstream_git/spearman/'

def format_small_number(v):
    if 'e' in str(v):
        val, x, powerof = str(v).replace("e", " x 10^").split(' ')
        formatted = ''.join([f"{float(val):.2f}", x, powerof])
    else:
        formatted = f"{v:.2f}"
    return formatted


def annotpval(pv):
    if float(f'{pv:.2e}') == 0.00e+00:
        return 'p<10^-323'
    else:
        return f'p={pv:.2e}'.replace("e", "x10^")

def style_axis(ax):
    ax.tick_params(length=2.5)
    ax.set_yticklabels(ax.get_yticklabels(), x=0)
    ax.set_xticklabels(ax.get_xticklabels(), ha='center', y=0)
    ax.set_xlabel(ax.get_xlabel(), labelpad=0)
    ax.set_ylabel(ax.get_ylabel(), labelpad=0)

def pval_to_asterisks(p):
    if p <= 0.001:
        return '***'
    elif p <= 0.01:
        return '**'
    elif p <= 0.05:
        return '*'
    else:
        return ''

def scale_with_MinMax(df,smpl, minval = 10**-2):
    scaler = MinMaxScaler(feature_range=(minval, 1))
    plot_data = df[[smpl]]
    plot_data = plot_data[plot_data > 0].dropna().sort_values(plot_data.columns[0], ascending=False)

    plot_data = plot_data.assign(minmax = scaler.fit_transform(plot_data[[smpl]]))
    plot_data = plot_data.assign(percentiles = plot_data['minmax'].rank(pct=True, method='first') * 100).sort_values('minmax', ascending=True)
    return plot_data


gene_sets_list = ['GO_Biological_Process_2025', 'GO_Cellular_Component_2025', 'GO_Molecular_Function_2025', 'Reactome_Pathways_2024', 'Human_Phenotype_Ontology', 'KEGG_2021_Human']    
def do_enrichr(gene_list:list, background_genelist:list, trend:str):
    enr = gseapy.enrichr(gene_list=gene_list, # or "./tests/data/gene_list.txt",
                gene_sets= gene_sets_list,
                background=background_genelist,
                outdir=None) # don't write to disk
    df = enr.results
    df['-log10(Adjusted P-value)'] = -np.log10(df['Adjusted P-value'])
    df['-log10(P-value)'] = -np.log10(df['P-value'])
    df['QuerynGenes'] = len(gene_list)
    df['BGnGenes'] = len(background_genelist)
    df['nGenesTerm'] = df.Genes.str.count(';') + 1
    df['trend'] = trend
    return df

def categorise_gene_trend(row, m,fdr_colname, coef_colname, coef_cutoff=0, sign_cutoff=0.05):
    label_dict = {'div': 'I-DIV', 'tpm': 'TPM'}
    if (row[coef_colname] > coef_cutoff) & (row[fdr_colname] <= sign_cutoff):
        return label_dict[m] + '+'

    elif (row[coef_colname] < -coef_cutoff) & (row[fdr_colname] <= sign_cutoff):
        return label_dict[m] + '-'
    else:
        return label_dict[m] + '0'
    
    
def count_gene_significantcoeffs(plotdf, titlestr,filename=''):
    fig, ax = plt.subplots(figsize=(0.75,0.75))

    sns.scatterplot(data = plotdf, x='coef_AGE', y='coef_SEX', alpha=0.2, color='silver', edgecolor='None', s=3)
    sns.scatterplot(data = plotdf.loc[plotdf.fdr_AGE < 0.05], x='coef_AGE', y='coef_SEX', alpha=1, color='Royalblue', edgecolor='None', s=3)
    sns.scatterplot(data = plotdf.loc[plotdf.fdr_SEX < 0.05], x='coef_AGE', y='coef_SEX', alpha=1, color='Gold', edgecolor='None', s=3)
    sns.scatterplot(data = plotdf.loc[(plotdf.fdr_SEX < 0.05) & (plotdf.fdr_AGE < 0.05)], x='coef_AGE', y='coef_SEX', alpha=1, color='limegreen', edgecolor='None', s=5)
    ax.set(title=titlestr, xlabel = 'AGE coefficients', ylabel = 'SEX coefficients')
    
    ax.tick_params(length=2.5)
    ax.set_yticklabels(ax.get_yticklabels(), x=0)

    ax.set_xticklabels(ax.get_xticklabels(), ha='right', rotation = 25, y=0)
    ax.set_xlabel(ax.get_xlabel(), labelpad=0)
    ax.set_ylabel(ax.get_ylabel(), labelpad=0)

    ngenes = {}
    ngenes['AGE only'] = len(plotdf.loc[(plotdf.fdr_AGE < 0.05) & (plotdf.fdr_SEX >= 0.05)])
    ngenes['SEX only'] = len(plotdf.loc[(plotdf.fdr_SEX < 0.05) & (plotdf.fdr_AGE >= 0.05)])
    ngenes['Both'] = len(plotdf.loc[(plotdf.fdr_SEX < 0.05) & (plotdf.fdr_AGE < 0.05)])
    ngenes['Neither'] = len(plotdf.loc[(plotdf.fdr_SEX >= 0.05) & (plotdf.fdr_AGE >= 0.05)])

    s = f'''
    # significant genes AGE: {ngenes['AGE only']}
    # significant genes SEX: {ngenes['SEX only']}
    # significant genes Both: {ngenes['Both']}
    # significant genes Neither: {ngenes['Neither']}
    '''

    print(s)
    
    
    
    if len(filename) > 1:
        
        if os.path.isdir('compare_coeffs') == False:
            os.mkdir('compare_coeffs')
        
        with open(figsoutput_path + 'compare_coeffs/' + filename + '.txt', 'w') as f:
            f.write(s)
            f.close()
        fig.savefig(figsoutput_path + 'compare_coeffs/' + filename + '.png', dpi=300, bbox_inches='tight')
        
def correct_pvalues_with_nan(pvals, method='fdr_bh', alpha=0.05):
    pvals = np.asarray(pvals, dtype=float)
    mask = ~np.isnan(pvals)

    corrected = np.full_like(pvals, np.nan)
    reject = np.full(pvals.shape, False)

    if mask.sum() > 0:
        rej, pvals_corr, _, _ = multipletests(
            pvals[mask], alpha=alpha, method=method
        )
        corrected[mask] = pvals_corr
        reject[mask] = rej

    return reject, corrected
        
    
def make_spearmanres_dfs(m):
    coefficients = gene_info[[]]  
    spearmandfrs = gene_info[[]]  
    for tissue_name in list(set([f.split('.')[0] for f in os.listdir(spearman_path + 'SMTSD/') if (f.startswith('Brain - ')) & (f.endswith('.spearmanr.div.csv.gz'))   ])):
        coefficients = coefficients.join(pd.read_csv(spearman_path + f'SMTSD/{tissue_name}.spearmanr.{m}.csv.gz', 
                                                    usecols = ['Unnamed: 0','coef'], index_col = 'Unnamed: 0').rename(columns = {'coef': tissue_name})  )
        spearmandfrs = spearmandfrs.join(pd.read_csv(spearman_path + f'SMTSD/{tissue_name}.spearmanr.{m}.csv.gz', 
                                                    usecols = ['Unnamed: 0','p_val'], index_col = 'Unnamed: 0').rename(columns = {'p_val': tissue_name})  )

    coefficients = coefficients.join(gene_info.gene_name).set_index('gene_name')
    spearmandfrs = spearmandfrs.join(gene_info.gene_name).set_index('gene_name')
    
    for col in spearmandfrs.columns:
        spearmandfrs[col] = correct_pvalues_with_nan(spearmandfrs[col], alpha=0.05 ,method = 'fdr_bh')[1]
    return coefficients, spearmandfrs


def make_swarmplot_df(subt, gene):
    if gene[:4] == 'ENSG': # if gene is geneid
        gn = id_to_name[gene]
        gid = gene
    else: # if gene is gene name
        gn = gene
        gid = name_to_id[gn]    
    
    div_gene = pd.read_csv(data_path + 'QN.Brain.IDIV.csv', index_col='accession_numbers', usecols = ['accession_numbers'] + metadata[metadata.SMTSD == subt].index.to_list())
    div_gene = div_gene.loc[gid]
    return metadata.loc[metadata.SMTSD == subt, ['AGE_grouped']].join(div_gene).rename(columns={gid:gn})


def get_isoforms_from_gene(gene):
    ''' This function gets a gene name or gene id and returns the isoforms of this gene with their expression values in the readtpm_t.csv.gz file.
    Args:
        gene: A string which is either a gene name or a gene id.
    Returns:
        A dataframe which contains the isoforms of the gene with their expression values in the readtpm_t.csv.gz file.
    '''
    if gene[:4] != 'ENSG':
        gid = gene_info.loc[gene_info['gene_name'] == gene].index[0]
    else:
        gid = gene
    rows = pd.read_csv(data_path + 'Brain.I.csv.gz', usecols=['accession_number'])
    rows_to_skip = rows[~rows['accession_number'].str.contains(gid)].index + 1
    isof_df = pd.read_csv(data_path + 'Brain.I.csv.gz', skiprows=rows_to_skip, index_col='transcript_id')
    return isof_df
    

def get_isoforms_limma(gene, subt):
    if gene[:4] != 'ENSG':
        gid = gene_info.loc[gene_info['gene_name'] == gene].index[0]
    else:
        gid = gene
    
    isoforms_limma = pd.read_csv(multiple_linear_regression_path + f'/{gid}/{gid}.brain.isoform.limma.MLR.csv.gz').merge(gene_info.gene_name, left_on='gene_id', right_index=True)
    gene_pvals = isoforms_limma.loc[isoforms_limma.SMTSD == subt , ['transcript_id', 'SMTSD', 'age_int_p']].pivot(index='transcript_id', columns='SMTSD', values='age_int_p')
    gene_astrc = gene_pvals.applymap(pval_to_asterisks)
    return gene_astrc

    
def isoform_expression_barplots(gene,isoforms_df, subt, savefig=False, path_and_name = ''):
    melted_df = isoforms_df.reset_index()[['transcript_id'] + metadata.loc[metadata.SMTSD == subt].index.to_list()].melt(id_vars=['transcript_id'])
    melted_df['tpm_log2'] = melted_df['value'].apply(lambda x: math.log2(x+1))
    melted_df = melted_df.merge(metadata.age_grouped, left_on='variable', right_index=True)
    oorder = melted_df.groupby('transcript_id')['value'].mean().sort_values(ascending=False).index

    fig, ax = plt.subplots(figsize=(len(isoforms_df)/4,0.75))
    sns.barplot(data = melted_df,  x = 'transcript_id', y ='tpm_log2', errorbar='se', 
                hue='age_grouped', hue_order = sorted(melted_df.age_grouped.unique()), palette=age_color_palette,
                order = oorder,
                err_kws={"color": "dimgrey", "linewidth": 0.75}) # errorbars of standard deviation
    ax.legend(loc='center left', bbox_to_anchor=(1,1), frameon=False)
    
    xticks = [x.get_text().split('.')[0] for x in ax.get_xticklabels()]
    
    ax.tick_params(length=2.5, bottom=False)
    ax.set_yticklabels(ax.get_yticklabels(), x=0)
    ax.set_xticklabels(xticks, rotation=35, ha='right', y=0)
    ax.set_xlabel(ax.get_xlabel(), labelpad=0)
    ax.set_ylabel(ax.get_ylabel(), labelpad=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)     
    ax.patch.set_facecolor('None')
    
    ax.set(xlabel='', ylabel='Isoform level\nTPM [log2]')
    
    
    gene_astrc = get_isoforms_limma(gene=gene, subt=subt)
    
    
    # add asteriks 
    asterisk_y = 4.7
    line_y = 4.9
    for idx, transcript in enumerate(oorder):
        significance = gene_astrc.loc[transcript, subt]
        if significance:
            ax.text(x=idx,y=asterisk_y,s=significance,ha='center',va='bottom')
            ax.hlines(y=line_y,xmin=idx - 0.3,xmax=idx + 0.3,color='black',linewidth=0.7)
        else:
            ax.text(x=idx,y=asterisk_y + 0.12,s='ns',ha='center',va='bottom')
            ax.hlines(y=line_y ,xmin=idx - 0.3,xmax=idx + 0.3,color='black',linewidth=0.7)  
    
    
    
    if (savefig==True) and (len(path_and_name) > 1):
        fig.savefig(path_and_name + '.png', format='png', dpi=300, bbox_inches='tight')
        plt.close()
    
    
def isoform_meanexpression_lineplot(isoforms_df,subt, savefig=False, path_and_name = ''):
    ''' mean isoform expression per age group '''

    plot_data = isoforms_df[[]]
    plot_data = plot_data.assign(**{ag: isoforms_df[metadata.loc[(metadata.SMTSD == subt) & (metadata.AGE_grouped==ag)].index.to_list()].mean(axis='columns') \
        for ag in sorted(metadata.AGE_grouped.unique()) })

    plot_data_catdict = {}
    for ag in sorted(metadata.AGE_grouped.unique()):
        plot_data_catdict[ag] = scale_with_MinMax(df = plot_data, smpl=ag)
        plot_data_catdict[ag] = plot_data_catdict[ag].assign(AGE_grouped = ag)
    plot_data = pd.concat([plot_data_catdict[k][['minmax', 'percentiles','AGE_grouped']]for k in plot_data_catdict.keys()])
    
    fig, ax = plt.subplots(1, 1, figsize=(0.75, 0.75))
    sns.lineplot(data = plot_data, x='percentiles', y='minmax', hue = 'AGE_grouped', palette=age_color_palette, 
                marker='h', markersize=3, linewidth = 0.25, alpha = 1, markeredgewidth=0)
    ax.set_yscale('log')
    ax.legend(loc='center left', bbox_to_anchor=(1,1), frameon=False)
    ax.set(ylabel='Isoform level TPM\nnorm.', xlabel='Isoform level TPM\npercentile')
    style_axis(ax)
    sns.despine()

    if savefig:
        fig.savefig(path_and_name + '.png', format='png', dpi=300, bbox_inches='tight')
        plt.close()