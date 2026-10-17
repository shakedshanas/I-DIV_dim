# %%
import pandas as pd
import os
import glob
from functools import reduce
import gseapy
import gffpandas.gffpandas as gffpd
import time
# %%
start_time = time.perf_counter()
input_files_path = '/'.join(os.getcwd().split('/')[:-1]) + '/' + 'input_files_git/'
gene_info = gffpd.read_gff3(input_files_path + 'gencode.v28.annotation.gff3').filter_feature_of_type(['gene']).attributes_to_columns()[['seq_id', 'gene_id', 'gene_name', 'gene_type']]
gene_info = gene_info[(gene_info.gene_type=='protein_coding') & (~gene_info.seq_id.isin(['chrX', 'chrY', 'chrM']))]

''' RP '''
go_bp = gseapy.parser.download_library('GO_Cellular_Component_2025', organism='Human')
# Filter for ribosomal-related gene sets
ribosomal_sets = {k: v for k, v in go_bp.items() if 'ribosom' in k.lower()}
ribosomal_genes = list(set([gene for genes_list in ribosomal_sets.values() for gene in genes_list]))
gene_info.index = gene_info.gene_id.str.split('.',expand=True)[0]
gene_info = gene_info.rename_axis("accession_number", axis="columns")
# gene_info['RP'] = gene_info.gene_name.isin(ribosomal_genes)

bm = gseapy.Biomart()
RP_results = bm.query(dataset='hsapiens_gene_ensembl',
                attributes=['ensembl_gene_id', 'external_gene_name', 'external_synonym'], # , 'go_id'
                filters={'external_gene_name': ribosomal_genes})

RP_ids = list(RP_results.ensembl_gene_id.unique())
gene_info['RP'] = gene_info.index.isin(RP_ids)



ribosomal_only_df = gene_info.loc[(gene_info.RP == True) & (gene_info.gene_name.str[:3] == 'RPL') | (gene_info.gene_name.str[:3] == 'RPS')]

gene_info.to_csv(os.getcwd() + '/' +'gene_annotation.csv.gz', compression='gzip', index=True)




isoform_files = glob.glob(input_files_path + '*.isoforms.results.gz')

dfs = []
for s in isoform_files:
    smplname = s[s.find('day'):s.find('day') + len('day0_0')]
    df = pd.read_csv(s, sep='\t', usecols=['transcript_id', 'gene_id', 'TPM'])
    df = df.rename(columns={'TPM': smplname})
    dfs.append(df)

# Join on both transcript_id and gene_id so gene_id doesn't get duplicated as _x/_y columns
isoform_merged = reduce(lambda left, right: pd.merge(left, right, on=['transcript_id', 'gene_id'], how='outer'),dfs)

sample_cols = sorted(c for c in isoform_merged.columns if c not in ('transcript_id', 'gene_id'))
isoform_merged['accession_number'] = isoform_merged.gene_id.str.split('.', expand = True)[0]
isoform_merged = isoform_merged[['transcript_id', 'accession_number'] + sample_cols]
isoform_merged = isoform_merged.loc[isoform_merged.accession_number.isin(gene_info.index)]
isoform_merged.to_csv(os.getcwd() + '/' + 'transcript_tpm.csv.gz', index=True)
del isoform_merged

gene_files = glob.glob(input_files_path + '*.genes.results.gz')
gene_dfs = []
for s in gene_files:
    smplname = s[s.find('day'):s.find('day') + len('day0_0')]
    df = pd.read_csv(s, sep='\t', usecols=['gene_id', 'TPM'])
    df = df.rename(columns={'TPM': smplname})
    gene_dfs.append(df)

gene_merged = reduce(lambda left, right: pd.merge(left, right, on='gene_id', how='outer'),gene_dfs)
sample_cols = sorted(c for c in gene_merged.columns if c != 'gene_id')
gene_merged['accession_number'] = gene_merged.gene_id.str.split('.', expand = True)[0]
gene_merged = gene_merged[['accession_number'] + sample_cols]
gene_merged = gene_merged.loc[gene_merged.accession_number.isin(gene_info.index)]
gene_merged.to_csv(os.getcwd() + '/' + 'gene_tpm.csv.gz')

end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Script executed in {elapsed_time:.4f} seconds")