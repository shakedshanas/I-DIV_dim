# %%
import pandas as pd
import numpy as np
import gffpandas.gffpandas as gffpd
import gseapy
import re
import scanpy as sc
import argparse
import sys
# %%
''' 
map gene names from multiple datasets to a single set of stable IDs (Ensemble) 
combine the list of housekkeping gene names, ribosomal genes, and other annotation datasets to keep accurate comparisons
'''

parser = argparse.ArgumentParser(description='filter genes')
parser.add_argument('-gv', type=str, choices=['19', '24', '26', '28', '35','45', '49', 'scGTEX', 'sccovid_multiomics', 'ensembl_v93','ensembl_v75', 'ensembl_v100'], help='' )
args = parser.parse_args()


bm = gseapy.Biomart()
# bm.get_datasets() 
# bm.get_marts()
# bm.get_attributes()

v19_p = 'gencode.v19.annotation.gff3' # housekeeping genes erez levanon, iPSC
v24_p = 'gencode.v24.annotation.gff3' # Tcell
v26_p = 'gencode.v26.annotation.gff3' # bulk_GTEx & sc_GTEx
v28_p = 'gencode.v28.annotation.gff3' # Ulicevic 2024
v35_p = 'gencode.v35.annotation.gff3' # first trimetster brain development 
v45_p = 'gencode.v45.annotation.gff3' # Cerebellum Atlas (single cell)
v49_p = 'gencode.v49.annotation.gff3' # GSEA 2025


if args.gv == 'scGTEX':
    h5ad = sc.read_h5ad('/lustre1/home/martinmikl/sshanas/scGTEx/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad')
    annot = h5ad.var
    annot = annot[(~annot.Chromosome.isin(['Y', 'X', 'M'])) & (annot.gene_biotype == 'protein_coding')].reset_index().set_index('gene_ids')
    annot.rename(columns={'index': 'gene_name'}, inplace=True)
    del h5ad

elif args.gv == 'sccovid_multiomics':
    h5ad = sc.read_h5ad('/lustre1/home/martinmikl/sshanas/sc_covid_multiomics/covid_portal_210320_with_raw.h5ad')    
    # annot = h5ad.var
    h5ad.var.to_csv(f'gene_names_to_id.{args.gv}.temp.csv.gz', compression='gzip')
    
    sys.exit()
    # annot = annot[(~annot.Chromosome.isin(['Y', 'X', 'M']))].reset_index().set_index('gene_ids')
    # annot.rename(columns={'index': 'gene_name'}, inplace=True)
    # del h5ad

elif args.gv == 'ensembl_v93':
    ''' this is for COVID vax (covid2) dataset, which was mapped to Ensembl v93, so we need to use the corresponding annotation file '''
    annot = gffpd.read_gff3('Homo_sapiens.GRCh38.93.gff3.gz').filter_feature_of_type(['gene']).attributes_to_columns()[['seq_id', 'gene_id', 'Name', 'biotype']]
    annot = annot[(annot.biotype=='protein_coding') & (~annot.seq_id.isin(['X', 'Y', 'M']))].rename(columns={'Name': 'gene_name', 'biotype': 'gene_type'}).set_index('gene_id')

elif args.gv == 'ensembl_v100':
    ''' this is for sc_organogenesis dataset '''
    annot = gffpd.read_gff3('Homo_sapiens.GRCh38.100.gff3.gz').filter_feature_of_type(['gene']).attributes_to_columns()[['seq_id', 'gene_id', 'Name', 'biotype']]
    annot = annot[(annot.biotype=='protein_coding') & (~annot.seq_id.isin(['X', 'Y', 'M']))].rename(columns={'Name': 'gene_name', 'biotype': 'gene_type'}).set_index('gene_id')

elif args.gv == 'ensembl_v75':
    ''' this is for differentiation to endoderm (single cell) dataset, which was mapped to Ensembl v75, so we need to use the corresponding annotation file '''
    annot = gffpd.read_gff3('Homo_sapiens.GRCh37.75.gtf.gz')


    def parse_gtf_attributes(attr_str):
        attrs = {}
        for field in attr_str.strip().split(';'):
            field = field.strip()
            if not field:
                continue
            parts = field.split(' ', 1)
            if len(parts) == 2:
                attrs[parts[0]] = parts[1].strip().strip('"')
        return attrs


    df = annot.df  # get the underlying dataframe
    attr_df = df['attributes'].apply(parse_gtf_attributes).apply(pd.Series)
    annot = pd.concat([df.drop(columns='attributes'), attr_df], axis=1)

    annot = annot[(annot.gene_biotype=='protein_coding') & (~annot.seq_id.isin(['X', 'Y', 'M'])) & (annot['type'] == 'gene')].rename(columns={'gene_biotype': 'gene_type'}).set_index('gene_id')
    
elif (args.gv != 'scGTEX') and (args.gv != 'sccovid_multiomics') and (args.gv != 'ensembl_v93') and (args.gv != 'ensembl_v75') and (args.gv != 'ensembl_v100'):
    annot = gffpd.read_gff3(f'gencode.v{args.gv}.annotation.gff3').filter_feature_of_type(['gene']).attributes_to_columns()[['seq_id', 'gene_id', 'gene_name', 'gene_type']]
    annot['accession_number'] = annot.gene_id.str.split('.', expand=True)[0]
    annot = annot[(annot.gene_type=='protein_coding') & (~annot.seq_id.isin(['chrX', 'chrY', 'chrM']))]

    ''' count exons and gene length '''
    exon_annotation = annot[annot['type']=='exon']
    exon_annotation['gene_name_start_end'] = exon_annotation.apply(lambda x: '_'.join(sorted([x['gene_name'], str(x['start']), str(x['end'])])), axis='columns')
    exon_annotation.drop_duplicates('gene_name_start_end', inplace=True)
    exon_annotation = exon_annotation.merge(pd.Series(exon_annotation.gene_name.value_counts(), name = 'nexons'), right_index=True, left_on='gene_name').drop_duplicates('gene_name')
    annot = annot.merge(exon_annotation[['gene_name', 'nexons']], on = 'gene_name').drop(columns=['transcript_type'])
    annot = annot.assign(length = annot['end'] - annot['start'])
    annot.set_index('accession_number', inplace=True)    


''' get information on housekeeping genes '''
house_keeping_genes_refseq = pd.read_csv("housekeeping_genes.txt", header=None, sep='\t')[1].tolist()

chunks = np.array_split(house_keeping_genes_refseq, 8)
hk_df_refseq = pd.DataFrame(columns = ['ensembl_gene_id', 'external_gene_name', 'refseq_mrna'])
for i, chunk in enumerate(chunks):
    try:
        results = bm.query(dataset='hsapiens_gene_ensembl',
                        attributes=['ensembl_gene_id', 'external_gene_name', 'external_synonym', 'refseq_mrna'], # , 'go_id'
                        filters={'refseq_mrna':chunk})
        hk_df_refseq = pd.concat([hk_df_refseq, results], axis = 'rows')
    except:
        continue    
housekeeping_ids = list(hk_df_refseq.ensembl_gene_id.unique())


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



annot['HK'] = annot.index.isin(housekeeping_ids)
annot['RP'] = annot.index.isin(RP_ids)
annot['HK_RP_intersect'] = annot.index.isin(np.intersect1d(housekeeping_ids, RP_ids))

HK_only = [g for g in housekeeping_ids if g not in RP_ids]
RP_only = [g for g in RP_ids if g not in housekeeping_ids]

annot['HK_only'] = annot.index.isin(HK_only)
annot['RP_only'] = annot.index.isin(RP_only)

annot.to_csv(f'gene_names_to_id.{args.gv}.csv.gz', compression='gzip', index=True)

# %%
