# %%
''' How input_files/ribosomal_genes.csv was made (kept for documentation).

Ribosome-related genes = every gene of the Enrichr library GO_Cellular_Component_2025 (human)
whose term name contains "ribosom", mapped to Ensembl gene IDs through Ensembl Biomart.
The list was retrieved once and saved, because both online resources change over time:
re-running this script today may return a slightly different list. The pipelines read the saved file.
'''
import pandas as pd
import gseapy

go_cc = gseapy.parser.download_library('GO_Cellular_Component_2025', organism='Human')
ribosomal_sets = {k: v for k, v in go_cc.items() if 'ribosom' in k.lower()}
ribosomal_genes = list(set([gene for genes_list in ribosomal_sets.values() for gene in genes_list]))

RP_results = gseapy.Biomart().query(dataset='hsapiens_gene_ensembl',
                attributes=['ensembl_gene_id', 'external_gene_name', 'external_synonym'],
                filters={'external_gene_name': ribosomal_genes})
RP_ids = sorted(RP_results.ensembl_gene_id.unique())

# the saved file keeps the protein-coding autosomal genes (GENCODE v28) among RP_ids, with their gene names
pd.Series(RP_ids, name='ensembl_gene_id').to_csv('ribosomal_genes.downloaded_now.csv', index=False)
