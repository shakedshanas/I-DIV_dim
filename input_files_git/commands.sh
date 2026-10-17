#!/bin/bash

### download GTEx data V8 ###
# matrices: #
wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz .
wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz

# metadata: #
wget https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt .
wget https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt .

# GENCODE version 26 # 
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gff3.gz .
gunzip gencode.v26.annotation.gff3.gz

### Ulicevic data ###
# Mass-Spec: #
sbatch --wrap "wget https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/04/PXD046084/isoform_txt.zip ." 
# RNA-seq: #
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245325/suppl/GSE245325%5FRAW%2Etar .
tar -xvf GSE245325_RAW.tar

# GENCODE version 28 # 
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gff3.gz .
gunzip gencode.v28.annotation.gff3.gz