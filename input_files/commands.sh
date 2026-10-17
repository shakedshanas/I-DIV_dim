#!/bin/bash
# Download all public input data into this directory (run it from input_files/).

### download GTEx data V8 ###
# matrices: #
wget -nc https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz
wget -nc https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz

# metadata: #
wget -nc https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
wget -nc https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt

# GENCODE version 26 (the annotation GTEx v8 was quantified with) #
wget -nc https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gff3.gz
gunzip -k gencode.v26.annotation.gff3.gz

### Ulicevic data ###
# Mass-Spec (PRIDE PXD046084): #
wget -nc https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/04/PXD046084/isoform_txt.zip
# RNA-seq (GEO GSE245325), RSEM *.genes.results.gz and *.isoforms.results.gz per sample: #
wget -nc https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245325/suppl/GSE245325%5FRAW%2Etar
tar -xvf GSE245325_RAW.tar

# GENCODE version 28 #
wget -nc https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gff3.gz
gunzip -k gencode.v28.annotation.gff3.gz

### supplementary tables of published papers (see README.md) ###
# Donovan2020/41467_2020_14561_MOESM{5..18}_ESM.csv  - included in this repository (CC BY 4.0): Donovan et al. 2020, Nat Commun 11:955, Supplementary Data 4-17 (MOESM5 = Supp. Data 4)
# mele2015.AgeDiffExpGenes.S12.csv                   - download by hand: Mele et al. 2015, Science 348:660, Table S12
