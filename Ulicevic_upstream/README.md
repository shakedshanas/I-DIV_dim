# Ulicevic_upstream: neuronal differentiation time course (Ulicevic et al. 2024)

Step 4 of the pipeline. Needs, in `../input_files/`, the RSEM files of GEO GSE245325
(`*.genes.results.gz`, `*.isoforms.results.gz`) and `gencode.v28.annotation.gff3`.
Independent of the GTEx steps. Its outputs are read by `../downstream/Ulicevic.ipynb`.

Run from this folder, in order, waiting for each job to finish:

| Step | Command | What it does | Output |
|------|---------|--------------|--------|
| 1 | `bash 1.activate.conbine.data.sh` | Protein-coding autosomal genes of GENCODE v28 with the ribosomal flag (`../input_files/ribosomal_genes.csv`); merges the per-sample RSEM TPM files into one gene and one transcript matrix. Sample names are taken from the file names (`day<N>_<replicate>`) | `gene_annotation.csv.gz`, `gene_tpm.csv.gz`, `transcript_tpm.csv.gz` |
| 2 | `bash 2.activate.calculate_div.sh` | I-DIV (1% spreading) and number of expressed isoforms for every gene and sample | `div1.csv.gz`, `nisoforms_expressed.csv.gz` |
| 3 | `bash activate.limma_d0_d5.per_isoform.sh` | For DLL3, DHX15 and RNPS1: limma on log2(TPM+1) of each isoform, day 0 vs day 5 (BH FDR across the isoforms of the gene) | `isoform_day0_day5_de/<gene>/<gene>.day0.day5.isoform.limma.csv.gz` |

Step 3 runs `module load R/4.3.3` as on our cluster; replace it with however R 4.3.3 is made available
on yours. To test another gene: `Rscript limma_d0_d5.per_isoform.R --gene=<Ensembl gene ID>`.

The proteomics data (`../input_files/isoform_txt.zip`, PRIDE PXD046084) are read directly by the notebook.
