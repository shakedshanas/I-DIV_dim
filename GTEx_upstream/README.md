# GTEx_upstream: from GTEx v8 matrices to diversity matrices

Step 2 of the pipeline. Needs the GTEx v8 and GENCODE v26 files in `../input_files/`.
The outputs of step 10 are read by `../GTEx_pre_downstream/` and `../downstream/`.

Run everything from this folder. Before the first job, create the folder for the SLURM logs:

```bash
cd GTEx_upstream
mkdir -p outputs
```

Wait for each step to finish (`squeue -u $USER`) before starting the next.

| Step | Command | What it does | Main output |
|------|---------|--------------|-------------|
| 1 | `bash 1.activate.make_metadata_files.sh` | Gene annotation from GENCODE v26: protein-coding genes, not on chrX/Y/M; number of exons, length, ribosomal flag (`../input_files/ribosomal_genes.csv`). Sample metadata: death class 1 or 2, RIN ≥ 6 (missing RIN and ischemic time imputed by the SMTSD mean), SMTSD with ≥ 15 samples, samples present in both TPM matrices; age groups 20–49, 50–59, 60–79 | `1.outputs/gene_annotation.csv.gz`, `1.outputs/sample_metadata.csv.gz` |
| 2 | `bash 2.activate_firstfilter.sh` | Keeps the retained genes and samples in the gene and transcript TPM matrices (2 jobs) | `2.outputs/gene_tpm.csv.gz`, `2.outputs/transcript_tpm.csv.gz` |
| 3 | `sbatch 3.activate_calculate_div.sh` | I-DIV (1% spreading) and number of expressed isoforms for every gene and sample; one array task per gene (tasks 0–10000, each also takes gene + 10000) | `3.outputs/div1.csv`, `3.outputs/nisoforms_expressed.csv` |
| 4 | `bash 4.activate_genefilter.sh` | Keeps genes with I-DIV > 0 in more than 1% of samples (15,834 genes); z-scores gene TPM per sample | `4.outputs/` (annotation, `div1`, `nisoforms_expressed`, `gene_tpm`, `gene_tpm.zs`, `transcript_tpm`) |
| 5 | `bash 5.activate_calculate_transcriptome_div.sh` | Per sample: T-DIV on non-ribosomal genes (`div`), on ribosomal genes (`RP_div`) and on all genes (`15834_div`), number of expressed genes, µ I-DIV (`meanisofdiv`), mean NIS; Shannon diversity of the Donovan 2020 cell-type estimates and its Pearson correlation with T-DIV per SMTSD; tissue abbreviations | `4.outputs/metadata.csv.gz`, `4.outputs/CDIV_TDIV_pearosnr.csv.gz` |
| 6 | `sbatch 6.qn.activate_qn.sh` | Quantile normalisation within each tissue, separately for SMTS and SMTSD: log2(TPM+1) for genes, linear TPM for transcripts. One task per tissue (tasks 0–47; indices beyond the 25 SMTS / 47 SMTSD exit with an error and can be ignored) | `qn_data/{SMTS,SMTSD}/<tissue>.{G,I}.csv.gz` |
| 7 | `bash 7.qn.activate_jointables.sh` | Joins the per-tissue matrices into one matrix per level: all SMTS, all SMTSD, and the 13 brain SMTSD | `QN.{SMTS,SMTSD,Brain}.{G,I}.csv.gz` |
| 8 | `bash 8.make_headers.sh` | Writes the header line of the I-DIV matrices of step 9 | `QN.{SMTS,SMTSD,Brain}.IDIV.csv` (header only) |
| 9 | `sbatch 9.qn.activate_calc_qndiv.sh` | I-DIV from the quantile-normalised transcript matrices, one array task per gene | rows appended to `QN.{SMTS,SMTSD,Brain}.IDIV.csv` |
| 10 | `bash 10.clean.sh` | Moves the step-4/5 outputs here, deletes the intermediate folders and logs, gzips the CSV files | see below |

The array jobs of steps 3 and 9 append one row per gene to a shared CSV, so the row order differs between
runs. Later steps index the matrices by gene ID, so the order does not matter.

## Outputs (after step 10)

| File | Content |
|------|---------|
| `gene_annotation.csv.gz` | the 15,834 retained genes and their features |
| `metadata.csv.gz` | the 4,814 retained samples with covariates, T-DIV, µ I-DIV and cell-composition diversity |
| `gene_tpm.csv.gz`, `gene_tpm.zs.csv.gz` | gene TPM, and TPM z-scored per sample |
| `transcript_tpm.csv.gz` | transcript TPM of the retained genes |
| `div1.csv.gz` | I-DIV, genes × samples |
| `nisoforms_expressed.csv.gz` | number of isoforms with TPM > 0, genes × samples |
| `CDIV_TDIV_pearosnr.csv.gz` | Pearson r of T-DIV with cell-composition diversity, per SMTSD |
| `QN.<level>.G.csv.gz`, `QN.<level>.I.csv.gz`, `QN.<level>.IDIV.csv.gz` | quantile-normalised gene log2(TPM+1), transcript TPM and I-DIV; `<level>` = `SMTS`, `SMTSD` or `Brain` |

## Notes

- `3.calculate_div.py` and `9.qn.calc_qndiv.py` import `DIV()` from `../div.py`.
- The `#SBATCH --partition` values are those of our cluster; change them to yours.
