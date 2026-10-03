# GTEx Isoform Diversity Pipeline (upstream folder)

A SLURM-based Python pipeline for processing GTEx v8 gene/transcript expression data, computing per-gene isoform (splicing) diversity, correcting for technical covariates, and quantile-normalizing the results.

Each numbered step has a `.py` script and a matching `activate` `.sh` wrapper that submits it (or a family of related jobs) to a SLURM cluster.

## Pipeline overview

| Step | Script | Activation script | Purpose | Output dir |
|---|---|---|---|---|
| — | `commands.sh` | — | Downloads the GENCODE v26 GFF3 gene annotation used for filtering | — |
| 1 | `1.firstfilter.py` | `1.activate_firstfilter.sh` | Builds sample & gene metadata; filters samples (death type, RIN, tissue group size) and genes (protein-coding, non-sex/mito chromosomes); extracts the four raw GTEx matrices | `outputdir.1/` |
| 2 | `2.calculate_div.py` | `2.activate_calculate_div.sh` | Computes per-gene isoform diversity, total expression, and number of expressed isoforms, one gene at a time | `outputdir.2/` |
| 3 | `3.genefilter.py` | `3.activate_genefilter.sh` | Filters all downstream tables down to genes with detectable isoform diversity in ≥1% of samples | `outputdir.3/` |
| 4 | `4.correct_data.py` | `4.activate_correct_data.sh` | Regresses out technical covariates (batch, RIN, ischemic time, age, sex, death type) from expression and diversity data, globally or per tissue | `outputdir.4/` |
| 5 | `5.calculate_transcriptome_div.py` | `5.activate_calculate_transcriptome_div.sh` | Computes transcriptome-wide diversity/expression summary statistics per sample, for different housekeeping/ribosomal gene groupings | `outputdir.5/` |
| 6 | `6.quantile_norm.py` | `6.activate_quantile_norm.sh` | Quantile-normalizes the batch-corrected (residual) matrices from step 4 | `outputdir.6/` |
| — | `div.py` | — | Shared module: `DIV`/`DIV_min` (Shannon's evenness-based diversity index) and a row-skipping helper used across steps | — |

Run the steps in order (1 → 6); each step reads from the previous step's output directory.

## Requirements

- Access to a SLURM cluster (`sbatch`) with a `hive1d` partition (or edit the `#SBATCH` headers to match your cluster)
- Raw GTEx v8 files (gene reads, gene TPM, transcript reads, transcript TPM `.gct.gz` files) and the GTEx v8 sample/subject annotation files
- GENCODE v26 gene annotation (`gencode.v26.annotation.gff3.gz`, downloaded via `commands.sh`) and a derived `gene_names_to_id.26.csv.gz` mapping file

### Environment

All scripts are run inside the `spl_mpra_map` conda environment defined in [`BulkRNAseq.yaml`](BulkRNAseq.yaml), which pins Python 3.12 plus the packages the pipeline depends on (`pandas`, `numpy`, `scikit-learn`, `scipy`, `statsmodels` and Jupyter). Create and activate it with:

```bash
conda env create -f BulkRNAseq.yaml
conda activate BulkRNAseq
```

## Step details

### 0. `commands.sh`
Downloads the raw input data the pipeline expects one directory above the working directory:
- `gencode.v26.annotation.gff3.gz` (via `sbatch`) — GENCODE v26 human gene annotation, used as the basis for gene filtering in step 1
- `GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz` — transcript-level TPM
- `GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz` — transcript-level read counts
- `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz` — gene-level TPM
- `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz` — gene-level read counts
- `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt` — sample attributes (tissue, batch, RIN, ischemic time, etc.)
- `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` — subject phenotypes (sex, age, death type)

The four `GTEx_Analysis_*.gct.gz` files are the raw matrices `1.firstfilter.py` extracts from in step 1 (`--data 0-3`); the two `Annotations_*DS.txt` files are the sample/subject metadata `1.firstfilter.py` reads to build `sample_metadata.csv.gz`.

### 1. Sample & gene filtering — `1.firstfilter.py`
Run as a job array (`--array=0-3`), where the array index selects which raw matrix to extract:
- `0`: gene read counts
- `1`: gene TPM
- `2`: transcript read counts
- `3`: transcript TPM

On first run it also builds and caches `sample_metadata.csv.gz` and `gene_metadata.csv.gz` (filtering samples by death type (1 and 2), RNA quality (RIN ≥ 6), and tissue group size, and filtering genes to protein-coding, non-sex, non-mitochondrial genes present in all four raw matrices).

Outputs for this step are found at `outputdir.1/`.

### 2. Per-gene diversity calculation — `2.calculate_div.py`
Run as a large job array for each gene index, loads only that gene's row from the transcript-level matrix and computes:
- `div1`: isoform diversity (from `div.py`)
- `total_expression`: summed expression across isoforms
- `nisoforms_expressed`: count of isoforms with non-zero expression

Results are appended gene-by-gene into shared CSVs in `outputdir.2/`.

### 3. Gene filtering — `3.genefilter.py`
Submitted once per dataframe that needs filtering (TPM, gene metadata, diversity, isoform count, total expression). Keeps only genes with non-zero diversity in more than 1% of samples, reducing the gene set for downstream, more expensive steps.  
Outputs for this step are found at `outputdir.3/`.


### 4. Covariate correction — `4.correct_data.py`
Run as a job array (`--array=0-1`, mapped to `lgtype` 1 or 2) that in turn launches further `sbatch` jobs covering:
- `--dataname`: `readtpm_g` (log2-transformed TPM) or `div1` (isoform diversity)
- `--grouptype`: `i` (single regression across all samples) or `ii` (separate regression per tissue group, split by `SMTS` or `SMTSD`)
- `--lgtype`: `1` regresses out and adds back biological covariates (age, sex, death type) while removing only technical noise; `2` removes only technical covariates

Uses `scikit-learn`'s `LinearRegression` per gene against a design matrix of batch, RIN, ischemic time, and (for `lgtype=1`) age/sex/death type. Outputs residual and predicted matrices.

Outputs for this step are found at `outputdir.4/`.

### 5. Transcriptome-wide diversity — `5.calculate_transcriptome_div.py`
Run as a job array (`--array=0-6`), where each index selects a gene grouping (all genes, ribosomal protein genes, housekeeping genes, combinations/exclusions thereof). For each grouping, computes per-sample summary statistics — number of genes expressed, mean isoform diversity, mean TPM — overall and restricted to the top 1/5/10 percent most highly expressed genes.

Outputs for this step are found at `outputdir.5/`.

### 6. Quantile normalization — `6.quantile_norm.py`
Quantile-normalizes each batch-corrected residual matrix produced in step 4 (sample-group + tissue-specific residuals for TPM and diversity), so that each sample's value distribution matches a common reference distribution.

Outputs for this step are found at `outputdir.6/`.

## Directory structure produced

```
outputs/            # SLURM stdout/stderr logs
outputdir.1/         # filtered metadata + raw extracted matrices
outputdir.2/         # per-gene diversity, expression, isoform count
outputdir.3/         # gene-filtered versions of the above
outputdir.4/         # covariate-corrected (residual/predicted) matrices
outputdir.5/         # transcriptome-wide per-sample summary statistics
outputdir.6/         # quantile-normalized corrected matrices
```
