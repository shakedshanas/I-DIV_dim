# Isoform Diversity — Differentiation Time Course

Scripts for processing RSEM output from a time-course RNA-seq experiment (samples named `day{N}_{replicate}`, e.g. `day0_1`) into per-sample isoform/gene expression matrices and per-gene isoform diversity scores. Unlike the GTEx cluster pipeline, these are run locally/interactively (no SLURM job arrays, aside from the proteomics download) directly on a folder of RSEM result files.

## Scripts

| Script | Purpose |
|---|---|
| `commands_for_proteomics.sh` | Downloads an orthogonal proteomics isoform dataset (PRIDE project PXD046084) via a SLURM job |
| `combine_isoform_data.py` | Combines all isoform-level RSEM TPM results into a single matrix, filtered to a gene list |
| `make_isoform_data.py` | Builds separate per-sample TPM and transcript-length matrices from isoform-level RSEM results |
| `process.py` | Computes per-gene isoform diversity (and supporting stats) from gene- and isoform-level RSEM results |

## Input data

All three Python scripts expect RSEM output files in the working directory, one pair per sample:
- `*.genes.results.gz` — gene-level RSEM quantification (tab-separated, with `gene_id` and `TPM` columns)
- `*.isoforms.results.gz` — isoform-level RSEM quantification (tab-separated, with `transcript_id`, `gene_id`, `TPM`, and `length` columns)

Sample names are parsed out of the filenames (the `day{N}_{replicate}` portion) so that each sample becomes its own column in the combined output matrices.

## Script details

### `commands_for_proteomics.sh`
Submits a SLURM job to download `isoform_txt.zip` from PRIDE project [PXD046084](https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/04/PXD046084/), a proteomics dataset used alongside the RNA-seq isoform data.

### `combine_isoform_data.py`
Reads every `*isoforms.results.gz` file in the current directory, restricts to genes present in an external `genes.csv` gene list, and joins each sample's isoform TPM values into one wide table (one row per transcript, one column per sample). Writes `isoforms_tpm_combined.csv`.

> Note: the gene list path is currently hardcoded to a local Google Drive path (`Ulicevic2024/RNAseq/genes.csv`) — update this before reuse.

### `make_isoform_data.py`
Similar to `combine_isoform_data.py` but builds two separate matrices by merging all samples' isoform-level results on `transcript_id` + `gene_id`: one for `TPM`, one for transcript `length`. Writes `TPM.isoform.csv` and `length.isoform.csv`.

### `process.py`
For each gene, groups its isoforms across all samples and computes, using the shared `DIV`/`DIV_min` functions (from `div.py`, the same diversity module used in the main GTEx pipeline):
- `evenness0` / `evenness1` — isoform diversity with 0% and 1% spreading
- `evenness1_norm` / `evenness1_norm2` — diversity normalized against the gene's theoretical minimum diversity
- `sum_pergene` — total TPM per gene per sample
- `sum_nisoforms` — number of expressed isoforms per gene per sample

Also writes gene-level (`tpm.csv`) and isoform-level (merged across samples, TPM only) intermediate tables. Each resulting dataframe is saved to its own CSV named after the dictionary key (e.g. `evenness1.csv`, `sum_pergene.csv`).

## Requirements

Same Python environment as the main GTEx pipeline (`pandas`, `numpy`); `process.py` additionally requires the `div.py` module to be importable (update the hardcoded `sys.path.append(...)` to wherever `div.py` lives in your checkout).
