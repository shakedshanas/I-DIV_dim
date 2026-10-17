# Transcriptome and isoform diversity in GTEx

Code for the analyses and figures of the paper *[ADD: title, authors, journal, year, DOI]*.

## The idea

A gene can be expressed as one dominant isoform or as many isoforms at similar levels.
We measure this with **isoform diversity (I-DIV)**: for one gene in one sample, the Shannon evenness
of its isoform TPM values (natural log, normalised by the log of the number of annotated isoforms),
after adding 1% of the gene's total TPM spread evenly over its isoforms. I-DIV is 0 when at most one
isoform is expressed and approaches 1 when all isoforms are expressed equally.

Applying the same index to the gene-level TPM vector of a whole sample gives **transcriptome diversity
(T-DIV)**; T-DIV is computed over the non-ribosomal genes. The mean I-DIV of a sample over all genes is **µ I-DIV**. The function is `DIV()` in
[`div.py`](div.py).

We compute I-DIV, T-DIV and µ I-DIV for 15,834 protein-coding genes in 4,814 GTEx v8 samples from 267
donors (47 detailed tissues, SMTSD; 25 tissues, SMTS; samples filtered for death class 1–2, RIN ≥ 6 and
tissues with ≥ 15 samples) and ask how isoform diversity relates to

- gene expression and the number of expressed isoforms (Fig. 1),
- differences between tissues: limma differential and equivalence tests per tissue pair (Fig. 2),
- the cell-type composition of the tissue (Donovan et al. 2020) (Fig. 3),
- donor age: Spearman correlation, limma regression per tissue and mixed models across tissues (Fig. 4),

and validate the approach in a neuronal differentiation time course with matched RNA-seq and proteomics
(Ulicevic et al. 2024) (Fig. 5).

## Repository layout and the order to run it

```
input_files/          1. download the public data (GTEx v8, GENCODE, Ulicevic et al. 2024)
GTEx_upstream/        2. GTEx: filter samples and genes, I-DIV / T-DIV, quantile normalisation (steps 1-10)
GTEx_pre_downstream/  3. GTEx: statistics per figure (correlations, limma, Spearman, mixed models)
Ulicevic_upstream/    4. Ulicevic et al. 2024: combine samples, I-DIV, isoform-level limma
downstream/           5. notebooks that make the figures, source data and supplementary tables
demo/                 a small self-contained example of the I-DIV and T-DIV calculation
div.py                the diversity function used by every step
```

Steps 2 → 3 → 5 must run in this order. Step 4 only needs step 1 and can run at any time before
`downstream/Ulicevic.ipynb`. Every folder has its own README with the exact commands.

| Figure | Statistics (folder 3 or 4)                                     | Notebook (folder 5)                  |
|--------|----------------------------------------------------------------|--------------------------------------|
| 1, S2  | `GTEx_pre_downstream/fig1_idiv_tpm_pearson`                    | `downstream/intro.ipynb`             |
| 2, S3  | `GTEx_pre_downstream/fig2_tissue_pairs`                        | `downstream/limma_equivalence.ipynb` |
| 3, S4  | `GTEx_pre_downstream/fig3_tdiv_tpm_pearson`                    | `downstream/transcriptome_div.ipynb` |
| 4, S5  | `GTEx_pre_downstream/fig4_spearman`, `fig4_mlr`, `fig4_lme`    | `downstream/lme.ipynb`               |
| 5, S6  | `Ulicevic_upstream`                                            | `downstream/Ulicevic.ipynb`          |
| S1     | none                                                           | `downstream/metrics_comparisons.ipynb` |

## System requirements

- **Operating system:** Linux. The analysis was run on Red Hat Enterprise Linux 9.1 on an HPC cluster.
  The demo also runs on macOS or Windows with the Python environment below.
- **Scheduler:** the GTEx steps are SLURM array jobs (`sbatch`), one task per gene, tissue or tissue pair.
  The partition names in the `#SBATCH` lines (`hive1d`, `hive7d`, `queen`) are those of our cluster;
  change them to yours.
- **Hardware:** a standard HPC cluster; no non-standard hardware. The GTEx input matrices are large
  (the transcript TPM file is ~200,000 transcripts × ~17,000 samples), so the steps that load a whole
  matrix (`GTEx_upstream` steps 1, 2, 4, 5, 7 and the R models) need a node with tens of GB of RAM;
  the quantile-normalisation tasks request 10 GB (`--mem=10G`). Memory use was not measured.
  The demo runs on a laptop.
- **Software:**
  - Python 3.12.2 with the packages and versions in [`environment.yml`](environment.yml)
  - R 4.3.3 with limma 3.58.1, nlme 3.1.171, data.table 1.18.6.1, optparse 1.8.2, R.utils 2.13.0
    ([`R_requirements.txt`](R_requirements.txt))

## Installation

```bash
git clone https://github.com/shakedshanas/transcriptome_and_isoform_diversity_GTEx.git
cd transcriptome_and_isoform_diversity_GTEx
conda env create -f environment.yml
conda activate isoform_div
# R packages, inside R 4.3.3:
#   install.packages(c("data.table", "nlme", "optparse", "R.utils"))
#   install.packages("BiocManager"); BiocManager::install("limma")
```

Creating the conda environment usually takes 5–20 minutes, depending on the network (not timed).

## Demo

```bash
cd demo
jupyter nbconvert --to notebook --execute demo.ipynb --output demo.executed.ipynb
```

It computes I-DIV for the genes and T-DIV for the 9 GTEx samples in `demo/*.csv`, and plots the isoform
abundances behind them. It took 41 seconds on a login node of our cluster. See [`demo/README.md`](demo/README.md) for the
expected output.

## Instructions for use

1. `input_files/`: run `bash commands.sh` and add the two supplementary tables named in its README.
2. `GTEx_upstream/`: run steps 1–10 in order (README there).
3. `GTEx_pre_downstream/`: run the folder of each figure (README there).
4. `Ulicevic_upstream/`: run steps 1, 2 and the isoform limma script (README there).
5. `downstream/`: run `bash commands.sh`, or open each notebook and run it from top to bottom.

All paths are relative to the repository, so the folders must keep their names and positions.
Run every script from its own folder (`cd` into it first).

## Data

Raw data are public and are not included: GTEx v8 (dbGaP phs000424.v8, open-access TPM matrices and
annotations), GENCODE v26 and v28, GEO GSE245325 and PRIDE PXD046084 (Ulicevic et al. 2024).
The repository includes three small files: the fixed list of ribosome-related genes
(`input_files/ribosomal_genes.csv`), the cell-composition estimates of Donovan et al. 2020
(`input_files/Donovan2020/`, CC BY 4.0), and the tissue colour palettes (`downstream/palettes.json`).

## Nature Research code and software submission checklist

| Checklist item                                     | Where                                              |
|----------------------------------------------------|----------------------------------------------------|
| Source code                                        | this repository                                    |
| Compiled standalone software / installation guide  | not applicable (scripts); see Installation         |
| System requirements, dependencies with versions    | System requirements, `environment.yml`, `R_requirements.txt` |
| Tested on                                          | RHEL 9.1, Python 3.12.2, R 4.3.3                   |
| Typical install time                               | Installation                                       |
| Demo: data, instructions, expected output, run time | `demo/`, `demo/README.md`                         |
| Instructions for use on the paper's data           | Instructions for use and the README of each folder |
| Open-source license                                | MIT, [`LICENSE`](LICENSE)                          |
| Repository with a persistent identifier            | *[ADD: Zenodo DOI of the release used in the paper]* |

## Citation

*[ADD: citation of the paper.]*

Data and methods used: GTEx Consortium (2020) *Science* 369:1318; Donovan et al. (2020)
*Nat Commun* 11:955; Melé et al. (2015) *Science* 348:660; Ulicevic et al. (2024), PMID 38755290.
