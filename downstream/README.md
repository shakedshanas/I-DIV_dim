# downstream: figures, source data and supplementary tables

Step 5 of the pipeline. One Jupyter notebook per figure. The notebooks read the outputs of
`../GTEx_upstream/`, `../GTEx_pre_downstream/` and `../Ulicevic_upstream/`, and import `DIV()` and helper
functions from `../div.py` (and, for `lme.ipynb`, from `fig4funcs.py` in this folder).

## Run

All notebooks, as SLURM jobs:

```bash
cd downstream
bash commands.sh
```

or one notebook interactively: `jupyter lab`, open it, *Run All*. Notebooks must be run from this folder.
`commands.sh` first creates `source_data/`, `supplementary_data/` and `tables_main/`; create them by hand
if you run a notebook interactively.

| Notebook | Figures | Needs (besides `../GTEx_upstream/`) |
|----------|---------|-------------------------------------|
| `metrics_comparisons.ipynb` | S1: I-DIV compared with other isoform-diversity indices (`isoplex`) on simulated data | nothing (simulated data) |
| `intro.ipynb` | 1, S2: I-DIV examples, I-DIV vs TPM and NIS, k-means clusters of genes and their enrichment | `../GTEx_pre_downstream/fig1_idiv_tpm_pearson/` |
| `limma_equivalence.ipynb` | 2, S3: tissue-pair silhouette scores, differential and equivalent genes | `../GTEx_pre_downstream/fig2_tissue_pairs/` |
| `transcriptome_div.ipynb` | 3, S4: T-DIV and µ I-DIV per tissue, cell-composition diversity, gene correlations with T-DIV | `../GTEx_pre_downstream/fig3_tdiv_tpm_pearson/`, `../input_files/Donovan2020/`, `palettes.json` |
| `lme.ipynb` | 4, S5: age effects (Spearman, limma, mixed models), comparison with Melé et al. 2015, isoform examples | `../GTEx_pre_downstream/fig4_*/`, `../input_files/mele2015.AgeDiffExpGenes.S12.csv`, `fig4funcs.py` |
| `Ulicevic.ipynb` | 5, S6: I-DIV and T-DIV over neuronal differentiation, proteomics | `../Ulicevic_upstream/`, `../input_files/isoform_txt.zip` |

Gene-set enrichment uses the Enrichr web service through `gseapy.enrichr`, so these cells need internet
access, and their results depend on the library versions on the Enrichr server at the time of running
(libraries named in the notebooks, e.g. GO_Biological_Process_2025).

## Outputs

Each notebook writes its figures to `figs_<name>/` and its intermediate tables to `tables_<name>/`.
Source data of every figure panel go to `source_data/<panel>.csv.gz` (e.g. `1b.csv.gz`, `S5c_top.csv.gz`),
and supplementary tables to `supplementary_data/`.

`palettes.json` holds the fixed colours of the 25 SMTS and 47 SMTSD tissues (RGBA), so that every figure
uses the same colour for a tissue.
