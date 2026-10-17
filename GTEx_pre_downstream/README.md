# GTEx_pre_downstream: statistics behind Figures 1–4

Step 3 of the pipeline. Every subfolder reads the outputs of `../GTEx_upstream/` (after its step 10)
and writes tables that the notebooks in `../downstream/` read. The subfolders are independent of each
other and can run in any order or at the same time.

Run each command from inside its subfolder. Jobs whose logs go to `outputs/` need that folder first:
`mkdir -p outputs`.

## fig1_idiv_tpm_pearson: I-DIV, expression and number of isoforms (Fig. 1, S2)

```bash
cd fig1_idiv_tpm_pearson && mkdir -p outputs
sbatch activate.sh
```

Pearson correlation for every gene (one array task per gene):

- per SMTS tissue: I-DIV with log2(TPM+1), and I-DIV with the number of expressed isoforms (NIS);
- all samples pooled, the six sex-specific tissues (uterus, breast, testis, ovary, vagina, prostate)
  excluded (`-b all_wo_sext`): I-DIV with NIS, NIS with log2(TPM+1), I-DIV with log2(TPM+1).

Only samples where the first measure is > 0 (I-DIV or NIS) are used; with fewer than 10 such samples the
gene gets r = 0, p = 1.
Outputs: `div_tpm/SMTS.pearson.csv`, `div_nis/SMTS.pearson.csv` (per tissue `<SMTS>_coeff`, `<SMTS>_pval`)
and `{div_nis,nis_tpm,div_tpm}/all_wo_sext.pearson.csv` (`coeff`, `pval`).

## fig2_tissue_pairs: differences between tissues (Fig. 2, S3)

```bash
cd fig2_tissue_pairs
bash 1.makedirs.sh                      # SMTS_limma/, SMTSD_limma/, outputs/
sbatch 2.activate_limma.SMTS.sh         # 300 SMTS pairs
sbatch 2.activate_limma.SMTSD.sh        # 1081 SMTSD pairs
# when both have finished:
sbatch 3.activate_limma_process_res.sh
```

- `DE_limma.R` (step 2), per tissue pair, separately for log2(TPM+1) and for I-DIV: limma moderated
  t-test (`lmFit`, `eBayes`, BH FDR) and a two one-sided equivalence test on the same fit with margin
  δ = 0.5. A gene is *differentially expressed* at FDR < 0.05, *equivalent* when the equivalence FDR < 0.05
  and the 90% CI of the difference lies within ±δ, and *inconclusive* otherwise.
  Output: `<level>_limma/<pair>.<tissue1>.<tissue2>.{tpm,div}.limma_DE_equiv.csv.gz`.
- `3.limma_process_res.py` (step 3), per pair: PCA (2 components, standardised features) and the
  silhouette score of the two tissues, using all genes, the TPM-equivalent genes or the I-DIV-equivalent
  genes. Output: one row per pair in `<level>_limma/overview_cluster_results.csv`.

R needs `TMPDIR` on a disk with space; the job scripts set it to `$SCRATCH_DIR` (default `/tmp`).

## fig3_tdiv_tpm_pearson: transcriptome diversity and gene expression (Fig. 3, S4)

```bash
cd fig3_tdiv_tpm_pearson
sbatch 3.activate_cor_sample_metrices.sh
```

For every gene, the Pearson correlation of its expression (TPM z-scored per sample) with T-DIV
(`T-DIV.py`) and with µ I-DIV (`mI-DIV.py`), across all samples, per SMTSD and per SMTS. Only samples with
TPM > 0 for the gene are used, tissues need more than 10 samples, and a correlation needs ≥ 10 samples
(otherwise r = 0, p = 1). For T-DIV, samples below the 1st percentile of the number of expressed genes are
excluded. Outputs: `pearson.{T-DIV,mI-DIV}.{all,SMTSD,SMTS}.csv`.

## fig4_spearman: Spearman correlation with age per tissue (Fig. 4, S5)

```bash
cd fig4_spearman && mkdir -p outputs
sbatch 5.activate_spearman.sh
```

Per tissue with more than 10 samples (SMTS and SMTSD), for every gene: Spearman correlation of age group
with gene expression and with I-DIV, from the within-tissue quantile-normalised matrices
`QN.<level>.{G,IDIV}.csv.gz`. Outputs: `<level>/<tissue>.spearmanr.{tpm,div}.csv.gz` (and `.raw.` input tables).

## fig4_mlr: limma regression of age per tissue (Fig. 4, S5)

```bash
cd fig4_mlr/per_gene && mkdir -p outputs
sbatch activate.limma.sh                # one task per tissue, SMTSD and SMTS
cd ../per_isoform
bash commands.sh                        # PPARD, FAM161B, KNDC1
```

- `per_gene/limma_age_per_tissue.R`: per tissue (> 10 samples), y ~ age + sex + RIN + ischemic time
  (+ SMTSD for multi-subtissue SMTS), with age coded 1, 2, 3 for 20–49, 50–59, 60–79; `eBayes(trend = TRUE)`,
  BH FDR per term. y is log2(TPM+1) or I-DIV, both quantile-normalised within SMTSD (`QN.SMTSD.*`) for
  both levels. Output: `<level>/<tissue>.limma.{tpm,div}.coefficients.csv.gz`.
- `per_isoform/isoform_age_mlr_limma.R`: the same model per isoform of one gene, in each brain SMTSD,
  on the quantile-normalised transcript TPM (`QN.Brain.I.csv.gz`).
  Output: `<gene>/<gene>.brain.isoform.limma.MLR.csv.gz`.

## fig4_lme: mixed models of age across tissues (Fig. 4, S5)

```bash
cd fig4_lme
bash activate.lme.sh
```

Per gene, `nlme::lme` (REML) with a random intercept per donor, separately for log2(TPM+1) and I-DIV:

| Script | Samples | Fixed effects | Output folder |
|--------|---------|---------------|---------------|
| `lme.across_SMTS.R` | all SMTS (`QN.SMTS.*`) | SMTS + sex + age + RIN + ischemic time | `lme_across_SMTS/` |
| `lme_across_SMTS_exRINISCH.R` | all SMTS (`QN.SMTS.*`) | SMTS + sex + age | `lme_across_SMTS_exRINISCH/` |
| `lme.across_SMTSDbrain.R` | 13 brain SMTSD (`QN.Brain.*`) | SMTSD + sex + age + RIN + ischemic time | `lme_across_brainSMTSD/` |

Each folder holds a summary table (age and sex coefficients with BH FDR), all coefficients in long and
wide format, per-gene summaries and residuals, and the fitted models (`.rds`).
