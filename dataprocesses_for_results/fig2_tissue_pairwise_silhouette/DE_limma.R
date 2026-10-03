x#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# R translation of gtex pairwise-tissue t-test script.
#
# The per-gene independent t-test is replaced with:
#   1) a limma lmFit/eBayes moderated t-test for differential
#      expression (analogous to a two-group comparison), run
#      directly on log2(TPM+1) values (no voom/counts needed here,
#      since the input is already a continuous expression matrix);
#   2) a TOST-style equivalence test per gene (build on the same
#      lmFit object) to flag genes whose effect size is small enough
#      to call "non-differentially expressed" / equivalent, using a
#      practical-equivalence margin `delta`.
#
# Both tests are run once for the TPM matrix and once for the
# diversity matrix, mirroring the two t-tests in the original script.
#
# Note on indexing: R's combn()/list indexing is 1-based, Python's
# combinations()/list indexing is 0-based. If reusing the same array
# index as the Python job, add 1 to --pair.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(limma)
})

# ---------------------------- arguments ----------------------------
option_list <- list(
  make_option("--pair", type = "integer",
              help = "pair number (1-indexed) out of tissue_pairs: 1080 subtissue pairs (SMTSD) or 300 tissue pairs (SMTS) in GTEx v8"),
  make_option("--basedon", type = "character",
              help = "pairwise on the subtissue level (SMTSD) or the organ level (SMTS)")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$pair) || is.null(opt$basedon)) {
  stop("Both --pair and --basedon are required.")
}
if (!opt$basedon %in% c("SMTSD", "SMTS")) {
  stop("--basedon must be one of 'SMTSD' or 'SMTS'")
}

# ------------------------- metadata & pairs -------------------------
metadata_full <- fread(
  "/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.5/no_filter/metadata.csv.gz",
  select = c("SAMPID", opt$basedon)
)
metadata_full <- as.data.frame(metadata_full)
rownames(metadata_full) <- metadata_full$SAMPID

tissue_levels <- sort(unique(metadata_full[[opt$basedon]]))
tissue_pairs  <- combn(tissue_levels, 2, simplify = FALSE)

if (opt$pair < 1 || opt$pair > length(tissue_pairs)) {
  stop(sprintf("--pair must be between 1 and %d for --basedon %s",
               length(tissue_pairs), opt$basedon))
}
ppair <- tissue_pairs[[opt$pair]]

# --------------- load data: metadata, tpm, diversity matrix ---------------
metadata <- metadata_full[metadata_full[[opt$basedon]] %in% (ppair), , drop = FALSE] # metadata_full[metadata_full$SMTSD %in% ppair[[1]], , drop = FALSE]
metadata <- metadata[order(metadata[[opt$basedon]]), ]
paired_samples <- rownames(metadata)

d <- list(t1 = ppair[1], t2 = ppair[2]) 

tpm_path <- "/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/readtpm_g.csv.gz"
gtex_tpm <- fread(tpm_path, select = c("accession_number", paired_samples))
gtex_tpm <- as.data.frame(gtex_tpm)
rownames(gtex_tpm) <- gtex_tpm$accession_number
# gtex_tpm$accession_number <- NULL
gtex_tpm <- gtex_tpm[, paired_samples]
gtex_tpm <- log2(gtex_tpm + 1)

div_path <- "/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/div1.csv.gz"
# gtex_div <- fread(div_path, select = c("Unnamed: 0" = "V1", paired_samples))
gtex_div <- fread(div_path, select = c("Unnamed: 0", paired_samples))
gtex_div <- as.data.frame(gtex_div)
rownames(gtex_div) <- gtex_div[[1]]
# gtex_div[[1]] <- NULL
gtex_div <- gtex_div[, paired_samples]

# # make identical genes in both matrices, same order of samples
# gene_list <- intersect(rownames(gtex_tpm), rownames(gtex_div))
# print(length(gene_list))
# gtex_div <- gtex_div[gene_list, , drop = FALSE]
# print(length(colnames(gtex_div)))
# gtex_tpm <- gtex_tpm[gene_list, , drop = FALSE]
print(length(rownames(gtex_tpm)))
print(length(rownames(gtex_div)))
# ---------------------- condition / design matrix ----------------------
metadata$condition <- factor(metadata[[opt$basedon]], levels = ppair)
coef_name <- paste0("condition", ppair[2])
design <- model.matrix(~ condition, data = metadata)

# ------------- DE + equivalence test, run once per matrix -------------
run_de_and_equivalence <- function(expr_mat, design, coef_name, delta = 0.5) {

  fit <- lmFit(expr_mat, design)
  fit <- eBayes(fit)

  de_rs <- topTable(
    fit,
    coef = coef_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )

  estimate <- fit$coefficients[, coef_name]
  se <- fit$stdev.unscaled[, coef_name] * fit$sigma
  df <- fit$df.total

  p_lower <- pt((estimate + delta) / se, df = df, lower.tail = FALSE)
  p_upper <- pt((estimate - delta) / se, df = df, lower.tail = TRUE)
  p_equiv <- pmax(p_lower, p_upper)

  FDR_equiv <- p.adjust(p_equiv, method = "BH")

  CI90_lower <- estimate - qt(0.95, df) * se
  CI90_upper <- estimate + qt(0.95, df) * se

  TPM_equivalent <-
    FDR_equiv < 0.05 &
    CI90_lower > -delta &
    CI90_upper < delta

  res <- de_rs
  res$estimate       <- estimate[rownames(res)]
  res$se             <- se[rownames(res)]
  res$df             <- df[rownames(res)]
  res$p_equiv        <- p_equiv[rownames(res)]
  res$FDR_equiv      <- FDR_equiv[rownames(res)]
  res$CI90_lower     <- CI90_lower[rownames(res)]
  res$CI90_upper     <- CI90_upper[rownames(res)]
  res$equivalent     <- TPM_equivalent[rownames(res)]

  # single classification column: differentially expressed / equivalent / inconclusive
  res$classification <- ifelse(
    res$adj.P.Val < 0.05, "differentially_expressed",
    ifelse(res$equivalent, "equivalent", "inconclusive")
  )

  res <- res[order(res$adj.P.Val), ]
  res
}

delta_tpm <- 0.5

tpm_rs <- run_de_and_equivalence(gtex_tpm, design, coef_name, delta = delta_tpm)
div_rs <- run_de_and_equivalence(gtex_div, design, coef_name, delta = delta_tpm)

d$n_sig_genes_tpm <- sum(tpm_rs$classification == "differentially_expressed")
d$n_equiv_genes_tpm <- sum(tpm_rs$classification == "equivalent")
d$sig_genes_tpm  <- rownames(tpm_rs)[tpm_rs$classification == "differentially_expressed"]
d$equiv_genes_tpm <- rownames(tpm_rs)[tpm_rs$classification == "equivalent"]

d$n_sig_genes_div <- sum(div_rs$classification == "differentially_expressed")
d$n_equiv_genes_div <- sum(div_rs$classification == "equivalent")
d$sig_genes_div  <- rownames(div_rs)[div_rs$classification == "differentially_expressed"]
d$equiv_genes_div <- rownames(div_rs)[div_rs$classification == "equivalent"]

# ------------------------------ save results ------------------------------
fwrite(
  data.table(gene = rownames(tpm_rs), tpm_rs),
  paste0(opt$basedon,'_limma',"/",opt$pair,'.',d$t1,'.',d$t2,".tpm.limma_DE_equiv.csv.gz"),
  compress = "gzip"
)
fwrite(
  data.table(gene = rownames(div_rs), div_rs),
  paste0(opt$basedon,'_limma',"/",opt$pair,'.',d$t1,'.',d$t2,".div.limma_DE_equiv.csv.gz"),
  compress = "gzip"
)

cat(sprintf(
  "Pair %d: %s vs %s | TPM: %d DE, %d equivalent | DIV: %d DE, %d equivalent\n",
  opt$pair, d$t1, d$t2,
  d$n_sig_genes_tpm, d$n_equiv_genes_tpm,
  d$n_sig_genes_div, d$n_equiv_genes_div
))
cat("Analysis complete!\n")