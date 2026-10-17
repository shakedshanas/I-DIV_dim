#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Joint Brain-region mixed-effects model (Brain SMTSD subregions only)
#
# Usage:
#   Rscript corr_ages_brain_lme_allcoefs.R --inferred-covs yes   # include InferredCov1/2
#   Rscript corr_ages_brain_lme_allcoefs.R --inferred-covs no    # exclude InferredCov1/2
#   (also accepted: --inferred-covs=yes, values yes/no/true/false/1/0)
#
# Fits, per gene, across ALL Brain SMTSD subregions at once:
#   y ~ SMTSD + SEX + age_int + SMRIN + SMTSISCH
#   random intercept: ~1 | donor
#
# For EVERY fixed-effect coefficient (intercept, each SMTSD level, SEX,
# age_int, SMRIN, SMTSISCH) the script
# extracts coefficient, std. error, df, t-value, p-value and FDR
# (Benjamini-Hochberg across genes, computed separately for each term).
#
# Outputs per matrix, in <group_label>/ (all carry the cov_tag):
#   <group_label>.lme.<matrix>.<tag>.csv.gz                 main table: age/sex columns
#                                                           (+ age_fdr, sex_fdr) PLUS
#                                                           coef/p_val/fdr for every term
#   <group_label>.lme.<matrix>.<tag>.allcoefs.long.csv.gz   one row per gene x term
#   <group_label>.lme.<matrix>.<tag>.allcoefs.wide.csv.gz   one row per gene; <stat>.<term>
#   summary_tables/, residuals/, raw table, and the .rds with all fitted models
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
})



fwrite_gz <- function(df, path) {
  con <- gzfile(path, "w")
  on.exit(close(con))
  utils::write.csv(df, con, row.names = FALSE)
}

# --------------------------- PATHS ------------------------------------------
metadata_path <- "../../GTEx_upstream/metadata.csv.gz"
div_path <- "../../GTEx_upstream/QN.Brain.IDIV.csv.gz"
tpm_path <- "../../GTEx_upstream/QN.Brain.G.csv.gz"

group_label <- "lme_across_brainSMTSD"  # used for output directory/filenames
if (!dir.exists(group_label)) dir.create(group_label, recursive = TRUE)

# --------------------------- LOAD METADATA -----------------------------------
metadata <- fread(cmd = paste("zcat", metadata_path), data.table = FALSE)
rownames(metadata) <- metadata$SAMPID

brain_pattern <- "^Brain"
metadata <- metadata[grepl(brain_pattern, metadata$SMTSD), , drop = FALSE]
if (nrow(metadata) == 0) {
  stop(sprintf("No samples with SMTSD matching '%s' were found in the metadata.",
               brain_pattern))
}
metadata$SMTSD <- droplevels(factor(metadata$SMTSD))


metadata$age_int <- match(metadata[['AGE_int']], c(2, 5, 6))
cat("nsamples (Brain only):", nrow(metadata), "\n")
cat("SMTSD levels:", paste(levels(metadata$SMTSD), collapse = ", "), "\n")

# --------------------------- MODEL SPECIFICATION -----------------------------
fixed_terms <- c("SMTSD", "SEX", "age_int", "SMRIN", "SMTSISCH")

cov_cols <- c(fixed_terms, "donor")
fixed_formula <- stats::as.formula(paste("y ~", paste(fixed_terms, collapse = " + ")))
cat("Model: ", deparse(fixed_formula), " | random = ~1|donor\n", sep = "")

missing_cov <- setdiff(cov_cols, colnames(metadata))
if (length(missing_cov) > 0) {
  stop(paste("Missing required metadata column(s):", paste(missing_cov, collapse = ", ")))
}

metadata$SEX <- factor(metadata$SEX, levels = c(1, 2), labels = c("Male", "Female"))

numeric_covs <- c("SMRIN", "SMTSISCH")
for (v in numeric_covs) {
  metadata[[v]] <- suppressWarnings(as.numeric(as.character(metadata[[v]])))
}

sample_ids <- rownames(metadata)

# --------------------------- LOAD EXPRESSION DATA -----------------------------
read_matrix <- function(path) {
  header <- fread(file = path, nrows = 0, header = TRUE)
  id_col <- colnames(header)[1]
  keep_cols <- intersect(colnames(header), c(id_col, sample_ids))
  df <- fread(cmd = paste("zcat", path), select = keep_cols, data.table = FALSE,
              colClasses = "character", header = TRUE)
  df <- df[!duplicated(df[[id_col]]), , drop = FALSE]
  rownames(df) <- df[[id_col]]
  df$id <- NULL

  df
}

clean_matrix <- function(df) {
  id_col <- colnames(df)[1]
  sample_cols <- setdiff(colnames(df), id_col)

  vals <- df[, sample_cols, drop = FALSE]
  parsed <- suppressWarnings(lapply(vals, as.numeric))
  parsed <- as.data.frame(parsed, stringsAsFactors = FALSE, check.names = FALSE)

  is_bad_cell <- is.na(parsed) & !(is.na(vals) | vals == "")
  row_is_artifact <- apply(is_bad_cell, 1, any)

  n_dropped <- sum(row_is_artifact)
  if (n_dropped > 0) {
    message(sprintf("clean_matrix: dropping %d header-artifact row(s)", n_dropped))
  }
  df <- df[!row_is_artifact, , drop = FALSE]
  parsed <- parsed[!row_is_artifact, , drop = FALSE]

  mat <- as.matrix(parsed)
  rownames(mat) <- df[[id_col]]
  mat
}

gtex_list <- list(
  div = read_matrix(div_path),
  tpm = read_matrix(tpm_path)
)

for (m in names(gtex_list)) {
  gtex_list[[m]] <- clean_matrix(gtex_list[[m]])
}

for (metric in c("tpm", "div")) {
  mat <- gtex_list[[metric]]
  mat[is.infinite(mat)] <- NA
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  cat(metric, "genes remaining after Inf/NA drop:", nrow(mat), "\n")
  gtex_list[[metric]] <- t(mat)  # rows = samples, columns = genes
}

common_genes <- intersect(colnames(gtex_list[["div"]]), colnames(gtex_list[["tpm"]]))
cat("common genes (div n tpm):", length(common_genes), "\n")

if (length(common_genes) == 0) {
  stop(sprintf(
    paste(
      "No common genes between tpm (%d genes) and div (%d genes).",
      "This usually means one of the input files failed to parse correctly",
      "(e.g. a header-detection issue, or the two files use non-overlapping",
      "gene ID sets). Check the per-file 'genes remaining after Inf/NA drop'",
      "counts above and inspect the raw .csv.gz files."
    ),
    ncol(gtex_list[["tpm"]]), ncol(gtex_list[["div"]])
  ))
}

gtex_list[["div"]] <- gtex_list[["div"]][, common_genes, drop = FALSE]
gtex_list[["tpm"]] <- gtex_list[["tpm"]][, common_genes, drop = FALSE]

# --------------------------- MODEL FITTING ------------------------------------
fit_ages_all <- function(matrix_name) {
  expr_mat <- gtex_list[[matrix_name]]

  common_samples <- intersect(rownames(expr_mat), rownames(metadata))
  expr_mat <- expr_mat[common_samples, , drop = FALSE]
  covs <- metadata[common_samples, cov_cols, drop = FALSE]

  tempdf <- cbind(as.data.frame(expr_mat), covs)

  raw_path <- file.path(
    group_label,
    paste(group_label, "lme", matrix_name,  "raw", "csv", "gz", sep = ".")
  )
  fwrite_gz(cbind(SAMPID = rownames(tempdf), tempdf), raw_path)

  genes <- colnames(expr_mat)
  n_genes <- length(genes)

  corr_df <- data.frame(
    gene        = genes,
    age_coef    = NA_real_,
    age_p_val   = NA_real_,
    sex_coef    = NA_real_,
    sex_p_val   = NA_real_,
    row.names = genes,
    stringsAsFactors = FALSE
  )

  summary_dir <- file.path(group_label, "summary_tables")
  if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)

  resid_dir <- file.path(group_label, "residuals")
  if (!dir.exists(resid_dir)) dir.create(resid_dir, recursive = TRUE)

  fit_list <- vector("list", n_genes)
  names(fit_list) <- genes

  # NEW: every fixed-effect row of every gene's fit, in long format
  coef_list <- vector("list", n_genes)
  names(coef_list) <- genes

  pb <- txtProgressBar(min = 0, max = n_genes, style = 3)
  for (i in seq_along(genes)) {
    gene <- genes[i]

    dat <- cbind(
      data.frame(SAMPID = rownames(tempdf), y = tempdf[[gene]],
                 stringsAsFactors = FALSE),
      tempdf[, cov_cols, drop = FALSE]
    )
    rownames(dat) <- NULL
    dat <- dat[stats::complete.cases(dat), , drop = FALSE]
    dat <- droplevels(dat)

    fit <- tryCatch(
      nlme::lme(
        fixed  = fixed_formula,
        random = ~ 1 | donor,
        data   = dat,
        method = "REML"
      ),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      s <- summary(fit)$tTable

      if ("age_int" %in% rownames(s)) {
        corr_df[gene, "age_coef"]  <- s["age_int", "Value"]
        corr_df[gene, "age_p_val"] <- s["age_int", "p-value"]
      }

      sex_row <- grep("^SEX", rownames(s), value = TRUE)
      if (length(sex_row) == 1) {
        corr_df[gene, "sex_coef"]  <- s[sex_row, "Value"]
        corr_df[gene, "sex_p_val"] <- s[sex_row, "p-value"]
      }

      # NEW: keep all coefficients of this gene
      coef_list[[gene]] <- data.frame(
        gene   = gene,
        term   = rownames(s),
        coef   = s[, "Value"],
        se     = s[, "Std.Error"],
        df     = s[, "DF"],
        t_val  = s[, "t-value"],
        p_val  = s[, "p-value"],
        stringsAsFactors = FALSE,
        row.names = NULL
      )

      summary_out_path <- file.path(
        summary_dir,
        paste(group_label, "lme", gene, matrix_name,  "csv", "gz", sep = ".")
      )
      fwrite_gz(cbind(term = rownames(s), as.data.frame(s)), summary_out_path)

      resid_df <- data.frame(
        SAMPID     = dat$SAMPID,
        donor      = dat$donor,
        SMTSD      = dat$SMTSD,
        cond_resid = resid(fit),
        marg_resid = resid(fit, level = 0)
      )
      resid_out_path <- file.path(
        resid_dir,
        paste(group_label, "lme", gene, matrix_name,  "resid", "csv", "gz", sep = ".")
      )
      fwrite_gz(resid_df, resid_out_path)

      fit_list[[gene]] <- fit
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)

  # ---- NEW: FDR (BH across genes) separately for every coefficient ----------
  coef_long <- data.table::rbindlist(coef_list, use.names = TRUE)
  coef_long[, fdr := p.adjust(p_val, method = "BH"), by = term]

  long_path <- file.path(
    group_label,
    paste(group_label, "lme", matrix_name,  "allcoefs", "long", "csv", "gz", sep = ".")
  )
  fwrite_gz(as.data.frame(coef_long), long_path)

  # wide version: one row per gene, columns "<stat>.<term>", grouped by term
  coef_wide <- data.table::dcast(
    coef_long, gene ~ term,
    value.var = c("coef", "p_val", "fdr"),
    sep = "."
  )
  terms_all <- unique(coef_long$term)
  ordered_cols <- c("gene", unlist(lapply(terms_all, function(tm) {
    paste(c("coef", "p_val", "fdr"), tm, sep = ".")
  })))
  ordered_cols <- intersect(ordered_cols, colnames(coef_wide))
  coef_wide <- coef_wide[, ..ordered_cols]

  wide_path <- file.path(
    group_label,
    paste(group_label, "lme", matrix_name,  "allcoefs", "wide", "csv", "gz", sep = ".")
  )
  fwrite_gz(as.data.frame(coef_wide), wide_path)

  # FDR columns for the age / sex summary
  corr_df$age_fdr <- p.adjust(corr_df$age_p_val, method = "BH")
  corr_df$sex_fdr <- p.adjust(corr_df$sex_p_val, method = "BH")

  # main per-matrix table: age/sex columns PLUS coef / p_val / fdr for every term
  corr_df_all <- merge(corr_df, as.data.frame(coef_wide), by = "gene",
                       all.x = TRUE, sort = FALSE)

  out_path <- file.path(
    group_label,
    paste(group_label, "lme", matrix_name,  "csv", "gz", sep = ".")
  )
  fwrite_gz(corr_df_all, out_path)

  saveRDS(
    fit_list,
    file.path(group_label,
              paste(group_label, "lme", matrix_name,  "allmodels", "rds", sep = "."))
  )
}

fit_ages_all("tpm")
fit_ages_all("div")