#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# limma per tissue on gene expression (TPM) and isoform diversity (I-DIV),
# GTEx, run as a SLURM array job (one tissue per task).
#
# Model per gene:   y ~ age_int + SEX + SMRIN + SMTSISCH (+ SMTSD)
#   - age_int : first digit of the grouped age (20-49 -> 2, 50-59 -> 5,
#               60-79 -> 6), as requested
#   - SMTSD   : only enters the model when --basedon=SMTS and the broad tissue
#               contains several detailed tissues; within a single SMTSD it is
#               constant and dropped automatically
#   - any covariate without variation in the tissue (e.g. one sex) is dropped
#
# Output (one file per tissue and matrix), in ./<basedon>/:
#   <tissue>.limma.<matrix>.coefficients.csv.gz
#     one row per gene; for every model term: coef, SE, t, p, FDR (BH within
#     the tissue and term), plus AveExpr, n_samples and the formula used
#   <tissue>.limma.<matrix>.design.csv.gz
#     the design matrix actually used (for checking)
#
# Usage:   Rscript limma_age_per_tissue.R --tindex=0 --basedon=SMTSD
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

fwrite_gz <- function(df, path) {
  con <- gzfile(path, "w")
  on.exit(close(con))
  utils::write.csv(df, con, row.names = FALSE)
}

# --------------------------- SETTINGS ----------------------------------------
metadata_path <- "../../GTEx_upstream/metadata.csv.gz"
div_path <- "../../GTEx_upstream/QN.SMTSD.IDIV.csv.gz"
tpm_path <- "../../GTEx_upstream/QN.SMTSD.G.csv.gz"

MIN_SAMPLES <- 10        # tissues need MORE than this many samples

# --------------------------- CLI ARGS ----------------------------------------
raw_args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  hit <- raw_args[grepl(paste0("^--", flag, "="), raw_args)]
  if (length(hit) == 0) return(NULL)
  sub(paste0("^--", flag, "="), "", hit[1])
}

tindex_str <- get_arg("tindex")
basedon    <- get_arg("basedon")

if (is.null(tindex_str)) stop("--tindex is required, e.g. --tindex=0")
if (is.null(basedon))    stop("--basedon is required: SMTS or SMTSD")
if (!basedon %in% c("SMTS", "SMTSD")) stop("--basedon must be one of: SMTS, SMTSD")

tindex <- suppressWarnings(as.integer(tindex_str))
if (is.na(tindex)) stop("--tindex must be an integer, e.g. --tindex=0")

# --------------------------- LOAD METADATA -----------------------------------
metadata <- fread(cmd = paste("zcat", metadata_path), data.table = FALSE)
rownames(metadata) <- metadata$SAMPID
metadata$SAMPID <- NULL

required <- unique(c("SMTSD", "AGE", "SMTSISCH", "SMRIN", "SEX", basedon))
missing_cols <- setdiff(required, colnames(metadata))
if (length(missing_cols) > 0)
  stop("metadata is missing column(s): ", paste(missing_cols, collapse = ", "))

# grouped age and numeric age
metadata$age_int <- match(metadata[['AGE_int']], c(2,5,6))   # 1, 2, 3 (as in the gene-level script)


# covariate types
metadata$SEX      <- factor(metadata$SEX, levels = c(1, 2), labels = c("Male", "Female"))
metadata$SMRIN    <- suppressWarnings(as.numeric(as.character(metadata$SMRIN)))
metadata$SMTSISCH <- suppressWarnings(as.numeric(as.character(metadata$SMTSISCH)))

# drop samples with missing model variables
model_vars <- c("age_int", "SEX", "SMRIN", "SMTSISCH", "SMTSD")
ok <- stats::complete.cases(metadata[, model_vars, drop = FALSE])
if (any(!ok)) message(sprintf("Dropping %d sample(s) with missing covariates.", sum(!ok)))
metadata <- metadata[ok, , drop = FALSE]

# --------------------------- SELECT TISSUE -----------------------------------
tissue_counts <- sort(table(metadata[[basedon]]), decreasing = TRUE)
tissues_for_process <- names(tissue_counts[tissue_counts > MIN_SAMPLES])

if (tindex < 0 || tindex >= length(tissues_for_process)) {
  message(sprintf("tindex %d out of range [0, %d) for %d qualifying %s tissues -- exiting.",
                  tindex, length(tissues_for_process), length(tissues_for_process), basedon))
  quit(save = "no", status = 0)
}

t <- tissues_for_process[tindex + 1]  # R is 1-indexed
cat(tindex, t, "\n")

metadata <- metadata[metadata[[basedon]] == t, , drop = FALSE]
cat("nsamples:", nrow(metadata), "\n")
age_counts <- table(metadata[['AGE_int']])
cat("AGE group distribution:", paste(names(age_counts), age_counts, sep = "=", collapse = ", "), "\n")

sample_ids <- rownames(metadata)

# --------------------------- LOAD MATRICES (tissue samples only) -------------
read_matrix <- function(path) {
  header <- fread(cmd = paste("zcat", path), nrows = 0, header = TRUE)
  id_col <- colnames(header)[1]
  keep_cols <- intersect(colnames(header), c(id_col, sample_ids))
  df <- fread(cmd = paste("zcat", path), select = keep_cols, data.table = FALSE,
              colClasses = "character", header = TRUE)
  df <- df[!duplicated(df[[id_col]]), , drop = FALSE]
  rownames(df) <- df[[id_col]]
  df[[id_col]] <- NULL          # remove the ID column itself
  df
}

# character -> numeric matrix; drop embedded header/artifact rows
to_numeric_matrix <- function(df) {
  parsed <- suppressWarnings(as.data.frame(lapply(df, as.numeric), check.names = FALSE))
  bad <- is.na(parsed) & !(is.na(df) | df == "")
  row_bad <- apply(bad, 1, any)
  if (any(row_bad)) message(sprintf("Dropping %d non-numeric (header-artifact) row(s).", sum(row_bad)))
  mat <- as.matrix(parsed[!row_bad, , drop = FALSE])
  rownames(mat) <- rownames(df)[!row_bad]
  mat
}

mats <- list(
  tpm = to_numeric_matrix(read_matrix(tpm_path)),
  div = to_numeric_matrix(read_matrix(div_path))
)

# Genes are NOT removed from the output. Each matrix is handled on its own
# (no intersection between TPM and I-DIV), missing values are allowed (limma
# fits each gene on its available samples), and genes that cannot be tested
# are kept in the results table with NA statistics and a 'status' column:
#   tested          - model fitted
#   too_few_values  - fewer than MIN_OBS non-missing samples in this tissue
#   no_variance     - constant in this tissue (e.g. I-DIV = 0 in all samples,
#                     i.e. a single expressed isoform; or TPM = 0 everywhere)
MIN_OBS <- 3
gene_status <- list()
for (m in names(mats)) {
  x <- mats[[m]]
  x[is.infinite(x)] <- NA
  n_obs <- rowSums(!is.na(x))
  v <- apply(x, 1, function(r) if (sum(!is.na(r)) > 1) stats::var(r, na.rm = TRUE) else NA)
  status <- ifelse(n_obs < MIN_OBS, "too_few_values",
            ifelse(is.na(v) | v == 0, "no_variance", "tested"))
  names(status) <- rownames(x)
  cat(m, ":", nrow(x), "genes x", ncol(x), "samples; status:",
      paste(names(table(status)), table(status), sep = "=", collapse = ", "), "\n")
  gene_status[[m]] <- status
  mats[[m]] <- x
}

# --------------------------- LIMMA -------------------------------------------
out_dir <- basedon
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

run_limma <- function(matrix_name) {
  y_all  <- mats[[matrix_name]]
  status <- gene_status[[matrix_name]]
  samples <- intersect(colnames(y_all), rownames(metadata))
  y <- y_all[status == "tested", samples, drop = FALSE]
  meta <- droplevels(metadata[samples, , drop = FALSE])
  meta$SMTSD <- factor(meta$SMTSD)

  # keep only covariates that vary in this tissue
  candidate <- c("age_int", "SEX", "SMRIN", "SMTSISCH", "SMTSD")
  terms <- candidate[vapply(candidate, function(v) {
    x <- meta[[v]]
    if (is.factor(x)) nlevels(x) >= 2 else length(unique(x)) > 1
  }, logical(1))]
  if (!"age_int" %in% terms) stop(sprintf("Only one age group in tissue '%s'.", t))

  form <- as.formula(paste("~", paste(terms, collapse = " + ")))
  design <- model.matrix(form, data = meta)
  ne <- nonEstimable(design)
  if (!is.null(ne)) {
    if ("age_int" %in% ne) stop("age_int is not estimable (confounded with covariates).")
    message("Dropping non-estimable column(s): ", paste(ne, collapse = ", "))
    design <- design[, !colnames(design) %in% ne, drop = FALSE]
  }
  cat(sprintf("[%s] %s: %d samples, formula %s, residual df %d\n", matrix_name, t,
              ncol(y), deparse(form), nrow(design) - ncol(design)))

  fit <- eBayes(lmFit(y, design), trend = TRUE)

  # one block of columns per coefficient: coef, SE, t, p, FDR
  res <- data.frame(gene = rownames(y), tissue = t,
                    n_samples = rowSums(!is.na(y)),
                    formula = deparse(form), AveExpr = fit$Amean,
                    stringsAsFactors = FALSE)
  for (cf in colnames(design)) {
    tt <- topTable(fit, coef = cf, number = Inf, adjust.method = "BH", sort.by = "none")
    nm <- safe_name(cf)
    res[[paste0(nm, "_coef")]] <- tt$logFC
    res[[paste0(nm, "_SE")]]   <- fit$stdev.unscaled[, cf] * sqrt(fit$s2.post)
    res[[paste0(nm, "_t")]]    <- tt$t
    res[[paste0(nm, "_p")]]    <- tt$P.Value
    res[[paste0(nm, "_FDR")]]  <- tt$adj.P.Val
  }

  # add back the untested genes (NA statistics) so every input gene is listed
  untested <- names(status)[status != "tested"]
  if (length(untested) > 0) {
    pad <- res[rep(NA_integer_, length(untested)), , drop = FALSE]
    pad$gene <- untested
    pad$tissue <- t
    pad$n_samples <- rowSums(!is.na(y_all[untested, samples, drop = FALSE]))
    pad$formula <- deparse(form)
    res <- rbind(res, pad)
  }
  res$status <- unname(status[res$gene])
  res <- res[match(rownames(y_all), res$gene), , drop = FALSE]   # input order
  rownames(res) <- NULL
  cat(sprintf("[%s] %s: %d genes in output (%d tested)\n", matrix_name, t,
              nrow(res), sum(res$status == "tested")))

  prefix <- file.path(out_dir, paste(safe_name(t), "limma", matrix_name, sep = "."))
  fwrite_gz(res, paste0(prefix, ".coefficients.csv.gz"))
  fwrite_gz(cbind(SAMPID = rownames(design), as.data.frame(design)),
            paste0(prefix, ".design.csv.gz"))
  cat(sprintf("[%s] %s: %d genes with age_int FDR < 0.05\n", matrix_name, t,
              sum(res$age_int_FDR < 0.05, na.rm = TRUE)))
  invisible(res)
}

run_limma("tpm")
run_limma("div")