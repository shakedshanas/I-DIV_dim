#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# isoform_age_mlr_limma.R
#
# Isoform-level multiple linear regression with limma, mirroring the
# gene-level script (limma_age_per_tissue.R): age enters as an INTEGER
# (age_int), together with SEX, SMRIN (RIN) and SMTSISCH (ischemic time).
#
# Same structure as isoform_age_de_limma.R: one --gene per run, all of its
# isoforms, looping over every SMTSD inside the chosen SMTS (default Brain).
#
# Model per isoform, per SMTSD:
#     y ~ age_int + SEX + SMRIN + SMTSISCH
#   - age_int = match(AGE_grouped, c("20-49","50-59","60-79")) -> 1, 2, 3
#     (identical to the gene-level script; ALL three age groups are used)
#   - covariates without variation inside a tissue are dropped automatically
#   - samples with missing covariates are dropped
#   - isoforms that cannot be tested stay in the output with NA statistics
#     and a 'status' column
#
# Output: one gzipped CSV per gene, one row per isoform x tissue. For every
# model term: <term>_coef, _SE, _t, _p, _FDR (BH within tissue and term).
#
# Usage:
#   Rscript isoform_age_mlr_limma.R --gene=ENSG00000000419
#   Rscript isoform_age_mlr_limma.R --gene=ENSG00000000419 --smts=Brain
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})


MIN_TISSUE_SAMPLES <- 10   # tissue needs MORE than this many samples
MIN_OBS <- 3               # isoform needs at least this many non-missing samples

COVARIATES <- c("SEX", "SMRIN", "SMTSISCH")

fwrite_gz <- function(df, path) {
  con <- gzfile(path, "w")
  on.exit(close(con))
  utils::write.csv(df, con, row.names = FALSE)
}
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# --------------------------- CLI ARGS ---------------------------------------
raw_args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  hit <- raw_args[grepl(paste0("^--", flag, "="), raw_args)]
  if (length(hit) == 0) return(NULL)
  sub(paste0("^--", flag, "="), "", hit[1])
}

gene_query <- get_arg("gene")
if (is.null(gene_query)) stop("--gene is required, e.g. --gene=ENSG00000000419")
smts <- get_arg("smts")
if (is.null(smts)) smts <- "Brain"

# --------------------------- PATHS -------------------------------------------
metadata_path <- "../../GTEx_upstream/metadata.csv.gz"
isoform_matrix_path <- "../../GTEx_upstream/QN.Brain.I.csv.gz"

outdir <- file.path(gene_query)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# --------------------------- LOAD METADATA -----------------------------------
metadata <- fread(cmd = paste("zcat", metadata_path), data.table = FALSE)
rownames(metadata) <- metadata$SAMPID
metadata$SAMPID <- NULL

required <- c("AGE", "SMTS", "SMTSD", "SEX", "SMRIN", "SMTSISCH")
missing_cols <- setdiff(required, colnames(metadata))
if (length(missing_cols) > 0)
  stop("metadata is missing column(s): ", paste(missing_cols, collapse = ", "),
       " -- check that this metadata file contains the covariates.")

metadata$age_int <- match(metadata[['AGE_int']], c(2,5,6))   # 1, 2, 3 (as in the gene-level script)


# covariate types (SEX coded 1/2 as in the gene-level script; labels kept if already text)
sex_raw <- as.character(metadata$SEX)
sex_raw[sex_raw %in% c("1", "1.0")] <- "Male"
sex_raw[sex_raw %in% c("2", "2.0")] <- "Female"
metadata$SEX      <- factor(sex_raw)
metadata$SMRIN    <- suppressWarnings(as.numeric(as.character(metadata$SMRIN)))
metadata$SMTSISCH <- suppressWarnings(as.numeric(as.character(metadata$SMTSISCH)))

# restrict to the chosen SMTS (all three age groups are kept)
metadata <- metadata[metadata$SMTS == smts, , drop = FALSE]

# drop samples with missing model variables
model_vars <- c("age_int", COVARIATES)
ok <- stats::complete.cases(metadata[, model_vars, drop = FALSE])
if (any(!ok)) message(sprintf("Dropping %d sample(s) with missing age/covariates.", sum(!ok)))
metadata <- metadata[ok, , drop = FALSE]

tissues <- names(sort(table(metadata$SMTSD), decreasing = TRUE))
message(sprintf("Found %d %s SMTSD tissues (post filters): %s",
                length(tissues), smts, paste(tissues, collapse = ", ")))

sample_ids_needed <- rownames(metadata)

# --------------------------- LOAD ISOFORM MATRIX -----------------------------
# 3 leading id columns: [1] gene_id (no version) [2] transcript_id 
header <- fread(file = isoform_matrix_path, nrows = 0, header = TRUE)
if (ncol(header) < 4) stop("isoform matrix does not have the expected >=4 columns")

id_cols <- colnames(header)[1:2]
transcript_col <- id_cols[1]
gene_col <- id_cols[2]


sample_cols_wanted <- intersect(colnames(header)[3:ncol(header)], sample_ids_needed)
if (length(sample_cols_wanted) == 0)
  stop("None of the selected SAMPIDs were found as columns in the isoform matrix -- check SAMPID naming.")

message(sprintf("Reading isoform matrix: %d id columns + %d sample columns",
                length(id_cols), length(sample_cols_wanted)))

expr_raw <- fread(
  file = isoform_matrix_path,
  select = c(id_cols, sample_cols_wanted),
  colClasses = "character",
  data.table = FALSE,
  header = TRUE
)

# defensive numeric parse (drop embedded header/artifact rows)
vals <- expr_raw[, sample_cols_wanted, drop = FALSE]
parsed <- suppressWarnings(lapply(vals, as.numeric))
parsed <- as.data.frame(parsed, stringsAsFactors = FALSE, check.names = FALSE)
is_bad_cell <- is.na(parsed) & !(is.na(vals) | vals == "")
row_is_artifact <- apply(is_bad_cell, 1, any)
if (any(row_is_artifact))
  message(sprintf("Dropping %d header-artifact/non-numeric row(s) from isoform matrix.", sum(row_is_artifact)))
expr_raw <- expr_raw[!row_is_artifact, , drop = FALSE]
parsed <- parsed[!row_is_artifact, , drop = FALSE]

# --------------------------- SUBSET TO THE QUERIED GENE -----------------------
gene_mask <- expr_raw[[gene_col]] == gene_query
if (sum(gene_mask) == 0)
  stop(sprintf("Gene '%s' not found in column '%s' of the isoform matrix.", gene_query, gene_col))

gene_mat <- as.matrix(parsed[gene_mask, , drop = FALSE])
rownames(gene_mat) <- expr_raw[[transcript_col]][gene_mask]
colnames(gene_mat) <- sample_cols_wanted

message(sprintf("Gene %s: %d transcripts, %d samples available.",
                gene_query, nrow(gene_mat), ncol(gene_mat)))

# # isoform matrix is raw TPM -> log2(TPM + 1)
# gene_mat <- log2(gene_mat + 1)

# --------------------------- PER-TISSUE LIMMA ---------------------------------
classify_rows <- function(x) {
  n_obs <- rowSums(!is.na(x))
  v <- apply(x, 1, function(r) if (sum(!is.na(r)) > 1) stats::var(r, na.rm = TRUE) else NA)
  ifelse(n_obs < MIN_OBS, "too_few_values",
         ifelse(is.na(v) | v == 0, "no_variance", "tested"))
}

run_one_tissue <- function(tissue) {
  meta_sub <- metadata[metadata$SMTSD == tissue & rownames(metadata) %in% colnames(gene_mat), , drop = FALSE]

  if (nrow(meta_sub) <= MIN_TISSUE_SAMPLES) {
    message(sprintf("  [skip] %s: only %d samples (need > %d)", tissue, nrow(meta_sub), MIN_TISSUE_SAMPLES))
    return(NULL)
  }
  if (length(unique(meta_sub$age_int)) < 2) {
    message(sprintf("  [skip] %s: only one age group present", tissue))
    return(NULL)
  }

  meta_sub <- droplevels(meta_sub)
  y_all <- gene_mat[, rownames(meta_sub), drop = FALSE]
  status <- classify_rows(y_all)
  names(status) <- rownames(y_all)
  if (!any(status == "tested")) {
    message(sprintf("  [skip] %s: no testable isoforms", tissue))
    return(NULL)
  }

  # raw joined dump (expression + age + covariates)
  raw_df <- cbind(as.data.frame(t(y_all)),
                  age_int = meta_sub$age_int, AGE = meta_sub[['AGE_int']],
                  meta_sub[, COVARIATES, drop = FALSE])
  raw_path <- file.path(outdir, paste(gene_query, tissue, "isoform", "raw", "csv", "gz", sep = "."))
  fwrite_gz(cbind(SAMPID = rownames(raw_df), raw_df), raw_path)

  # keep only covariates that vary in this tissue
  candidate <- c("age_int", COVARIATES)
  terms <- candidate[vapply(candidate, function(v) {
    x <- meta_sub[[v]]
    if (is.factor(x)) nlevels(x) >= 2 else length(unique(x)) > 1
  }, logical(1))]

  form <- as.formula(paste("~", paste(terms, collapse = " + ")))
  design <- model.matrix(form, data = meta_sub)

  ne <- nonEstimable(design)
  if (!is.null(ne)) {
    if ("age_int" %in% ne) {
      message(sprintf("  [skip] %s: age_int not estimable (confounded)", tissue))
      return(NULL)
    }
    message("  Dropping non-estimable column(s): ", paste(ne, collapse = ", "))
    design <- design[, !colnames(design) %in% ne, drop = FALSE]
  }
  if (nrow(design) - ncol(design) < 1) {
    message(sprintf("  [skip] %s: no residual degrees of freedom", tissue))
    return(NULL)
  }
  cat(sprintf("[%s] %d samples, formula %s, residual df %d\n",
              tissue, nrow(design), deparse(form), nrow(design) - ncol(design)))

  y <- y_all[status == "tested", , drop = FALSE]
  fit <- eBayes(lmFit(y, design))

  # one block per coefficient: coef, SE, t, p, FDR
  res <- data.frame(transcript_id = rownames(y), AveExpr = fit$Amean,
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

  # add back untested isoforms (NA statistics), keep input order
  untested <- names(status)[status != "tested"]
  if (length(untested) > 0) {
    pad <- res[rep(NA_integer_, length(untested)), , drop = FALSE]
    pad$transcript_id <- untested
    res <- rbind(res, pad)
  }
  res <- res[match(rownames(y_all), res$transcript_id), , drop = FALSE]
  res$status <- unname(status[res$transcript_id])
  res$gene_id <- gene_query
  res$SMTSD <- tissue
  res$n_samples <- nrow(design)
  res$formula <- deparse(form)
  rownames(res) <- NULL
  res
}

results_list <- lapply(tissues, function(tissue) {
  message(sprintf("Processing %s ...", tissue))
  run_one_tissue(tissue)
})
results_list <- Filter(Negate(is.null), results_list)

if (length(results_list) == 0)
  stop(sprintf("No %s SMTSD tissue had enough samples/age coverage for gene %s.", smts, gene_query))

# tissues can differ in which covariates survive -> fill missing columns with NA
final <- as.data.frame(rbindlist(results_list, fill = TRUE))
lead <- c("gene_id", "SMTSD", "transcript_id", "n_samples", "AveExpr", "status", "formula")
final <- final[, c(lead, setdiff(colnames(final), lead))]
if ("age_int_FDR" %in% colnames(final))
  final <- final[order(final$SMTSD, final$age_int_FDR), ]

out_path <- file.path(outdir, paste(gene_query, tolower(smts), "isoform", "limma", "MLR", "csv", "gz", sep = "."))
fwrite_gz(final, out_path)
message(sprintf("Wrote %s (%d rows across %d tissues).", out_path, nrow(final), length(results_list)))
