#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# isoform_day0_day5_de_limma.R
#
# Adapted from isoform_age_de_limma.R for the neurodifferentiation timecourse
# dataset. Instead of GTEx metadata (SAMPID/AGE/SMTS/SMTSD), this dataset
# encodes the timepoint directly in each sample column name, e.g.
# "day0_5" = day 0, replicate 5. There is no tissue dimension to loop over,
# so this version runs a single limma DAY comparison (day0 vs day5) on one
# gene's isoforms at a time.
#
# Input file is expected to have the layout shown in your sample:
#   transcript_id,gene_id,day3_3,day2_1,day5_5,day1_4,day0_5, ...
# i.e. 2 leading id columns (transcript_id, gene_id -- gene_id may carry a
# version suffix like ENSG00000000419.12) followed by one column per sample,
# each named "day<N>_<replicate>". Values are assumed to be raw (non-logged)
# TPM, matching the isoform matrix convention in isoform_age_de_limma.R.
#
# Usage:
#   Rscript isoform_day0_day5_de_limma.R --gene=ENSG00000000419
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

DAY_LOW  <- "0"   # baseline
DAY_HIGH <- "5"   # comparison timepoint

# Skip the gene entirely if either day group ends up with fewer than this
# many replicates after column matching.
MIN_GROUP_SAMPLES <- 2

fwrite_gz <- function(df, path) {
  con <- gzfile(path, "w")
  on.exit(close(con))
  utils::write.csv(df, con, row.names = FALSE)
}

# strip a "gene_id" or "gene_id.version" style string down to the
# version-less accession, e.g. "ENSG00000000419.12" -> "ENSG00000000419"
strip_version <- function(x) sub("\\..*$", "", x)

# --------------------------- CLI ARGS ---------------------------------------
# Usage: Rscript isoform_day0_day5_de_limma.R --gene=ENSG00000000419
raw_args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  hit <- raw_args[grepl(paste0("^--", flag, "="), raw_args)]
  if (length(hit) == 0) return(NULL)
  sub(paste0("^--", flag, "="), "", hit[1])
}

gene_query <- get_arg("gene")
if (is.null(gene_query)) stop("--gene is required, e.g. --gene=ENSG00000000419")
gene_query_clean <- strip_version(gene_query)

# --------------------------- PATHS -------------------------------------------
# EDIT THIS to point at your day0-day5 isoform matrix file.
isoform_matrix_path <- "transcript_tpm.csv.gz"  # written by 1.combine_data.py

outdir <- file.path("isoform_day0_day5_de", gene_query)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# --------------------------- INSPECT HEADER, PICK DAY0/DAY5 COLUMNS ----------
header <- fread(file = isoform_matrix_path, nrows = 0, header = TRUE)
if (ncol(header) < 3) stop("isoform matrix does not have the expected >=3 columns")

id_cols <- c("transcript_id", "accession_number")
transcript_col <- id_cols[1]
gene_col <- id_cols[2]
if (transcript_col != "transcript_id" || gene_col != "gene_id") {
  message(sprintf("Note: using columns '%s' and '%s' as transcript_id/gene_id (based on position, not name).",
                   transcript_col, gene_col))
}

sample_cols_available <- setdiff(colnames(header), c("V1", "", id_cols))

# sample col name convention: day<N>_<replicate>
day_of_sample <- sub("^day(\\d+)_.*$", "\\1", sample_cols_available)
is_day_col <- grepl("^day\\d+_", sample_cols_available)
if (any(!is_day_col)) {
  message(sprintf("Note: %d sample column(s) don't match 'day<N>_<rep>' and will be ignored: %s",
                   sum(!is_day_col), paste(sample_cols_available[!is_day_col], collapse = ", ")))
}

day0_cols <- sample_cols_available[is_day_col & day_of_sample == DAY_LOW]
day5_cols <- sample_cols_available[is_day_col & day_of_sample == DAY_HIGH]

message(sprintf("Found %d day%s samples and %d day%s samples.",
                 length(day0_cols), DAY_LOW, length(day5_cols), DAY_HIGH))

if (length(day0_cols) < MIN_GROUP_SAMPLES || length(day5_cols) < MIN_GROUP_SAMPLES) {
  stop(sprintf("Not enough replicates for day%s vs day%s (need >= %d per group; found %d and %d).",
               DAY_LOW, DAY_HIGH, MIN_GROUP_SAMPLES, length(day0_cols), length(day5_cols)))
}

sample_cols_wanted <- c(day0_cols, day5_cols)

# --------------------------- LOAD MATRIX (id cols + day0/day5 cols only) -----
expr_raw <- fread(
  file = isoform_matrix_path,
  select = c(id_cols, sample_cols_wanted),
  colClasses = "character",
  data.table = FALSE,
  header = TRUE
)

# defensive numeric parse: drop any row where a "sample" cell doesn't parse
# as numeric (e.g. an embedded header/artifact row), same idea as in
# isoform_age_de_limma.R
vals <- expr_raw[, sample_cols_wanted, drop = FALSE]
parsed <- suppressWarnings(lapply(vals, as.numeric))
parsed <- as.data.frame(parsed, stringsAsFactors = FALSE, check.names = FALSE)
is_bad_cell <- is.na(parsed) & !(is.na(vals) | vals == "")
row_is_artifact <- apply(is_bad_cell, 1, any)
n_dropped <- sum(row_is_artifact)
if (n_dropped > 0) {
  message(sprintf("Dropping %d header-artifact/non-numeric row(s) from isoform matrix.", n_dropped))
}
expr_raw <- expr_raw[!row_is_artifact, , drop = FALSE]
parsed <- parsed[!row_is_artifact, , drop = FALSE]

# --------------------------- SUBSET TO THE QUERIED GENE -----------------------
gene_mask <- strip_version(expr_raw[[gene_col]]) == gene_query_clean
if (sum(gene_mask) == 0) {
  stop(sprintf("Gene '%s' not found in column '%s' of the isoform matrix.", gene_query, gene_col))
}

transcript_ids <- expr_raw[[transcript_col]][gene_mask]
gene_mat <- as.matrix(parsed[gene_mask, , drop = FALSE])
rownames(gene_mat) <- transcript_ids
colnames(gene_mat) <- sample_cols_wanted

message(sprintf("Gene %s: %d transcripts, %d day%s + %d day%s samples.",
                 gene_query, nrow(gene_mat), length(day0_cols), DAY_LOW, length(day5_cols), DAY_HIGH))

# log2(TPM + 1), same convention as isoform_age_de_limma.R
gene_mat <- log2(gene_mat + 1)

# raw joined dump for this gene, mirroring the age script's raw csv
raw_df <- data.frame(
  SAMPID = colnames(gene_mat),
  day = c(rep(DAY_LOW, length(day0_cols)), rep(DAY_HIGH, length(day5_cols))),
  t(gene_mat),
  check.names = FALSE
)
raw_path <- file.path(outdir, paste(gene_query, "isoform", "raw", "csv", "gz", sep = "."))
fwrite_gz(raw_df, raw_path)

# --------------------------- LIMMA: DAY0 vs DAY5 ------------------------------
condition <- factor(
  c(rep(DAY_LOW, length(day0_cols)), rep(DAY_HIGH, length(day5_cols))),
  levels = c(DAY_LOW, DAY_HIGH),
  labels = c("Day0", "Day5")
)
design <- model.matrix(~condition)
rownames(design) <- colnames(gene_mat)

fit <- lmFit(gene_mat, design)
fit <- eBayes(fit)

coef_name <- colnames(design)[2]

res <- topTable(
  fit,
  coef = coef_name,
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)
res$transcript_id <- rownames(res)
res$gene_id <- gene_query
res$n_Day0 <- length(day0_cols)
res$n_Day5 <- length(day5_cols)
rownames(res) <- NULL

final <- res[, c("gene_id", "transcript_id", "logFC", "AveExpr", "t",
                  "P.Value", "adj.P.Val", "B", "n_Day0", "n_Day5")]
final <- final[order(final$adj.P.Val), ]

out_path <- file.path(outdir, paste(gene_query, "day0", "day5", "isoform", "limma", "csv", "gz", sep = "."))
fwrite_gz(final, out_path)
message(sprintf("Wrote %s (%d transcripts).", out_path, nrow(final)))