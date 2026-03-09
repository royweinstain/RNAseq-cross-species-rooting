#!/usr/bin/env Rscript
# Step 1: Per-species differential expression analysis
# Reads config.yaml, runs DESeq2 (counts) or limma-voom (TPM/FPKM) for each species,
# writes standardized de_results.csv.
#
# Usage: Rscript steps/01_de_analysis.R config.yaml [--species argan]

suppressPackageStartupMessages({
  library(yaml)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(readxl)
})

# ── helpers ──────────────────────────────────────────────────────────────────

detect_format <- function(mat) {
  # mat: numeric matrix, genes x samples (no ID column)
  vals <- as.vector(mat)
  vals <- vals[!is.na(vals) & is.finite(vals) & vals >= 0]
  if (length(vals) == 0) stop("Expression matrix is empty or all-NA.")

  pct_integer <- mean(vals == floor(vals))
  max_val     <- max(vals)
  col_sums    <- colSums(mat, na.rm = TRUE)
  near_million <- mean(abs(col_sums - 1e6) / 1e6 < 0.05)  # within 5% of 1M

  if (pct_integer > 0.95 && max_val > 500) {
    return("counts")
  } else if (near_million > 0.80) {
    return("tpm")
  } else {
    # Ambiguous — check if values look like FPKM (not summing to 1M, not integers)
    if (pct_integer < 0.50) {
      return("fpkm")
    }
    stop(paste0(
      "Could not auto-detect data format.\n",
      "  % integer values: ", round(pct_integer * 100, 1), "%\n",
      "  Max value: ", round(max_val, 1), "\n",
      "  Columns near 1M sum: ", round(near_million * 100, 1), "%\n",
      "Please set data_format explicitly in config.yaml (counts/tpm/fpkm)."
    ))
  }
}


run_deseq2 <- function(counts_mat, easy_cols, hard_cols, species_name) {
  message("  Running DESeq2 for ", species_name)

  # Build colData
  condition <- factor(
    c(rep("easy", length(easy_cols)), rep("hard", length(hard_cols))),
    levels = c("easy", "hard")
  )
  col_data <- data.frame(condition = condition,
                         row.names = c(easy_cols, hard_cols))

  # Subset and reorder columns
  mat <- counts_mat[, c(easy_cols, hard_cols), drop = FALSE]
  storage.mode(mat) <- "integer"

  # Pre-filter: keep genes with >= 10 counts in >= min(n_easy, n_hard) samples
  min_samples <- min(length(easy_cols), length(hard_cols))
  keep <- rowSums(mat >= 10) >= min_samples
  message("  Genes before filter: ", nrow(mat), "  after filter: ", sum(keep))
  mat <- mat[keep, , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData   = col_data,
    design    = ~condition
  )
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds,
                 contrast       = c("condition", "hard", "easy"),
                 pAdjustMethod  = "BH")
  res_df <- as.data.frame(res)
  res_df$gene_id    <- rownames(res_df)
  res_df$species    <- species_name
  res_df$de_method  <- "DESeq2"

  # Standardize columns
  out <- data.frame(
    gene_id    = res_df$gene_id,
    log2FC     = res_df$log2FoldChange,
    pvalue     = res_df$pvalue,
    padj       = ifelse(is.na(res_df$padj), 1.0, res_df$padj),
    mean_expr  = res_df$baseMean,
    species    = species_name,
    de_method  = "DESeq2",
    stringsAsFactors = FALSE
  )
  return(out)
}


run_limma_log2 <- function(expr_mat, easy_cols, hard_cols, species_name) {
  # For data that is already log2-transformed (e.g. log2 FPKM).
  # Uses plain limma (lmFit + eBayes) — no voom, no re-transformation.
  message("  Running limma (log2 input, no voom) for ", species_name)

  mat <- expr_mat[, c(easy_cols, hard_cols), drop = FALSE]

  # Remove rows that are all-NA
  keep <- rowSums(is.na(mat)) < ncol(mat)
  mat  <- mat[keep, , drop = FALSE]
  message("  Genes after NA-filter: ", nrow(mat))

  condition <- factor(
    c(rep("easy", length(easy_cols)), rep("hard", length(hard_cols))),
    levels = c("easy", "hard")
  )
  design   <- model.matrix(~0 + condition)
  colnames(design) <- levels(condition)

  fit      <- lmFit(mat, design)
  contrast <- makeContrasts(hard - easy, levels = design)
  fit2     <- contrasts.fit(fit, contrast)
  fit2     <- eBayes(fit2)

  tt <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "none")

  out <- data.frame(
    gene_id    = rownames(tt),
    log2FC     = tt$logFC,
    pvalue     = tt$P.Value,
    padj       = ifelse(is.na(tt$adj.P.Val), 1.0, tt$adj.P.Val),
    mean_expr  = tt$AveExpr,
    species    = species_name,
    de_method  = "limma-log2",
    stringsAsFactors = FALSE
  )
  return(out)
}


run_limma_voom <- function(expr_mat, easy_cols, hard_cols, species_name) {
  message("  Running limma-voom for ", species_name)

  mat <- expr_mat[, c(easy_cols, hard_cols), drop = FALSE]

  # Remove rows that are all-zero or all-NA
  keep <- rowSums(is.na(mat)) < ncol(mat) & rowSums(mat, na.rm = TRUE) > 0
  mat  <- mat[keep, , drop = FALSE]
  message("  Genes after zero-filter: ", nrow(mat))

  # Replace any remaining NAs with 0 (DGEList does not accept NAs)
  mat[is.na(mat)] <- 0

  condition <- factor(
    c(rep("easy", length(easy_cols)), rep("hard", length(hard_cols))),
    levels = c("easy", "hard")
  )

  # TMM normalization via edgeR DGEList
  dge <- DGEList(counts = mat, group = condition)
  dge <- calcNormFactors(dge, method = "TMM")

  design   <- model.matrix(~0 + condition)
  colnames(design) <- levels(condition)

  v        <- voom(dge, design, plot = FALSE)
  fit      <- lmFit(v, design)
  contrast <- makeContrasts(hard - easy, levels = design)
  fit2     <- contrasts.fit(fit, contrast)
  fit2     <- eBayes(fit2)

  tt <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "none")

  out <- data.frame(
    gene_id    = rownames(tt),
    log2FC     = tt$logFC,
    pvalue     = tt$P.Value,
    padj       = ifelse(is.na(tt$adj.P.Val), 1.0, tt$adj.P.Val),
    mean_expr  = tt$AveExpr,
    species    = species_name,
    de_method  = "limma-voom",
    stringsAsFactors = FALSE
  )
  return(out)
}


process_species <- function(sp, config_dir) {
  name         <- sp$name
  expr_file    <- file.path(config_dir, sp$expression_file)
  id_col       <- sp$gene_id_column
  easy_cols    <- sp$samples$easy
  hard_cols    <- sp$samples$hard
  de_output    <- file.path(config_dir, sp$de_output)
  fmt          <- if (is.null(sp$data_format)) "auto" else sp$data_format

  message("\n[", name, "] Reading expression file: ", expr_file)
  if (!file.exists(expr_file)) stop("Expression file not found: ", expr_file)

  ext <- tolower(tools::file_ext(expr_file))
  if (ext %in% c("xlsx", "xls")) {
    sheet_name <- if (!is.null(sp$expression_sheet)) sp$expression_sheet else 1
    skip_rows  <- if (!is.null(sp$header_row)) sp$header_row - 1 else 0
    message("  Reading Excel sheet='", sheet_name, "', skip=", skip_rows, " rows")
    expr_raw <- as.data.frame(
      readxl::read_excel(expr_file, sheet = sheet_name, skip = skip_rows,
                         col_names = TRUE, .name_repair = "minimal")
    )
  } else {
    expr_raw <- read.csv(expr_file, row.names = NULL, check.names = FALSE,
                         stringsAsFactors = FALSE)
  }

  # Validate columns
  missing_id <- !id_col %in% colnames(expr_raw)
  if (missing_id) stop("gene_id_column '", id_col, "' not found in ", expr_file,
                       "\nAvailable columns: ", paste(colnames(expr_raw), collapse=", "))

  missing_easy <- setdiff(easy_cols, colnames(expr_raw))
  missing_hard <- setdiff(hard_cols, colnames(expr_raw))
  if (length(missing_easy) > 0)
    stop("Easy sample columns not found: ", paste(missing_easy, collapse=", "))
  if (length(missing_hard) > 0)
    stop("Hard sample columns not found: ", paste(missing_hard, collapse=", "))

  # Build numeric matrix
  gene_ids <- expr_raw[[id_col]]
  mat      <- as.matrix(expr_raw[, c(easy_cols, hard_cols), drop = FALSE])
  class(mat) <- "numeric"
  rownames(mat) <- gene_ids

  # Detect format
  if (fmt == "auto") {
    fmt <- detect_format(mat)
    message("  Auto-detected format: ", fmt)
  } else {
    message("  Using specified format: ", fmt)
  }

  # Run appropriate DE method
  if (fmt == "counts") {
    de_df <- run_deseq2(mat, easy_cols, hard_cols, name)
  } else if (fmt == "log2fpkm") {
    # Already log2-transformed: use plain limma (no voom — voom would re-transform)
    de_df <- run_limma_log2(mat, easy_cols, hard_cols, name)
  } else {
    # tpm or fpkm: use limma-voom
    de_df <- run_limma_voom(mat, easy_cols, hard_cols, name)
  }

  # Write output
  out_dir <- dirname(de_output)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  write.csv(de_df, de_output, row.names = FALSE, quote = FALSE)
  message("  Written: ", de_output, "  (", nrow(de_df), " genes)")
  invisible(de_df)
}


# ── main ──────────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript 01_de_analysis.R config.yaml [--species NAME]")

config_file <- args[1]
if (!file.exists(config_file)) stop("Config file not found: ", config_file)

# Determine config directory (all relative paths are resolved from there)
config_dir <- dirname(normalizePath(config_file))

config <- yaml.load_file(config_file)

# Optional --species filter
species_filter <- NULL
if ("--species" %in% args) {
  idx <- which(args == "--species")
  if (idx < length(args)) species_filter <- args[idx + 1]
}

species_list <- config$species
if (!is.null(species_filter)) {
  species_list <- Filter(function(s) s$name == species_filter, species_list)
  if (length(species_list) == 0)
    stop("Species '", species_filter, "' not found in config.")
  message("Running DE for species: ", species_filter)
} else {
  message("Running DE for all species: ",
          paste(sapply(species_list, `[[`, "name"), collapse=", "))
}

for (sp in species_list) {
  tryCatch(
    process_species(sp, config_dir),
    error = function(e) {
      message("\nERROR processing species '", sp$name, "': ", conditionMessage(e))
      quit(status = 1)
    }
  )
}

message("\nStep 1 complete.")
