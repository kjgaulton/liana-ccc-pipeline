#!/usr/bin/env Rscript
# ==============================================================================
# 01_convert_seurat_to_h5ad.R
#
# LIANA+ (liana-py) is part of the scverse ecosystem and only reads
# AnnData/MuData objects — it cannot read a Seurat .rds directly. This script
# converts a Seurat single-cell RDS object into a .h5ad file that
# 02_run_liana_ccc_pipeline.py can consume.
#
# Requirements (install once):
#   install.packages(c("Seurat", "remotes"))
#   remotes::install_github("cellgeni/sceasy")
#   # sceasy calls reticulate under the hood, which needs a Python env with
#   # anndata installed:
#   #   reticulate::install_miniconda()          # if you have no Python env
#   #   reticulate::py_install("anndata")
#
# Usage:
#   Rscript 01_convert_seurat_to_h5ad.R <input.rds> <output.h5ad> [assay] [celltype_col]
#
# Example:
#   Rscript 01_convert_seurat_to_h5ad.R my_data.rds my_data.h5ad RNA cell_type
#
# Notes:
#   - LIANA+ expects log1p-normalized expression in AnnData's .X. This script
#     converts the Seurat "data" slot (Seurat's log-normalized layer) into .X.
#   - If your object stores raw counts only, it will be normalized with
#     Seurat's default NormalizeData() (log1p, CP10K) before conversion.
#   - `celltype_col` is optional here — it's just checked for existence — but
#     you MUST pass the same column name as --groupby to the Python script.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(sceasy)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 01_convert_seurat_to_h5ad.R <input.rds> <output.h5ad> [assay] [celltype_col]")
}

input_rds    <- args[1]
output_h5ad  <- args[2]
assay_name   <- if (length(args) >= 3) args[3] else "RNA"
celltype_col <- if (length(args) >= 4) args[4] else NA

cat(sprintf("[1/5] Reading Seurat object from: %s\n", input_rds))
seu <- readRDS(input_rds)

if (!inherits(seu, "Seurat")) {
  stop("The .rds file does not contain a Seurat object (class: ", paste(class(seu), collapse = ", "), ").")
}

cat(sprintf("[2/5] Using assay: %s\n", assay_name))
DefaultAssay(seu) <- assay_name

has_data_slot <- tryCatch({
  d <- GetAssayData(seu, assay = assay_name, layer = "data")
  !is.null(d) && max(d) > 0
}, error = function(e) FALSE)

if (!has_data_slot) {
  cat("[!] No log-normalized 'data' layer detected — running NormalizeData() (log1p, CP10K) first.\n")
  seu <- NormalizeData(seu, assay = assay_name, verbose = FALSE)
}

if (!is.na(celltype_col)) {
  if (!celltype_col %in% colnames(seu@meta.data)) {
    stop(sprintf(
      "Column '%s' not found in metadata. Available columns: %s",
      celltype_col, paste(colnames(seu@meta.data), collapse = ", ")
    ))
  }
  cat(sprintf("[3/5] Confirmed cell-type column exists: %s\n", celltype_col))
} else {
  cat("[3/5] No cell-type column specified for validation — make sure adata.obs has one before running the Python step.\n")
}

cat("[4/5] Converting Seurat -> AnnData with sceasy (this may take a while for large objects)...\n")
sceasy::convertFormat(
  seu,
  from       = "seurat",
  to         = "anndata",
  assay      = assay_name,
  main_layer = "data",   # log-normalized data becomes AnnData's .X, as LIANA+ expects
  outFile    = output_h5ad
)

cat(sprintf("[5/5] Done. Wrote: %s\n", output_h5ad))
