# =====================================================================
# Build one CellChat object per condition.
#
# This script only runs the inference and saves the objects. Plotting is
# separate - anything that reads cellchat_<condition>.rds.
#
# Resume-safe: skips any condition already on disk, and caches the h5ad
# conversion, so an interrupted run can be restarted cheaply.
#
#   Rscript cellchat_run.R
# =====================================================================

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(CellChat)
  library(future)
})

# ---- Config ---------------------------------------------------------
# EDIT THESE. R cannot read config.py, so the paths are set here.
DATA_FILE <- "/path/to/nerandomilast_10celltypes.h5ad"
OUT_DIR   <- "/path/to/CellChat/Nerandomilast_10Celltypes"

IDENT_COL <- "cell_type"   # what CellChat groups cells by
GROUP_COL <- "type"        # what splits the object into conditions

MIN_CELLS          <- 10   # a cell type below this in a condition is dropped there
MAX_CELLS_PER_TYPE <- 800  # larger populations are downsampled to this
NBOOT              <- 20   # permutations; CellChat's default is 100
N_WORKERS          <- 20   # set to the core count of the machine

dir.create(OUT_DIR, showWarnings = FALSE)
set.seed(42)               # makes the downsampling reproducible
plan("multisession", workers = N_WORKERS)
options(future.globals.maxSize = 16 * 1024^3)

ts <- function() format(Sys.time())

# ---- Load -----------------------------------------------------------
# the h5ad -> SingleCellExperiment conversion is slow, so it is cached
sce_file <- file.path(OUT_DIR, "sce.rds")
if (file.exists(sce_file)) {
  cat(ts(), "sce.rds exists, skipping load.\n")
  sce <- readRDS(sce_file)
} else {
  cat(ts(), "Reading h5ad...\n")
  sce <- readH5AD(DATA_FILE, verbose = TRUE)
  saveRDS(sce, sce_file)
  cat(ts(), "Saved sce.rds\n")
}

cat("\nDimensions:", dim(sce), "\n")
cat("Conditions:\n"); print(table(colData(sce)[[GROUP_COL]]))
cat("Cell types:\n"); print(table(colData(sce)[[IDENT_COL]]))

meta <- as.data.frame(colData(sce))
expr_full <- assay(sce, "X")
rm(sce); gc()

# ---- Build one object per condition ---------------------------------
all_conds <- sort(unique(meta[[GROUP_COL]]))
cat("\n", ts(), "Conditions to process:", paste(all_conds, collapse = ", "), "\n\n")

for (COND in all_conds) {
  out_file <- file.path(OUT_DIR, sprintf("cellchat_%s.rds", COND))
  if (file.exists(out_file)) {
    cat(ts(), sprintf("SKIP %s (already done)\n\n", COND))
    next
  }
  cat(ts(), sprintf("===== START %s =====\n", COND))

  idx <- which(meta[[GROUP_COL]] == COND)
  cat(ts(), sprintf("Total cells: %d\n", length(idx)))

  expr_sub <- expr_full[, idx, drop = FALSE]
  meta_sub <- meta[idx, , drop = FALSE]
  meta_sub[[IDENT_COL]] <- as.character(meta_sub[[IDENT_COL]])

  # both filters are per condition: a type common overall can still be too
  # rare in one condition to give a usable probability
  ct_counts <- table(meta_sub[[IDENT_COL]])
  keep_cts  <- names(ct_counts)[ct_counts >= MIN_CELLS]

  subsample_idx <- unlist(lapply(keep_cts, function(ct) {
    ct_idx <- which(meta_sub[[IDENT_COL]] == ct)
    if (length(ct_idx) > MAX_CELLS_PER_TYPE) sample(ct_idx, MAX_CELLS_PER_TYPE)
    else ct_idx
  }))

  expr_sub <- expr_sub[, subsample_idx, drop = FALSE]
  meta_sub <- meta_sub[subsample_idx, , drop = FALSE]
  meta_sub[[IDENT_COL]] <- droplevels(factor(meta_sub[[IDENT_COL]]))

  # barcodes repeat across samples, and createCellChat needs unique names
  uniq_names <- make.unique(colnames(expr_sub))
  colnames(expr_sub) <- uniq_names
  rownames(meta_sub) <- uniq_names

  cat(ts(), sprintf("After subsampling: %d cells, %d cell types\n",
                    ncol(expr_sub), length(keep_cts)))
  print(table(meta_sub[[IDENT_COL]]))

  cc <- createCellChat(object = expr_sub, meta = meta_sub, group.by = IDENT_COL)
  cc@DB <- CellChatDB.mouse

  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)

  cat(ts(), sprintf("computeCommunProb (nboot=%d, %d workers)...\n", NBOOT, N_WORKERS))
  cc <- computeCommunProb(cc, type = "triMean", nboot = NBOOT)

  cc <- filterCommunication(cc, min.cells = MIN_CELLS)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")

  saveRDS(cc, out_file)
  cat(ts(), sprintf("===== DONE %s (%d pathways) =====\n\n",
                    COND, length(cc@netP$pathways)))
  rm(cc, expr_sub, meta_sub); gc()
}

cat(ts(), "Done. Objects in:", OUT_DIR, "\n")
