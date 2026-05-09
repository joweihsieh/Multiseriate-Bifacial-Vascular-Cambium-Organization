#!/usr/bin/env Rscript

setwd("/home/woodydrylab/FileShare/20260121_Xenium/Merged_DEG_k10_5domains_all")


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 0) paths
# ============================================================
ROOT_DIR <- "/home/woodydrylab/FileShare/20260121_Xenium/Merged_DEG_k10_5domains_all"
SYN_PATH <- "/home/woodydrylab/FileShare/20260121_Xenium/Ptrichocarpa_533_v4.1.synonym_20241227.txt"

# ============================================================
# helper: normalize Potri ID to Potri.
# ============================================================
norm_potri <- function(x) {
  x <- as.character(x)
  x <- gsub("^Potri[-_]", "Potri.", x)  # Potri- / Potri_ -> Potri.
  x
}

# ============================================================
# 1) load synonym table
# ============================================================
syn <- fread(SYN_PATH, sep = "\t", header = TRUE)

stopifnot(all(c("v4_1", "v3") %in% colnames(syn)))

syn[, v4_1 := as.character(v4_1)]
syn[, v3   := as.character(v3)]

# normalize v4.1 gene id
syn[, v4_1_norm := norm_potri(v4_1)]

# Remove duplicates and keep the first record
map_v3 <- unique(syn[, .(v4_1_norm, v3)], by = "v4_1_norm")

# ============================================================
# 2) function: add v3 annotation
# ============================================================
add_v3 <- function(dt, gene_col = "gene") {
  dt <- as.data.table(dt)
  stopifnot(gene_col %in% colnames(dt))

  dt[, gene_norm := norm_potri(get(gene_col))]

  out <- merge(
    dt,
    map_v3,
    by.x = "gene_norm",
    by.y = "v4_1_norm",
    all.x = TRUE
  )

  setnames(out, "v3", "gene_v3")

  # Standardize the gene column to the Potri. format
  out[, (gene_col) := gene_norm]
  out[, gene_norm := NULL]

  # Move gene and gene_v3 columns to the front
  front <- intersect(c("gene", "gene_v3"), names(out))
  setcolorder(out, c(front, setdiff(names(out), front)))

  out
}

# ============================================================
# 3) find DEG files
# ============================================================
files <- list.files(
  ROOT_DIR,
  pattern = "^DEG_all_5class_one_vs_rest(_sig)?\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Found", length(files), "target DEG files\n")
if (length(files) == 0) {
  stop("No DEG_all_5class_one_vs_rest.tsv or DEG_all_5class_one_vs_rest_sig.tsv found under: ", ROOT_DIR)
}

# ============================================================
# 4) process each file
# ============================================================
for (f in files) {
  cat("\nProcessing:", f, "\n")

  dt <- fread(f)

  if (!"gene" %in% colnames(dt)) {
    warning("Skip file (no 'gene' column): ", f)
    next
  }

  dt_v3 <- add_v3(dt, gene_col = "gene")

  out_f <- sub("\\.tsv$", "_withV3.tsv", f)
  fwrite(dt_v3, out_f, sep = "\t")

  cat("  Mapping rate:", sum(!is.na(dt_v3$gene_v3)), "/", nrow(dt_v3), "\n")
  cat("  Output:", out_f, "\n")
}

cat("\nDone.\n")


######################################################################
suppressPackageStartupMessages({
  library(data.table)
})

# -------- helper: normalize Potri ID to Potri. --------
norm_potri <- function(x) {
  x <- as.character(x)
  x <- gsub("^Potri[-_]", "Potri.", x)  # Potri- / Potri_ -> Potri.
  x
}

# =========================
# 0) Paths
# =========================
ROOT_DIR  <- "/home/woodydrylab/FileShare/20260121_Xenium/Merged_DEG_k10_5domains_all"
FUNC_PATH <- "/home/woodydrylab/FileShare/20260121_Xenium/Ptrichocarpa_v4.1_complete table_v1.1.csv"

# =========================
# 1) Load functional annotation table
# =========================
func <- fread(FUNC_PATH, sep = ",")

# Check required columns
stopifnot("Gene.ID" %in% colnames(func))

# Gene.ID should originally follow the Potri.XXX.v4.1 format
func[, Gene.ID := as.character(Gene.ID)]

# Keep only commonly used columns to avoid overly large output files
keep_func_cols <- intersect(
  c("Gene.ID","TF","TF.family_v4.1","TF.family_v3.0","Common.name","SCW",
    "Best.hit.Arabidopsis","Annotated.function",
    "Leaves_Ave","Xylem_Ave","Phloem_Ave","YoungShoot_Ave",
    "Xylem.specific.gene","X217.Xylem.specific.TF",
    "Xylem.high.expression.gene","X282.xylem.high.expression.TF",
    "Sap.peptide.precursor","Xylem.peptide.precursor"),
  colnames(func)
)
func <- func[, ..keep_func_cols]

# Remove duplicates
setkey(func, Gene.ID)
func <- unique(func, by = "Gene.ID")

# =========================
# 2) Function: attach functional annotation by gene
# =========================
add_func_anno <- function(dt, gene_col = "gene") {
  dt <- as.data.table(dt)
  stopifnot(gene_col %in% colnames(dt))

  dt[, gene_norm := norm_potri(get(gene_col))]
  dt[, Gene.ID := paste0(gene_norm, ".v4.1")]

  out <- merge(dt, func, by = "Gene.ID", all.x = TRUE)

  # Standardize the gene column to start with Potri.
  out[, (gene_col) := gene_norm]
  out[, gene_norm := NULL]

  # Move commonly used columns to the front
  front <- intersect(
    c("gene", "Gene.ID", "Best.hit.Arabidopsis", "Annotated.function",
      "Common.name", "TF", "TF.family_v4.1", "SCW"),
    colnames(out)
  )
  setcolorder(out, c(front, setdiff(names(out), front)))

  out
}

# =========================
# 3) Find DEG files
# =========================
files <- list.files(
  ROOT_DIR,
  pattern = "^DEG_all_5class_one_vs_rest(_sig)?\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Found", length(files), "target DEG files\n")
if (length(files) == 0) {
  stop("No DEG_all_5class_one_vs_rest.tsv or DEG_all_5class_one_vs_rest_sig.tsv found under: ", ROOT_DIR)
}

# =========================
# 4) Process each DEG file
# =========================
for (f in files) {
  cat("\nProcessing:", f, "\n")

  dt <- fread(f)

  if (!"gene" %in% colnames(dt)) {
    warning("Skip file (no 'gene' column): ", f)
    next
  }

  dt_withFunc <- add_func_anno(dt, gene_col = "gene")

  out_f <- sub("\\.tsv$", "_withFunc.tsv", f)
  fwrite(dt_withFunc, out_f, sep = "\t")

  n_mapped <- sum(!is.na(dt_withFunc$Annotated.function))
  cat("  Mapped Annotated.function:", n_mapped, "/", nrow(dt_withFunc), "\n")
  cat("  Output:", out_f, "\n")
}

cat("\nDone.\n")


#rsync -av --relative \
#"YCL:/home/woodydrylab/FileShare/20260121_Xenium/Merged_DEG_k10_5domains_all/./*/*/DEG_all_5class_one_vs_rest_sig_withFunc.tsv" \
#./

#rsync -av --relative \
#"YCL:/home/woodydrylab/FileShare/20260121_Xenium/Merged_DEG_k10_5domains_all/./*/*/DEG_all_5class_one_vs_rest_sig_withV3.tsv" \
#./

##### Use Excel to select cambium genes for ShinyGO 7.7 analysis


