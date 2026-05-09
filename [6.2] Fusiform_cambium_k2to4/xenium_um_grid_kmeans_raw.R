#!/usr/bin/env Rscript
# xenium_grid_kmeans_raw.R
# Grid-based gene expression -> PCA -> K-means (NO normalization)
# Rscript xenium_5um_grid_kmeans_raw.R \
#  --matrix grid5um_out/counts_bins_by_genes_sparse.rds \
#  --binmeta grid5um_out/bin_metadata.tsv \
#  --outdir grid5um_out/kmeans_raw_out \
#  --k 4 \
#  --n_pcs 30


suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(irlba)
  library(optparse)
})

# ----------------------------
# CLI options
# ----------------------------
option_list <- list(
  make_option(c("-m","--matrix"), type="character",
              help="counts_bins_by_genes_sparse.rds"),
  make_option(c("-b","--binmeta"), type="character",
              help="bin_metadata.tsv"),
  make_option(c("-o","--outdir"), type="character", default="kmeans_raw_out"),
  make_option(c("--n_pcs"), type="integer", default=20,
              help="Number of PCs to use [default %default]"),
  make_option(c("--k"), type="integer", default=5,
              help="Number of K-means clusters [default %default]"),
  make_option(c("--top_genes"), type="integer", default=2000,
              help="Number of most variable genes to use [default %default]"),
  make_option(c("--seed"), type="integer", default=123)
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$matrix) || is.null(opt$binmeta)) {
  stop("Please provide --matrix and --binmeta")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

set.seed(opt$seed)

# ----------------------------
# Load data
# ----------------------------
cat("Loading matrix...\n")
mat <- readRDS(opt$matrix)   # bins x genes (dgCMatrix)

cat("Loading bin metadata...\n")
bin_meta <- fread(opt$binmeta)

stopifnot(nrow(bin_meta) == nrow(mat))

# ----------------------------
# Feature selection (raw counts)
# ----------------------------

# 每個 gene 出現在多少個 bins（非零）— 不會轉成 dense
detect_bins <- Matrix::colSums(mat > 0)

# 由大到小排序
ord <- order(detect_bins, decreasing = TRUE)

# 取前 top_genes 個（避免超過總 gene 數）
top_n <- min(opt$top_genes, length(ord))
top_genes <- colnames(mat)[ord[seq_len(top_n)]]

# 保險：移除可能的 NA（通常不會發生，但加了更穩）
top_genes <- top_genes[!is.na(top_genes)]

mat_use <- mat[, top_genes, drop = FALSE]

cat("Selected ", ncol(mat_use), " genes for PCA.\n", sep = "")

# ----------------------------
# PCA (NO normalization)
# ----------------------------
cat("Running PCA (raw counts)...\n")
pca <- prcomp_irlba(
  mat_use,
  n = opt$n_pcs,
  center = TRUE,
  scale. = FALSE
)

pcs <- pca$x   # bins x PCs

# Save PCA variance
pca_var <- pca$sdev^2 / sum(pca$sdev^2)
fwrite(
  data.table(PC = seq_along(pca_var), variance = pca_var),
  file = file.path(opt$outdir, "pca_variance.tsv"),
  sep = "\t"
)

# ----------------------------
# K-means clustering
# ----------------------------
cat("Running K-means clustering...\n")
km <- kmeans(
  pcs[, seq_len(opt$n_pcs), drop = FALSE],
  centers = opt$k,
  nstart = 50
)

cluster <- km$cluster

# ----------------------------
# Attach clusters to bins
# ----------------------------
bin_meta[, cluster_raw := factor(cluster)]

# ----------------------------
# Save outputs
# ----------------------------
saveRDS(
  list(
    pca = pca,
    kmeans = km,
    top_genes = top_genes
  ),
  file = file.path(opt$outdir, "kmeans_raw_results.rds")
)

fwrite(
  bin_meta,
  file = file.path(opt$outdir, "bin_metadata_with_cluster_raw.tsv"),
  sep = "\t"
)

cat("Done.\n")
cat("Outputs written to:", opt$outdir, "\n")
