###### 
###### GMM K=4
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mclust)
})

setwd("/home/woodydrylab/FileShare/20260117_MultiLayer_Unifacial/two_cambium")

# ============================================================
# Settings
# ============================================================
BASE_DIR <- "/home/woodydrylab/FileShare/20260121_Xenium"
SUMMARY_FILE <- file.path(BASE_DIR, "WeiLi_UMAP_XeniumCorr_20260409_all_cambium_summary.tsv")

OUTDIR <- file.path("GMM_on_mean_median_correlation")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

K_USE <- 4
PLOT_WIDTH <- 6
PLOT_HEIGHT <- 6
PLOT_DPI <- 300

# ============================================================
# Helper
# ============================================================
normalize_cell_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\.", "-", x)
  x <- sub("-1$", "", x)
  x
}

zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

plot_cluster <- function(dt, cluster_col, title_txt, out_png) {
  p <- ggplot(dt, aes(UMAP_1, UMAP_2, color = .data[[cluster_col]])) +
    geom_point(alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("#4A6B82", "#C44233", "#E29E36", "#00B0F6")) +
    theme_classic() +
    labs(title = title_txt, x = "UMAP_1", y = "UMAP_2") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      aspect.ratio = 1
    )
  ggsave(out_png, p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
}

plot_heat <- function(dt, colname, title_txt, out_png) {
  p <- ggplot(dt, aes(UMAP_1, UMAP_2, color = .data[[colname]])) +
    geom_point(alpha = 0.8, size = 1.2) +
    scale_color_gradient(low = "grey90", high = "firebrick3") +
    theme_classic() +
    labs(title = title_txt, x = "UMAP_1", y = "UMAP_2") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      aspect.ratio = 1
    )
  ggsave(out_png, p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
}

# ============================================================
# 1. Read summary
# ============================================================
summary_dt <- fread(SUMMARY_FILE)
valid_jobs <- summary_dt[status == "ok"]

if (nrow(valid_jobs) == 0) {
  stop("No valid jobs in summary file.")
}

# ============================================================
# 2. Merge all samples
# ============================================================
all_list <- lapply(seq_len(nrow(valid_jobs)), function(i) {
  f <- valid_jobs$out_tsv[i]

  if (!file.exists(f)) return(NULL)

  dt <- fread(f)

  if (!all(c("barcode", "UMAP_1", "UMAP_2", "cor_spearman") %in% names(dt))) return(NULL)

  dt <- dt[!is.na(barcode) & !is.na(UMAP_1) & !is.na(UMAP_2) & !is.na(cor_spearman)]

  dt[, .(
    cell = normalize_cell_id(barcode),
    UMAP_1 = as.numeric(UMAP_1),
    UMAP_2 = as.numeric(UMAP_2),
    cor_spearman = as.numeric(cor_spearman)
  )]
})

all_dt <- rbindlist(all_list, use.names = TRUE, fill = TRUE)

if (nrow(all_dt) == 0) {
  stop("No usable rows collected from valid jobs.")
}

fwrite(all_dt, file.path(OUTDIR, "all_cells_long.tsv"), sep = "\t")

# ============================================================
# 3. per-cell summary
# ============================================================
cell_dt <- all_dt[, .(
  UMAP_1 = mean(UMAP_1),
  UMAP_2 = mean(UMAP_2),
  mean_s = mean(cor_spearman),
  median_s = median(cor_spearman),
  sd_s = sd(cor_spearman),
  n = .N
), by = cell]

cell_dt[is.na(sd_s), sd_s := 0]

# ============================================================
# 4. GMM on mean
# ============================================================
cat("Running GMM on mean correlation...\n")

cell_dt[, z_mean := zscore_safe(mean_s)]
cell_dt[, z_x := zscore_safe(UMAP_1)]
cell_dt[, z_y := zscore_safe(UMAP_2)]

mat_mean <- as.matrix(cell_dt[, .(z_x, z_y, z_mean)])
fit_mean <- Mclust(mat_mean, G = K_USE, verbose = FALSE)
cell_dt[, cluster_mean := as.character(fit_mean$classification)]

# ============================================================
# 5. GMM on median
# ============================================================
cat("Running GMM on median correlation...\n")

cell_dt[, z_median := zscore_safe(median_s)]
mat_median <- as.matrix(cell_dt[, .(z_x, z_y, z_median)])
fit_median <- Mclust(mat_median, G = K_USE, verbose = FALSE)
cell_dt[, cluster_median := as.character(fit_median$classification)]

# ============================================================
# 6. cell summary
# ============================================================
fwrite(cell_dt, file.path(OUTDIR, "cell_summary_with_clusters.tsv"), sep = "\t")


cell_dt_small <- cell_dt[, .(
  cell, UMAP_1, UMAP_2, mean_s, median_s, sd_s, n,
  cluster_mean, cluster_median
)]
fwrite(cell_dt_small, file.path(OUTDIR, "cell_summary_with_clusters_small.tsv"), sep = "\t")

# ============================================================
# 7. cluster summary
# ============================================================
summary_mean <- cell_dt[, .(
  n = .N,
  mean_x = mean(UMAP_1),
  mean_y = mean(UMAP_2),
  mean_corr = mean(mean_s)
), by = cluster_mean]

summary_median <- cell_dt[, .(
  n = .N,
  mean_x = mean(UMAP_1),
  mean_y = mean(UMAP_2),
  mean_corr = mean(median_s)
), by = cluster_median]

fwrite(summary_mean, file.path(OUTDIR, "cluster_summary_mean.tsv"), sep = "\t")
fwrite(summary_median, file.path(OUTDIR, "cluster_summary_median.tsv"), sep = "\t")

# ============================================================
# 8. plots
# ============================================================
plot_cluster(
  cell_dt, "cluster_mean",
  "GMM on mean Spearman",
  file.path(OUTDIR, "GMM_mean_clusters.png")
)

plot_cluster(
  cell_dt, "cluster_median",
  "GMM on median Spearman",
  file.path(OUTDIR, "GMM_median_clusters.png")
)

plot_heat(
  cell_dt, "mean_s",
  "Mean Spearman",
  file.path(OUTDIR, "UMAP_mean_heat.png")
)

plot_heat(
  cell_dt, "median_s",
  "Median Spearman",
  file.path(OUTDIR, "UMAP_median_heat.png")
)

cat("Done.\n")
cat("Output:", normalizePath(OUTDIR), "\n")