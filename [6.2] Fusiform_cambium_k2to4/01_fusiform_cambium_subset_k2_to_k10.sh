#!/bin/bash
conda activate r_spatial

set -euo pipefail

# ============================================================
# user settings
# ============================================================
BASE_ROOT="/home/woodydrylab/FileShare/20260121_Xenium"
EXCEL_FILE="/home/woodydrylab/FileShare/20260121_Xenium/Cambium_fusiform.xlsx"
KMEANS_SCRIPT="/home/woodydrylab/FileShare/20260121_Xenium/xenium_um_grid_kmeans_raw.R"

K_START=02
K_END=10

TASK_TSV="${BASE_ROOT}/tmp_cambium_fusiform_tasks.tsv"

export BASE_ROOT EXCEL_FILE KMEANS_SCRIPT K_START K_END TASK_TSV

# ============================================================
# Step 0. read Excel and build task table
# ============================================================
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

base_root  <- Sys.getenv("BASE_ROOT")
excel_file <- Sys.getenv("EXCEL_FILE")
task_tsv   <- Sys.getenv("TASK_TSV")

if (!file.exists(excel_file)) {
  stop("Excel file not found: ", excel_file)
}

dt <- as.data.table(read_excel(excel_file))

required_cols <- c("sample_id", "label", "um", "k", "cambium_fusiform")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Excel file is missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Format grid directory name
dt[, grid_dir := fifelse(
  grepl("^grid[0-9]+um_out$", as.character(um)),
  as.character(um),
  paste0("grid", sprintf("%02d", as.integer(gsub("[^0-9]", "", as.character(um)))), "um_out")
)]

# k
dt[, k0 := as.integer(k)]

# target cluster
dt[, select_cluster := as.character(as.integer(cambium_fusiform))]
dt[, select_cluster_fmt := sprintf("%02d", as.integer(cambium_fusiform))]

# full paths
dt[, sample_dir := file.path(base_root, sample_id)]
dt[, src_dir := file.path(sample_dir, grid_dir, paste0("kmeans_k", k0, "_raw_out"))]
dt[, subdir_name := paste0("subset_from_k", k0, "_c", select_cluster_fmt)]
dt[, subdir := file.path(sample_dir, grid_dir, subdir_name)]

# keep only needed columns
out <- dt[, .(
  sample_id,
  label,
  grid_dir,
  k0,
  select_cluster,
  select_cluster_fmt,
  sample_dir,
  src_dir,
  subdir_name,
  subdir
)]

fwrite(out, task_tsv, sep = "\t")
cat("Task table written to:", task_tsv, "\n")
EOF

# ============================================================
# Step 1 + Step 2. loop all tasks
# ============================================================
while IFS=$'\t' read -r sample_id label grid_dir k0 select_cluster select_cluster_fmt sample_dir src_dir subdir_name subdir
do
  if [[ "$sample_id" == "sample_id" ]]; then
    continue
  fi

  echo "================================================="
  echo "Sample ID       : ${sample_id}"
  echo "Label           : ${label}"
  echo "Grid dir        : ${grid_dir}"
  echo "K0              : ${k0}"
  echo "Select cluster  : ${select_cluster}"
  echo "Subset dir      : ${subdir}"
  echo "================================================="

  mkdir -p "${subdir}"

  export SAMPLE_DIR="${sample_dir}"
  export GRID_DIR="${grid_dir}"
  export SRC_DIR="${src_dir}"
  export SUB_DIR="${subdir}"
  export SELECT_CLUSTER="${select_cluster}"

  # ----------------------------------------------------------
  # Step 1. subset one cluster only
  # ----------------------------------------------------------
  Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

sample_dir      <- Sys.getenv("SAMPLE_DIR")
grid_dir        <- Sys.getenv("GRID_DIR")
src_dir         <- Sys.getenv("SRC_DIR")
subdir          <- Sys.getenv("SUB_DIR")
select_cluster  <- Sys.getenv("SELECT_CLUSTER")

mat_file  <- file.path(sample_dir, grid_dir, "counts_bins_by_genes_sparse.rds")
meta_file <- file.path(sample_dir, grid_dir, "bin_metadata.tsv")
clu_file  <- file.path(src_dir, "bin_metadata_with_cluster_raw.tsv")

if (!file.exists(mat_file))  stop("Missing matrix: ", mat_file)
if (!file.exists(meta_file)) stop("Missing meta: ", meta_file)
if (!file.exists(clu_file))  stop("Missing cluster file: ", clu_file)

mat  <- readRDS(mat_file)
meta <- fread(meta_file)
clu  <- fread(clu_file)

stopifnot("bin_id" %in% names(meta), "bin_id" %in% names(clu))

meta2 <- merge(meta, clu[, .(bin_id, cluster_raw)], by = "bin_id")
meta2[, cluster_raw := as.character(cluster_raw)]

keep_ids <- meta2[cluster_raw == select_cluster, bin_id]
cat("Selected bins for cluster", select_cluster, ":", length(keep_ids), "\n")

if (length(keep_ids) == 0) {
  stop("No bins found for selected cluster: ", select_cluster)
}

stopifnot(!is.null(rownames(mat)))
mat_sub  <- mat[keep_ids, , drop = FALSE]
meta_sub <- meta2[bin_id %in% keep_ids]

dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
saveRDS(mat_sub, file.path(subdir, "counts_subset.rds"))
fwrite(meta_sub, file.path(subdir, "binmeta_subset.tsv"), sep = "\t")

cat("Subset written to:", subdir, "\n")
EOF

  # ----------------------------------------------------------
  # Step 2. run subset kmeans k=2..10 + plot
  # ----------------------------------------------------------
  for k in $(seq -w ${K_START} ${K_END})
  do
    outdir="${subdir}/kmeans_subset_k${k}_raw_out"
    mkdir -p "${outdir}"

    echo "Running subset k=${k} for ${sample_id}"

    Rscript "${KMEANS_SCRIPT}" \
      --matrix "${subdir}/counts_subset.rds" \
      --binmeta "${subdir}/binmeta_subset.tsv" \
      --outdir "${outdir}" \
      --k "${k}" \
      --n_pcs 30

    export OUTDIR="${outdir}"
    export K_NOW="${k}"

    Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

outdir <- Sys.getenv("OUTDIR")
k <- as.integer(Sys.getenv("K_NOW"))

f <- file.path(outdir, "bin_metadata_with_cluster_raw.tsv")
if (!file.exists(f)) stop("Missing: ", f)

bin_meta <- fread(f)

required_cols <- c("x_center", "y_center", "cluster_raw")
missing_cols <- setdiff(required_cols, names(bin_meta))
if (length(missing_cols) > 0) {
  stop("bin_metadata_with_cluster_raw.tsv is missing required columns: ", paste(missing_cols, collapse = ", "))
}

bin_meta <- bin_meta[!is.na(cluster_raw)]
bin_meta[, cluster_raw := factor(as.character(cluster_raw))]

clusters <- levels(bin_meta[["cluster_raw"]])
if (length(clusters) == 0) {
  cat("No cluster levels for k=", k, " at ", outdir, "\n", sep = "")
  quit(save = "no", status = 0)
}

# infer bin size from coordinates
infer_bin_size <- function(v) {
  u <- sort(unique(v))
  d <- diff(u)
  d <- d[d > 0]
  if (length(d) == 0) return(1)
  median(d)
}

bin_um_x <- infer_bin_size(bin_meta$x_center)
bin_um_y <- infer_bin_size(bin_meta$y_center)
bin_um <- median(c(bin_um_x, bin_um_y))

get_plot_limits <- function(dt, bin_um, x_col = "x_center", y_col = "y_center") {
  list(
    xlim = c(
      min(dt[[x_col]], na.rm = TRUE) - bin_um / 2,
      max(dt[[x_col]], na.rm = TRUE) + bin_um / 2
    ),
    ylim = c(
      min(dt[[y_col]], na.rm = TRUE) - bin_um / 2,
      max(dt[[y_col]], na.rm = TRUE) + bin_um / 2
    )
  )
}

plot_lim <- get_plot_limits(bin_meta, bin_um = bin_um)

black_theme <- function() {
  theme_void() +
    theme(
      plot.background   = element_rect(fill = "black", color = NA),
      panel.background  = element_rect(fill = "black", color = NA),
      legend.background = element_rect(fill = "black", color = NA),
      legend.key        = element_rect(fill = "black", color = NA),
      legend.text       = element_text(color = "white"),
      legend.title      = element_text(color = "white"),
      plot.title        = element_text(color = "white", hjust = 0.5),
      plot.margin       = margin(0, 0, 0, 0)
    )
}

pal <- hue_pal()(length(clusters))
names(pal) <- clusters

k_fmt <- sprintf("%02d", k)

# 1) all clusters
p_all <- ggplot(bin_meta, aes(x = x_center, y = y_center, fill = cluster_raw)) +
  geom_tile(width = bin_um, height = bin_um) +
  coord_equal(
    xlim = plot_lim$xlim,
    ylim = rev(plot_lim$ylim),
    expand = FALSE
  ) +
  scale_fill_manual(values = pal, drop = FALSE) +
  black_theme() +
  labs(title = paste0("Subset k=", k), fill = "Cluster")

ggsave(
  file.path(outdir, paste0("subset_k", k_fmt, "_all.png")),
  p_all, width = 6, height = 6, dpi = 300, bg = "black"
)

ggsave(
  file.path(outdir, paste0("subset_k", k_fmt, "_all.pdf")),
  p_all, width = 6, height = 6, bg = "black"
)

# 2) each cluster
for (cl in clusters) {
  dt_cl <- bin_meta[cluster_raw == cl]
  cl_fmt <- sprintf("%02d", as.integer(cl))

  p_one <- ggplot() +
    geom_tile(
      data = bin_meta,
      aes(x = x_center, y = y_center),
      width = bin_um, height = bin_um,
      fill = "black"
    ) +
    geom_tile(
      data = dt_cl,
      aes(x = x_center, y = y_center),
      width = bin_um, height = bin_um,
      fill = pal[[cl]]
    ) +
    coord_equal(
      xlim = plot_lim$xlim,
      ylim = rev(plot_lim$ylim),
      expand = FALSE
    ) +
    black_theme() +
    labs(title = paste0("Subset k=", k, " cluster=", cl))

  ggsave(
    file.path(outdir, paste0("subset_k", k_fmt, "_cluster_", cl_fmt, ".png")),
    p_one, width = 6, height = 6, dpi = 300, bg = "black"
  )

  ggsave(
    file.path(outdir, paste0("subset_k", k_fmt, "_cluster_", cl_fmt, ".pdf")),
    p_one, width = 6, height = 6, bg = "black"
  )
}
EOF

  done

done < "${TASK_TSV}"

echo "================================================="
echo "All samples finished."
echo "================================================="