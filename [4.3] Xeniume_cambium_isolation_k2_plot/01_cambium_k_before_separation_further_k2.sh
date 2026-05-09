#!/bin/bash
conda activate r_spatial

set -euo pipefail

# ============================================================
# user settings
# ============================================================
BASE_ROOT="/home/woodydrylab/FileShare/20260121_Xenium"
EXCEL_FILE="/home/woodydrylab/FileShare/20260121_Xenium/Cambium_k_single_furtherk2.xlsx"
KMEANS_SCRIPT="/home/woodydrylab/FileShare/20260121_Xenium/xenium_um_grid_kmeans_raw.R"


TASK_TSV="${BASE_ROOT}/tmp_cambium_single_tasks.tsv"

export BASE_ROOT EXCEL_FILE KMEANS_SCRIPT TASK_TSV

rm -f "${TASK_TSV}"

# ============================================================
# Step 0. read Excel and build task table
# 雙樣品 sample_id 一律展開成 stem1 / stem2
# 單樣品保留 single
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
  stop("找不到 Excel 檔案: ", excel_file)
}

dt <- as.data.table(read_excel(excel_file))

required_cols <- c("sample_id", "label", "um", "k", "cambium_single")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Excel 缺少必要欄位: ", paste(missing_cols, collapse = ", "))
}

dt[, grid_dir := fifelse(
  grepl("^grid[0-9]+um_out$", as.character(um)),
  as.character(um),
  paste0("grid", sprintf("%02d", as.integer(gsub("[^0-9]", "", as.character(um)))), "um_out")
)]

dt[, source_k := as.integer(k)]
dt[, source_k_fmt := sprintf("%02d", source_k)]
dt[, select_cluster := as.character(as.integer(cambium_single))]
dt[, select_cluster_fmt := sprintf("%02d", as.integer(cambium_single))]
dt[, sample_dir := file.path(base_root, sample_id)]

# 雙樣品資料夾判斷
dt[, is_double := grepl("TISSUE_[0-9]+_and_[0-9]+", sample_id)]

# 單樣品
dt_single <- copy(dt[is_double == FALSE])
if (nrow(dt_single) > 0) {
  dt_single[, stem_part := "single"]
}

# 雙樣品：直接複製兩份，不用 by 指派
dt_double_base <- copy(dt[is_double == TRUE])

dt_double_stem1 <- copy(dt_double_base)
if (nrow(dt_double_stem1) > 0) {
  dt_double_stem1[, stem_part := "stem1"]
}

dt_double_stem2 <- copy(dt_double_base)
if (nrow(dt_double_stem2) > 0) {
  dt_double_stem2[, stem_part := "stem2"]
}

dt2 <- rbindlist(
  list(dt_single, dt_double_stem1, dt_double_stem2),
  use.names = TRUE,
  fill = TRUE
)

dt2[, src_dir := file.path(
  sample_dir,
  grid_dir,
  paste0("kmeans_k", source_k_fmt, "_raw_out")
)]

dt2[, subdir_name := fifelse(
  stem_part == "single",
  paste0("subset_from_k", source_k_fmt, "_c", select_cluster_fmt),
  paste0("subset_", stem_part, "_from_k", source_k_fmt, "_c", select_cluster_fmt)
)]

dt2[, subdir := file.path(sample_dir, grid_dir, subdir_name)]

out <- dt2[, .(
  sample_id,
  label,
  grid_dir,
  source_k,
  source_k_fmt,
  select_cluster,
  select_cluster_fmt,
  stem_part,
  sample_dir,
  src_dir,
  subdir_name,
  subdir
)]

fwrite(out, task_tsv, sep = "\t", quote = FALSE, na = "")
cat("Task table written to:", task_tsv, "\n")
print(out)
EOF

# ============================================================
# Step 1 + Step 2. loop all tasks
# Step 1: split stem if needed, then subset one cluster only
# Step 2: run subset kmeans only K=2
# Step 3: plot using stem-specific subset only
# ============================================================
tail -n +2 "${TASK_TSV}" | while IFS=$'\t' read -r sample_id label grid_dir source_k source_k_fmt select_cluster select_cluster_fmt stem_part sample_dir src_dir subdir_name subdir
do
  if [[ -z "${sample_id}" ]]; then
    continue
  fi

  echo "================================================="
  echo "Sample ID       : ${sample_id}"
  echo "Label           : ${label}"
  echo "Grid dir        : ${grid_dir}"
  echo "Source k        : ${source_k}"
  echo "Source k fmt    : ${source_k_fmt}"
  echo "Select cluster  : ${select_cluster}"
  echo "Stem part       : ${stem_part}"
  echo "Sample dir      : ${sample_dir}"
  echo "Source dir      : ${src_dir}"
  echo "Subset dir      : ${subdir}"
  echo "================================================="

  mkdir -p "${subdir}"

  export SAMPLE_DIR="${sample_dir}"
  export GRID_DIR="${grid_dir}"
  export SRC_DIR="${src_dir}"
  export SUB_DIR="${subdir}"
  export SELECT_CLUSTER="${select_cluster}"
  export STEM_PART="${stem_part}"

  # ----------------------------------------------------------
  # Step 1. split stem if needed, then subset one cluster only
  # ----------------------------------------------------------
  Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(sp)
})

sample_dir      <- Sys.getenv("SAMPLE_DIR")
grid_dir        <- Sys.getenv("GRID_DIR")
src_dir         <- Sys.getenv("SRC_DIR")
subdir          <- Sys.getenv("SUB_DIR")
select_cluster  <- Sys.getenv("SELECT_CLUSTER")
stem_part       <- Sys.getenv("STEM_PART")

mat_file   <- file.path(sample_dir, grid_dir, "counts_bins_by_genes_sparse.rds")
meta_file  <- file.path(sample_dir, grid_dir, "bin_metadata.tsv")
clu_file   <- file.path(src_dir, "bin_metadata_with_cluster_raw.tsv")

# selection files 在 sample_dir，不在 grid_dir
upper_file  <- file.path(sample_dir, "Selection_Upper_coordinates.csv")
bottom_file <- file.path(sample_dir, "Selection_Bottom_coordinates.csv")

cat("mat_file   =", mat_file, "\n")
cat("meta_file  =", meta_file, "\n")
cat("clu_file   =", clu_file, "\n")
cat("upper_file =", upper_file, "\n")
cat("bottom_file=", bottom_file, "\n")
cat("stem_part  =", stem_part, "\n")
cat("select_cluster =", select_cluster, "\n")

if (!file.exists(mat_file))  stop("Missing matrix: ", mat_file)
if (!file.exists(meta_file)) stop("Missing meta: ", meta_file)
if (!file.exists(clu_file))  stop("Missing cluster file: ", clu_file)

mat  <- readRDS(mat_file)
meta <- fread(meta_file)
clu  <- fread(clu_file)

stopifnot("bin_id" %in% names(meta), "bin_id" %in% names(clu))
stopifnot("x_center" %in% names(meta), "y_center" %in% names(meta))

read_poly_xy <- function(f) {
  dt <- fread(f)

  nms <- names(dt)
  nms_low <- tolower(nms)

  x_idx <- which(nms_low %in% c("x", "x_center", "xcoord", "x_coord"))[1]
  y_idx <- which(nms_low %in% c("y", "y_center", "ycoord", "y_coord"))[1]

  if (!is.na(x_idx) && !is.na(y_idx)) {
    out <- data.table(
      x = as.numeric(dt[[x_idx]]),
      y = as.numeric(dt[[y_idx]])
    )
  } else {
    if (ncol(dt) < 2) stop("Selection file has fewer than 2 columns: ", f)
    out <- data.table(
      x = as.numeric(dt[[1]]),
      y = as.numeric(dt[[2]])
    )
  }

  out <- out[!is.na(x) & !is.na(y)]

  if (nrow(out) < 3) {
    stop("Selection polygon has fewer than 3 valid points: ", f)
  }

  if (!(out$x[1] == out$x[nrow(out)] && out$y[1] == out$y[nrow(out)])) {
    out <- rbind(out, out[1])
  }

  out
}

has_double_selection <- file.exists(upper_file) && file.exists(bottom_file)
cat("has_double_selection =", has_double_selection, "\n")

meta_use <- copy(meta)
clu_use  <- copy(clu)

if (stem_part %in% c("stem1", "stem2")) {
  if (!has_double_selection) {
    stop("Expected double-sample selection files but missing:\n",
         upper_file, "\n", bottom_file)
  }

  poly_dt <- if (stem_part == "stem1") read_poly_xy(upper_file) else read_poly_xy(bottom_file)

  pip <- point.in.polygon(
    point.x = meta$x_center,
    point.y = meta$y_center,
    pol.x   = poly_dt$x,
    pol.y   = poly_dt$y
  )

  keep_stem_ids <- meta[pip > 0, bin_id]
  cat("Bins inside", stem_part, "polygon:", length(keep_stem_ids), "\n")

  if (length(keep_stem_ids) == 0) {
    stop("No bins found inside polygon for ", stem_part)
  }

  meta_use <- meta[bin_id %in% keep_stem_ids]
  clu_use  <- clu[bin_id %in% keep_stem_ids]
} else {
  cat("No stem split applied.\n")
}

meta2 <- merge(meta_use, clu_use[, .(bin_id, cluster_raw)], by = "bin_id")
meta2[, cluster_raw := as.character(cluster_raw)]

cat("Cluster table within selected region:\n")
print(meta2[, .N, by = cluster_raw][order(suppressWarnings(as.integer(cluster_raw)), cluster_raw)])

keep_ids <- meta2[cluster_raw == select_cluster, bin_id]
cat("Selected bins for cluster", select_cluster, ":", length(keep_ids), "\n")

if (length(keep_ids) == 0) {
  stop("No bins found for selected cluster: ", select_cluster)
}

stopifnot(!is.null(rownames(mat)))

keep_ids2 <- intersect(keep_ids, rownames(mat))
cat("Selected bins present in matrix rownames:", length(keep_ids2), "\n")

if (length(keep_ids2) == 0) {
  stop("None of selected bins are present in matrix rownames.")
}

mat_sub  <- mat[keep_ids2, , drop = FALSE]
meta_sub <- meta2[bin_id %in% keep_ids2]

dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
saveRDS(mat_sub, file.path(subdir, "counts_subset.rds"))
fwrite(meta_sub, file.path(subdir, "binmeta_subset.tsv"), sep = "\t", quote = FALSE)

cat("Subset written to:", subdir, "\n")
cat("counts_subset.rds rows:", nrow(mat_sub), " cols:", ncol(mat_sub), "\n")
cat("binmeta_subset.tsv rows:", nrow(meta_sub), "\n")
EOF

  if [[ ! -f "${subdir}/counts_subset.rds" ]]; then
    echo "[ERROR] counts_subset.rds not found: ${subdir}/counts_subset.rds"
    echo "[SKIP] ${sample_id} ${stem_part}"
    continue
  fi

  if [[ ! -f "${subdir}/binmeta_subset.tsv" ]]; then
    echo "[ERROR] binmeta_subset.tsv not found: ${subdir}/binmeta_subset.tsv"
    echo "[SKIP] ${sample_id} ${stem_part}"
    continue
  fi

  # ----------------------------------------------------------
  # Step 2. run subset kmeans only k=2
  # ----------------------------------------------------------
  k=2
  k_fmt="02"
  outdir="${subdir}/kmeans_subset_k${k_fmt}_raw_out"
  mkdir -p "${outdir}"

  echo "Running subset k=${k_fmt} for ${sample_id} ${stem_part}"

  #Rscript "${KMEANS_SCRIPT}" \
    #--matrix "${subdir}/counts_subset.rds" \
    #--binmeta "${subdir}/binmeta_subset.tsv" \
    #--outdir "${outdir}" \
    #--k "${k}" \
    #--n_pcs 30

  if [[ ! -f "${outdir}/bin_metadata_with_cluster_raw.tsv" ]]; then
    echo "[ERROR] kmeans output missing: ${outdir}/bin_metadata_with_cluster_raw.tsv"
    echo "[SKIP PLOT] ${sample_id} ${stem_part} k=${k_fmt}"
    continue
  fi

  export OUTDIR="${outdir}"
  export SUBDIR_NOW="${subdir}"

  # ----------------------------------------------------------
  # Step 3. plot using stem-specific subset only
  # all_recolor.png: both clusters same color #FF62BC
  # major cluster: #A40000
  # minor cluster: #FFFB00
  # ----------------------------------------------------------
  Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

outdir <- Sys.getenv("OUTDIR")
subdir <- Sys.getenv("SUBDIR_NOW")

f_clu <- file.path(outdir, "bin_metadata_with_cluster_raw.tsv")
f_sub <- file.path(subdir, "binmeta_subset.tsv")

if (!file.exists(f_clu)) stop("Missing: ", f_clu)
if (!file.exists(f_sub)) stop("Missing: ", f_sub)

clu_dt <- fread(f_clu)
sub_dt <- fread(f_sub)

required_sub <- c("bin_id", "x_center", "y_center")
required_clu <- c("bin_id", "cluster_raw")
missing_sub <- setdiff(required_sub, names(sub_dt))
missing_clu <- setdiff(required_clu, names(clu_dt))

if (length(missing_sub) > 0) {
  stop("binmeta_subset.tsv 缺少欄位: ", paste(missing_sub, collapse = ", "))
}
if (length(missing_clu) > 0) {
  stop("bin_metadata_with_cluster_raw.tsv 缺少欄位: ", paste(missing_clu, collapse = ", "))
}

plot_dt <- merge(
  sub_dt[, .(bin_id, x_center, y_center)],
  clu_dt[, .(bin_id, cluster_raw)],
  by = "bin_id"
)

plot_dt <- plot_dt[!is.na(cluster_raw)]
plot_dt[, cluster_raw := as.character(cluster_raw)]

tab <- plot_dt[, .N, by = cluster_raw][order(-N, cluster_raw)]
if (nrow(tab) != 2) {
  stop("K=2 expected exactly 2 clusters, found: ", nrow(tab))
}

major_cl <- tab$cluster_raw[1]
minor_cl <- tab$cluster_raw[2]

infer_bin_size <- function(v) {
  u <- sort(unique(v))
  d <- diff(u)
  d <- d[d > 0]
  if (length(d) == 0) return(1)
  median(d)
}

bin_um_x <- infer_bin_size(plot_dt$x_center)
bin_um_y <- infer_bin_size(plot_dt$y_center)
bin_um <- median(c(bin_um_x, bin_um_y))

plot_lim <- list(
  xlim = c(min(plot_dt$x_center) - bin_um / 2, max(plot_dt$x_center) + bin_um / 2),
  ylim = c(min(plot_dt$y_center) - bin_um / 2, max(plot_dt$y_center) + bin_um / 2)
)

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

# all clusters together
p_all <- ggplot(plot_dt, aes(x = x_center, y = y_center)) +
  geom_tile(width = bin_um, height = bin_um, fill = "#FF62BC") +
  coord_equal(
    xlim = plot_lim$xlim,
    ylim = rev(plot_lim$ylim),
    expand = FALSE
  ) +
  scale_y_reverse() +
  black_theme() +
  labs(title = "Cambium (Before K)")

ggsave(
  file.path(outdir, "subset_Cambium_beforeK_recolor.png"),
  p_all, width = 6, height = 6, dpi = 300, bg = "black"
)

ggsave(
  file.path(outdir, "subset_Cambium_beforeK_recolor.pdf"),
  p_all, width = 6, height = 6, bg = "black"
)

# major cluster
dt_major <- plot_dt[cluster_raw == major_cl]
p_major <- ggplot() +
  geom_tile(
    data = plot_dt,
    aes(x = x_center, y = y_center),
    width = bin_um, height = bin_um,
    fill = "black"
  ) +
  geom_tile(
    data = dt_major,
    aes(x = x_center, y = y_center),
    width = bin_um, height = bin_um,
    fill = "#A40000"
  ) +
  coord_equal(
    xlim = plot_lim$xlim,
    ylim = rev(plot_lim$ylim),
    expand = FALSE
  ) +
  scale_y_reverse() +
  black_theme() +
  labs(title = paste0("Subset k=2 cluster=", major_cl, " (Fusiform Cambium)"))

ggsave(
  file.path(outdir, paste0("subset_k02_cluster_Fusiform_Cambium_only_recolor.png")),
  p_major, width = 6, height = 6, dpi = 300, bg = "black"
)

ggsave(
  file.path(outdir, paste0("subset_k02_cluster_Fusiform_Cambium_only_recolor.pdf")),
  p_major, width = 6, height = 6, bg = "black"
)

# minor cluster
dt_minor <- plot_dt[cluster_raw == minor_cl]
p_minor <- ggplot() +
  geom_tile(
    data = plot_dt,
    aes(x = x_center, y = y_center),
    width = bin_um, height = bin_um,
    fill = "black"
  ) +
  geom_tile(
    data = dt_minor,
    aes(x = x_center, y = y_center),
    width = bin_um, height = bin_um,
    fill = "#FFFB00"
  ) +
  coord_equal(
    xlim = plot_lim$xlim,
    ylim = rev(plot_lim$ylim),
    expand = FALSE
  ) +
  scale_y_reverse() +
  black_theme() +
  labs(title = paste0("Subset k=2 cluster=", minor_cl, " (Ray Cambium)"))

ggsave(
  file.path(outdir, paste0("subset_k02_cluster_Ray_Cambium_only_recolor.png")),
  p_minor, width = 6, height = 6, dpi = 300, bg = "black"
)

ggsave(
  file.path(outdir, paste0("subset_k02_cluster_Ray_Cambium_only_recolor.pdf")),
  p_minor, width = 6, height = 6, bg = "black"
)



# combined k=2 clusters with original colors
plot_dt[, cluster_group := fifelse(
  cluster_raw == major_cl, "Fusiform Cambium",
  fifelse(cluster_raw == minor_cl, "Ray Cambium", NA_character_)
)]

plot_dt[, cluster_group := factor(
  cluster_group,
  levels = c("Fusiform Cambium", "Ray Cambium")
)]

p_combined <- ggplot(plot_dt, aes(x = x_center, y = y_center, fill = cluster_group)) +
  geom_tile(width = bin_um, height = bin_um) +
  scale_fill_manual(
    values = c(
      "Fusiform Cambium" = "#A40000",
      "Ray Cambium" = "#FFFB00"
    ),
    drop = FALSE
  ) +
  coord_equal(
    xlim = plot_lim$xlim,
    ylim = rev(plot_lim$ylim),
    expand = FALSE
  ) +
  scale_y_reverse() +
  black_theme() +
  labs(title = "Subset k=2 combined (Fusiform + Ray Cambium)") +
  theme(
    legend.position = "none"
  )

ggsave(
  file.path(outdir, "subset_k02_cluster_Fusiform_Ray_Cambium_combined_recolor2.png"),
  p_combined, width = 6, height = 6, dpi = 300, bg = "black"
)

ggsave(
  file.path(outdir, "subset_k02_cluster_Fusiform_Ray_Cambium_combined_recolor2.pdf"),
  p_combined, width = 6, height = 6, bg = "black"
)

cat("Major cluster:", major_cl, " color=#F8A078\n")
cat("Minor cluster:", minor_cl, " color=#FA6464\n")
cat("All bins color=#FF62BC\n")
cat("Plotted bin count:", nrow(plot_dt), "\n")
EOF

done

echo "================================================="
echo "All samples finished."
echo "================================================="




