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

dt[, grid_dir := fifelse(
  grepl("^grid[0-9]+um_out$", as.character(um)),
  as.character(um),
  paste0("grid", sprintf("%02d", as.integer(gsub("[^0-9]", "", as.character(um)))), "um_out")
)]

dt[, k0 := as.integer(k)]
dt[, select_cluster := as.character(as.integer(cambium_fusiform))]
dt[, select_cluster_fmt := sprintf("%02d", as.integer(cambium_fusiform))]

dt[, sample_dir := file.path(base_root, sample_id)]
dt[, src_dir := file.path(sample_dir, grid_dir, paste0("kmeans_k", k0, "_raw_out"))]
dt[, subdir_name := paste0("subset_from_k", k0, "_c", select_cluster_fmt)]
dt[, subdir := file.path(sample_dir, grid_dir, subdir_name)]

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
  echo "Sample ID          : ${sample_id}"
  echo "Label              : ${label}"
  echo "Grid dir           : ${grid_dir}"
  echo "K0                 : ${k0}"
  echo "Select cluster     : ${select_cluster}"
  echo "Select cluster fmt : ${select_cluster_fmt}"
  echo "Subset dir         : ${subdir}"
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

    #Rscript "${KMEANS_SCRIPT}" \
    #  --matrix "${subdir}/counts_subset.rds" \
    #  --binmeta "${subdir}/binmeta_subset.tsv" \
    #  --outdir "${outdir}" \
    #  --k "${k}" \
    #  --n_pcs 30

    export OUTDIR="${outdir}"
    export K_NOW="${k}"
    export SAMPLE_ID_NOW="${sample_id}"
    export SAMPLE_DIR_NOW="${sample_dir}"
    export GRID_DIR_NOW="${grid_dir}"
    export LABEL_NOW="${label}"

    Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

outdir       <- Sys.getenv("OUTDIR")
k            <- as.integer(Sys.getenv("K_NOW"))
sample_id    <- Sys.getenv("SAMPLE_ID_NOW")
sample_dir   <- Sys.getenv("SAMPLE_DIR_NOW")
grid_dir     <- Sys.getenv("GRID_DIR_NOW")
sample_label <- Sys.getenv("LABEL_NOW")

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

custom_colors <- c(
  "1" = "#F2B700",
  "2" = "#377EB8",
  "3" = "#9467BD"
)

clusters_chr <- as.character(clusters)
clusters_num <- suppressWarnings(as.integer(clusters_chr))

if (k %in% c(2, 3) && all(!is.na(clusters_num))) {
  ord <- order(clusters_num)
  ordered_clusters <- clusters_chr[ord]
  pal <- custom_colors[seq_along(ordered_clusters)]
  names(pal) <- ordered_clusters
  pal <- pal[clusters_chr]
} else {
  pal <- hue_pal()(length(clusters_chr))
  names(pal) <- clusters_chr
}

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

find_id_col <- function(dt) {
  cands <- c("bin_id", "barcode", "cell_id", "id")
  hit <- cands[cands %in% names(dt)][1]
  if (is.na(hit)) return(NA_character_)
  hit
}

find_xy_cols <- function(dt) {
  x_cands <- c("x_center", "x", "x_centroid", "center_x", "pixel_x", "x_global")
  y_cands <- c("y_center", "y", "y_centroid", "center_y", "pixel_y", "y_global")
  x_col <- x_cands[x_cands %in% names(dt)][1]
  y_col <- y_cands[y_cands %in% names(dt)][1]
  if (is.na(x_col) || is.na(y_col)) {
    stop("Cannot find x/y columns. Available columns: ", paste(names(dt), collapse = ", "))
  }
  list(x = x_col, y = y_col)
}

split_subset_by_stem <- function(sub_dt, stem1_dt, stem2_dt) {
  sub_xy   <- find_xy_cols(sub_dt)
  stem1_xy <- find_xy_cols(stem1_dt)
  stem2_xy <- find_xy_cols(stem2_dt)

  sub_id   <- find_id_col(sub_dt)
  stem1_id <- find_id_col(stem1_dt)
  stem2_id <- find_id_col(stem2_dt)

  use_id_match <- !is.na(sub_id) && !is.na(stem1_id) && !is.na(stem2_id)

  if (use_id_match) {
    sub_map <- copy(sub_dt)[, .(join_id = as.character(get(sub_id)))]
    stem1_map <- unique(copy(stem1_dt)[, .(join_id = as.character(get(stem1_id)), stem = "stem1")])
    stem2_map <- unique(copy(stem2_dt)[, .(join_id = as.character(get(stem2_id)), stem = "stem2")])
    stem_map <- unique(rbindlist(list(stem1_map, stem2_map)))
    out <- cbind(copy(sub_dt), sub_map)
    out <- merge(out, stem_map, by = "join_id", all.x = TRUE)
  } else {
    sub_map <- copy(sub_dt)[, .(
      join_x = get(sub_xy$x),
      join_y = get(sub_xy$y)
    )]
    stem1_map <- unique(copy(stem1_dt)[, .(
      join_x = get(stem1_xy$x),
      join_y = get(stem1_xy$y),
      stem = "stem1"
    )])
    stem2_map <- unique(copy(stem2_dt)[, .(
      join_x = get(stem2_xy$x),
      join_y = get(stem2_xy$y),
      stem = "stem2"
    )])
    stem_map <- unique(rbindlist(list(stem1_map, stem2_map)))
    out <- cbind(copy(sub_dt), sub_map)
    out <- merge(out, stem_map, by = c("join_x", "join_y"), all.x = TRUE)
  }

  if (any(is.na(out$stem))) {
    message("  Some bins not matched in first pass; trying rounded x/y fallback...")

    sub2 <- copy(sub_dt)
    s12  <- copy(stem1_dt)
    s22  <- copy(stem2_dt)

    sub2[, rx := round(get(sub_xy$x), 4)]
    sub2[, ry := round(get(sub_xy$y), 4)]

    s12[, rx := round(get(stem1_xy$x), 4)]
    s12[, ry := round(get(stem1_xy$y), 4)]
    s12 <- unique(s12[, .(rx, ry, stem = "stem1")])

    s22[, rx := round(get(stem2_xy$x), 4)]
    s22[, ry := round(get(stem2_xy$y), 4)]
    s22 <- unique(s22[, .(rx, ry, stem = "stem2")])

    stem_map2 <- unique(rbindlist(list(s12, s22)))

    out2 <- copy(sub_dt)
    out2[, rx := round(get(sub_xy$x), 4)]
    out2[, ry := round(get(sub_xy$y), 4)]
    out2 <- merge(out2, stem_map2, by = c("rx", "ry"), all.x = TRUE)

    miss_idx <- which(is.na(out$stem))
    if (length(miss_idx) > 0 && "stem" %in% names(out2)) {
      out$stem[miss_idx] <- out2$stem[miss_idx]
    }
  }

  out[]
}

save_plot_set <- function(dt_plot, prefix, label_text) {
  if (nrow(dt_plot) == 0) {
    cat("[skip] no data for", prefix, "\n")
    return(invisible(NULL))
  }

  dt_plot <- copy(dt_plot)
  dt_plot[, cluster_raw := factor(as.character(cluster_raw), levels = clusters_chr)]

  # Zoom to the coordinate range of the current stem or sample
  plot_lim <- get_plot_limits(dt_plot, bin_um = bin_um)
  k_fmt <- sprintf("%02d", k)

  p_all <- ggplot(dt_plot, aes(x = x_center, y = y_center, fill = cluster_raw)) +
    geom_tile(width = bin_um, height = bin_um) +
    coord_equal(
      xlim = plot_lim$xlim,
      ylim = rev(plot_lim$ylim),
      expand = FALSE
    ) +
    scale_fill_manual(values = pal, drop = FALSE) +
    black_theme() +
    labs(title = paste0(label_text, " | subset k=", k_fmt), fill = "Cluster")

  ggsave(
    file.path(outdir, paste0(prefix, "subset_k", k_fmt, "_all.png")),
    p_all, width = 6, height = 6, dpi = 300, bg = "black"
  )
  ggsave(
    file.path(outdir, paste0(prefix, "subset_k", k_fmt, "_all.pdf")),
    p_all, width = 6, height = 6, bg = "black"
  )

  for (cl in clusters_chr) {
    dt_cl <- dt_plot[as.character(cluster_raw) == cl]
    if (nrow(dt_cl) == 0) next

    cl_fmt <- sprintf("%02d", as.integer(cl))

    p_one <- ggplot() +
      geom_tile(
        data = dt_plot,
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
      labs(title = paste0(label_text, " | subset k=", k_fmt, " cluster=", cl_fmt))

    ggsave(
      file.path(outdir, paste0(prefix, "subset_k", k_fmt, "_cluster_", cl_fmt, ".png")),
      p_one, width = 6, height = 6, dpi = 300, bg = "black"
    )
    ggsave(
      file.path(outdir, paste0(prefix, "subset_k", k_fmt, "_cluster_", cl_fmt, ".pdf")),
      p_one, width = 6, height = 6, bg = "black"
    )
  }
}

# ============================================================
# split if two stems exist
# ============================================================
k10_root <- file.path(sample_dir, grid_dir, "kmeans_k10_raw_out")
stem1_file <- file.path(k10_root, "bin_metadata_with_cluster_raw_stem1.tsv")
stem2_file <- file.path(k10_root, "bin_metadata_with_cluster_raw_stem2.tsv")

if (file.exists(stem1_file) && file.exists(stem2_file)) {
  cat("Detected dual-sample folder. Splitting stem1 / stem2 before plotting...\n")
  cat("  stem1:", stem1_file, "\n")
  cat("  stem2:", stem2_file, "\n")

  stem1_dt <- fread(stem1_file)
  stem2_dt <- fread(stem2_file)

  split_dt <- split_subset_by_stem(bin_meta, stem1_dt, stem2_dt)

  n_miss <- sum(is.na(split_dt$stem))
  if (n_miss > 0) {
    cat("[warning]", n_miss, "bins could not be assigned to stem1/stem2\n")
  }

  dt1 <- split_dt[stem == "stem1"]
  dt2 <- split_dt[stem == "stem2"]

  if (nrow(dt1) > 0) {
    save_plot_set(dt1, prefix = "stem1_", label_text = paste0(sample_label, " | stem1"))
  } else {
    cat("[skip] no bins assigned to stem1\n")
  }

  if (nrow(dt2) > 0) {
    save_plot_set(dt2, prefix = "stem2_", label_text = paste0(sample_label, " | stem2"))
  } else {
    cat("[skip] no bins assigned to stem2\n")
  }

} else {
  cat("Single-sample folder detected. Plotting whole subset directly...\n")
  save_plot_set(bin_meta, prefix = "", label_text = sample_label)
}
EOF

  done

done < "${TASK_TSV}"

echo "================================================="
echo "All samples finished."
echo "================================================="




