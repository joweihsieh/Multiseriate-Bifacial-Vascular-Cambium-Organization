#!/usr/bin/env bash
set -euo pipefail

for um in 5
do
  grid_dir=$(printf "grid%02dum_out" "$um")

  if [ ! -d "${grid_dir}" ]; then
    echo "[skip] ${grid_dir} not found"
    continue
  fi

  echo "=============================="
  echo "Grid = ${um} um  (${grid_dir})"
  echo "=============================="

  Rscript - <<EOF
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(readxl)
})

sample_id_run <- basename(getwd())
bin_um        <- ${um}
grid_name     <- sprintf("grid%02dum_out", bin_um)

EXCEL_XLSX <- "/home/woodydrylab/FileShare/20260121_Xenium/Cambium_fusiform.xlsx"

BACKGROUND_COLOR <- "black"

# assign colors
custom_colors <- c(
  "1" = "#F2B700",
  "2" = "#377EB8",
  "3" = "#9467BD",
  "4" = "#F8766D"
)

#SUBSET_K_VALUES <- 2:4
SUBSET_K_VALUES <- 4

# ============================================================
# helper
# ============================================================
std_colnames <- function(dt) {
  nms <- names(dt)
  nms2 <- tolower(nms)
  nms2 <- gsub("[[:space:]]+", "_", nms2)
  nms2 <- gsub("[^a-z0-9_]", "_", nms2)
  nms2 <- gsub("_+", "_", nms2)
  nms2 <- gsub("^_|_$", "", nms2)
  setnames(dt, old = nms, new = nms2)
  dt
}

rename_first_match <- function(dt, target, candidates) {
  hit <- intersect(candidates, names(dt))
  if (!(target %in% names(dt)) && length(hit) >= 1) {
    setnames(dt, hit[1], target)
  }
  dt
}

make_bg_theme <- function() {
  theme_void() +
    theme(
      plot.background   = element_rect(fill = "black", color = NA),
      panel.background  = element_rect(fill = "black", color = NA),
      legend.background = element_rect(fill = "black", color = NA),
      legend.key        = element_rect(fill = "black", color = NA),
      legend.text       = element_text(color = "white"),
      legend.title      = element_text(color = "white"),
      plot.title        = element_text(color = "white", hjust = 0.5)
    )
}

read_meta_safely <- function(meta_file) {
  if (!file.exists(meta_file)) stop("Missing metadata file: ", meta_file)

  dt <- fread(meta_file)

  required_cols <- c("x_center", "y_center", "cluster_raw")
  if (!all(required_cols %in% names(dt))) {
    stop(
      "Missing required columns in metadata. Required = ",
      paste(required_cols, collapse = ", "),
      " | Found = ",
      paste(names(dt), collapse = ", ")
    )
  }

  dt[, x_center := as.numeric(x_center)]
  dt[, y_center := as.numeric(y_center)]
  dt[, cluster_raw := trimws(as.character(cluster_raw))]
  unique(dt)
}

pick_excel_row <- function(excel_sub, mode_tag) {
  if (nrow(excel_sub) == 0) return(excel_sub)

  if (mode_tag == "single") {
    return(excel_sub[1])
  }

  if ("label" %in% names(excel_sub)) {
    lab <- trimws(as.character(excel_sub\$label))

    hit <- excel_sub[grepl(paste0(mode_tag, "\$"), lab, ignore.case = TRUE)]
    if (nrow(hit) > 0) return(hit[1])

    hit <- excel_sub[sapply(
      strsplit(lab, ",", fixed = TRUE),
      function(x) any(grepl(paste0(mode_tag, "\$"), trimws(x), ignore.case = TRUE))
    )]
    if (nrow(hit) > 0) return(hit[1])
  }

  excel_sub[1]
}

plot_cluster_overlay <- function(dt_bg, dt_fg, title_txt, out_png, bin_um, fill_values) {
  dt_fg <- copy(dt_fg)
  dt_fg[, cluster_plot := as.character(cluster_raw)]

  p <- ggplot() +
    geom_tile(
      data = dt_bg,
      aes(x = x_center, y = y_center),
      width = bin_um, height = bin_um,
      fill = BACKGROUND_COLOR
    ) +
    geom_tile(
      data = dt_fg,
      aes(x = x_center, y = y_center, fill = cluster_plot),
      width = bin_um, height = bin_um
    ) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    coord_equal(expand = FALSE) +
    scale_y_reverse() +
    labs(title = title_txt, fill = "cluster") +
    make_bg_theme()

  ggsave(out_png, p, width = 6, height = 6, dpi = 300, bg = "black")
}

plot_cluster_only <- function(dt_fg, title_txt, out_png, bin_um, fill_values) {
  dt_fg <- copy(dt_fg)
  dt_fg[, cluster_plot := as.character(cluster_raw)]

  p <- ggplot() +
    geom_tile(
      data = dt_fg,
      aes(x = x_center, y = y_center, fill = cluster_plot),
      width = bin_um, height = bin_um
    ) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    coord_equal(expand = FALSE) +
    scale_y_reverse() +
    labs(title = title_txt, fill = "cluster") +
    make_bg_theme()

  ggsave(out_png, p, width = 6, height = 6, dpi = 300, bg = "black")
}

get_fill_values <- function(cluster_ids) {
  cluster_ids <- sort(unique(as.character(cluster_ids)))
  miss <- setdiff(cluster_ids, names(custom_colors))
  if (length(miss) > 0) {
    stop(
      "Missing colors for cluster ids: ",
      paste(miss, collapse = ", "),
      " | custom_colors currently has: ",
      paste(names(custom_colors), collapse = ", ")
    )
  }
  custom_colors[cluster_ids]
}

# ============================================================
# read excel
# ============================================================
if (!file.exists(EXCEL_XLSX)) {
  stop("Excel not found: ", EXCEL_XLSX)
}

excel_dt <- as.data.table(read_excel(EXCEL_XLSX))
excel_dt <- std_colnames(excel_dt)

excel_dt <- rename_first_match(excel_dt, "sample_id",        c("sample_id", "sample", "sampleid"))
excel_dt <- rename_first_match(excel_dt, "label",            c("label"))
excel_dt <- rename_first_match(excel_dt, "um",               c("um", "grid", "grid_name"))
excel_dt <- rename_first_match(excel_dt, "k",                c("k"))
excel_dt <- rename_first_match(excel_dt, "cambium_fusiform", c("cambium_fusiform"))

required_excel_cols <- c("sample_id", "um", "k", "cambium_fusiform")
if (!all(required_excel_cols %in% names(excel_dt))) {
  stop(
    "Missing required columns in Excel. Required = ",
    paste(required_excel_cols, collapse = ", "),
    " | Found = ",
    paste(names(excel_dt), collapse = ", ")
  )
}

excel_dt[, sample_id := trimws(as.character(sample_id))]
excel_dt[, um := trimws(as.character(um))]
excel_dt[, k := as.integer(k)]
excel_dt[, cambium_fusiform := as.integer(cambium_fusiform)]
if ("label" %in% names(excel_dt)) {
  excel_dt[, label := trimws(as.character(label))]
}

excel_sub <- excel_dt[sample_id == sample_id_run & um == grid_name]

if (nrow(excel_sub) == 0) {
  message("[skip] No matching excel rows for sample_id = ", sample_id_run, ", um = ", grid_name)
  quit(save = "no", status = 0)
}

# ============================================================
# detect K10 mode
# ============================================================
outdir_k10 <- file.path(grid_name, "kmeans_k10_raw_out")
if (!dir.exists(outdir_k10)) {
  message("[skip] k10 output directory not found: ", outdir_k10)
  quit(save = "no", status = 0)
}

meta_k10_stem1  <- file.path(outdir_k10, "bin_metadata_with_cluster_raw_stem1.tsv")
meta_k10_stem2  <- file.path(outdir_k10, "bin_metadata_with_cluster_raw_stem2.tsv")
meta_k10_single <- file.path(outdir_k10, "bin_metadata_with_cluster_raw.tsv")

mode_dt <- NULL

if (file.exists(meta_k10_stem1) && file.exists(meta_k10_stem2)) {
  mode_dt <- data.table(
    mode_tag = c("stem1", "stem2"),
    meta_k10_file = c(meta_k10_stem1, meta_k10_stem2),
    out_suffix = c("_stem1", "_stem2")
  )
  message("Detected split mode: stem1/stem2")
} else if (file.exists(meta_k10_single)) {
  mode_dt <- data.table(
    mode_tag = "single",
    meta_k10_file = meta_k10_single,
    out_suffix = ""
  )
  message("Detected single-file mode")
} else {
  message("[skip] No usable k10 metadata found in: ", outdir_k10)
  quit(save = "no", status = 0)
}

summary_list <- list()

# ============================================================
# main loop
# ============================================================
for (ii in seq_len(nrow(mode_dt))) {

  mode_tag      <- mode_dt\$mode_tag[ii]
  meta_k10_file <- mode_dt\$meta_k10_file[ii]
  out_suffix    <- mode_dt\$out_suffix[ii]

  excel_one <- pick_excel_row(excel_sub, mode_tag)
  if (nrow(excel_one) == 0) {
    message("[skip] No excel row for ", sample_id_run, " ", mode_tag)
    next
  }

  k_base <- as.integer(excel_one\$k[1])
  cambium_fusiform_target <- as.integer(excel_one\$cambium_fusiform[1])

  if (is.na(k_base) || is.na(cambium_fusiform_target)) {
    message("[skip] Invalid k or cambium_fusiform for ", sample_id_run, " ", mode_tag)
    next
  }

  if (k_base != 10) {
    message("[warn] Excel k is not 10 for ", sample_id_run, " ", mode_tag, " | found: ", k_base)
  }

  subset_dir <- file.path(
    grid_name,
    sprintf("subset_from_k10_c%02d", cambium_fusiform_target)
  )

  if (!dir.exists(subset_dir)) {
    message("[skip] Cannot find subset dir: ", subset_dir)
    next
  }

  dt_k10 <- read_meta_safely(meta_k10_file)
  dt_k10_bins <- unique(dt_k10[, .(x_center, y_center, cluster_raw_k10 = cluster_raw)])

  result_dir <- file.path(subset_dir, paste0("projected_from_k10", out_suffix))
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

  # ------------------------------------------------------------
  # helper for each subset-k
  # ------------------------------------------------------------
  process_subset_k <- function(meta_subset_file, subset_k) {

    if (!file.exists(meta_subset_file)) return(NULL)

    dt_sub <- read_meta_safely(meta_subset_file)
    dt_sub2 <- unique(dt_sub[, .(x_center, y_center, cluster_raw)])

    dt_merge <- merge(
      dt_k10_bins,
      dt_sub2,
      by = c("x_center", "y_center"),
      all.x = TRUE
    )

    dt_hit <- dt_merge[!is.na(cluster_raw)]

    if (nrow(dt_hit) == 0) {
      message("[skip] No overlapping bins found for subset k=", subset_k, " in ", sample_id_run, " ", mode_tag)
      return(NULL)
    }

    dt_hit[, subset_cluster := as.character(cluster_raw)]
    dt_hit[, sample_id := sample_id_run]
    dt_hit[, mode := mode_tag]
    dt_hit[, um := grid_name]
    dt_hit[, k_base := sprintf("%02d", k_base)]
    dt_hit[, cambium_fusiform := sprintf("%02d", cambium_fusiform_target)]
    dt_hit[, subset_k := sprintf("%02d", subset_k)]

    out_tsv <- file.path(
      result_dir,
      sprintf("merged_k10_to_subset_k%02d%s.tsv", subset_k, out_suffix)
    )

    out_dt <- dt_hit[, .(
      sample_id, mode, um, k_base, cambium_fusiform, subset_k,
      x_center, y_center, cluster_raw_k10, subset_cluster
    )]
    fwrite(out_dt, out_tsv, sep = "\t")

    fill_values <- get_fill_values(dt_hit\$subset_cluster)

    # overlay
    title1 <- paste0(
      sample_id_run, out_suffix,
      " | K10 c", sprintf("%02d", cambium_fusiform_target),
      " -> subset k=", sprintf("%02d", subset_k)
    )

    out_png1 <- file.path(
      result_dir,
      sprintf("overlay_k10_to_subset_k%02d%s.png", subset_k, out_suffix)
    )

    plot_cluster_overlay(
      dt_bg = dt_k10_bins[, .(x_center, y_center)],
      dt_fg = dt_hit[, .(x_center, y_center, cluster_raw = subset_cluster)],
      title_txt = title1,
      out_png = out_png1,
      bin_um = bin_um,
      fill_values = fill_values
    )

    # subset only
    title2 <- paste0(
      sample_id_run, out_suffix,
      " | subset only | K10 c", sprintf("%02d", cambium_fusiform_target),
      " -> k=", sprintf("%02d", subset_k)
    )

    out_png2 <- file.path(
      result_dir,
      sprintf("subset_only_k%02d%s.png", subset_k, out_suffix)
    )

    plot_cluster_only(
      dt_fg = dt_hit[, .(x_center, y_center, cluster_raw = subset_cluster)],
      title_txt = title2,
      out_png = out_png2,
      bin_um = bin_um,
      fill_values = fill_values
    )

    cluster_summary <- dt_hit[, .N, by = .(subset_cluster)][order(as.integer(subset_cluster))]
    summary_file <- file.path(
      result_dir,
      sprintf("cluster_count_k%02d%s.tsv", subset_k, out_suffix)
    )
    fwrite(cluster_summary, summary_file, sep = "\t")

    message("Saved: ", out_png1)
    message("Saved: ", out_png2)
    message("Saved: ", out_tsv)
    message("Saved: ", summary_file)

    data.table(
      sample_id = sample_id_run,
      mode = mode_tag,
      um = grid_name,
      k_base = sprintf("%02d", k_base),
      cambium_fusiform = sprintf("%02d", cambium_fusiform_target),
      subset_k = sprintf("%02d", subset_k),
      n_k10_bins = nrow(dt_k10_bins),
      n_subset_bins = nrow(dt_hit),
      out_png_overlay = out_png1,
      out_png_only = out_png2,
      out_tsv = out_tsv,
      cluster_count_tsv = summary_file
    )
  }

  for (subset_k in SUBSET_K_VALUES) {

    meta_subset_file <- file.path(
      subset_dir,
      sprintf("kmeans_subset_k%02d_raw_out", subset_k),
      "bin_metadata_with_cluster_raw.tsv"
    )

    if (!file.exists(meta_subset_file)) {
      message("[skip] Missing subset metadata for k=", subset_k, ": ", meta_subset_file)
      next
    }

    res <- process_subset_k(meta_subset_file, subset_k)

    if (!is.null(res)) {
      summary_list[[length(summary_list) + 1]] <- res
    }
  }
}

if (length(summary_list) > 0) {
  summary_dt <- rbindlist(summary_list, fill = TRUE)
  summary_file <- file.path(outdir_k10, "projected_cambium_fusiform_summary.tsv")
  fwrite(summary_dt, summary_file, sep = "\t")
  message("Saved: ", summary_file)
}

cat("Done.\\n")
EOF

done


##############################
# Run each output-* file
##############################
#for d in /home/woodydrylab/FileShare/20260121_Xenium/output-*
#do
#  cd "$d"
#
#  if [ ! -d "grid05um_out" ]; then
#    echo "[skip] grid05um_out not found in $d"
#    continue
#  fi

#  bash /home/woodydrylab/FileShare/20260121_Xenium/plot_cambium_fusiform_based_on_k10_v3.sh
#done