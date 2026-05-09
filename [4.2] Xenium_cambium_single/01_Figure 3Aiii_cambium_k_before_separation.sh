####### Use K = 10 as the reference
# If bin_metadata_with_cluster_raw_stem1.tsv / stem2.tsv exist, use them for separate projection.
# If they do not exist, use the single file bin_metadata_with_cluster_raw.tsv instead.
# Always format the source K value with two digits.
# k08, k09
# Use sprintf("kmeans_k%02d_raw_out", k_source)
# Logic is correct.
# Use the K10 bin set.
# Merge it with bin_metadata_with_cluster_raw.tsv from K = 8 or K = 9.
# Identify bins that belong to cambium_single in the source K.
# Then plot them according to the K10 stem1 / stem2 or single-sample structure.

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

EXCEL_XLSX <- "/home/woodydrylab/FileShare/20260121_Xenium/Cambium_k_single_furtherk2.xlsx"

HILIGHT_COLOR    <- "#FF62BC"
BACKGROUND_COLOR <- "black"

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

plot_bin_tiles <- function(dt_bg, dt_fg, title_txt, out_png, bin_um) {
  p <- ggplot() +
    geom_tile(
      data = dt_bg,
      aes(x = x_center, y = y_center),
      width = bin_um, height = bin_um,
      fill = BACKGROUND_COLOR
    ) +
    geom_tile(
      data = dt_fg,
      aes(x = x_center, y = y_center),
      width = bin_um, height = bin_um,
      fill = HILIGHT_COLOR
    ) +
    coord_equal(expand = FALSE) +
    scale_y_reverse() +
    labs(title = title_txt) +
    make_bg_theme()

  ggsave(out_png, p, width = 6, height = 6, dpi = 300, bg = "black")
}

plot_only_tiles <- function(dt_fg, title_txt, out_png, bin_um) {
  p <- ggplot() +
    geom_tile(
      data = dt_fg,
      aes(x = x_center, y = y_center),
      width = bin_um, height = bin_um,
      fill = HILIGHT_COLOR
    ) +
    coord_equal(expand = FALSE) +
    scale_y_reverse() +
    labs(title = title_txt) +
    make_bg_theme()

  ggsave(out_png, p, width = 6, height = 6, dpi = 300, bg = "black")
}

# ============================================================
# read excel
# ============================================================
if (!file.exists(EXCEL_XLSX)) {
  stop("Excel not found: ", EXCEL_XLSX)
}

excel_dt <- as.data.table(read_excel(EXCEL_XLSX))
excel_dt <- std_colnames(excel_dt)

excel_dt <- rename_first_match(excel_dt, "sample_id",      c("sample_id", "sample", "sampleid"))
excel_dt <- rename_first_match(excel_dt, "label",          c("label"))
excel_dt <- rename_first_match(excel_dt, "um",             c("um", "grid", "grid_name"))
excel_dt <- rename_first_match(excel_dt, "k",              c("k"))
excel_dt <- rename_first_match(excel_dt, "cambium_single", c("cambium_single"))

required_excel_cols <- c("sample_id", "um", "k", "cambium_single")
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
excel_dt[, cambium_single := trimws(as.character(cambium_single))]
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

  k_source <- as.integer(excel_one\$k[1])
  cambium_single_target <- trimws(as.character(excel_one\$cambium_single[1]))

  if (is.na(k_source) || is.na(cambium_single_target) || cambium_single_target == "") {
    message("[skip] Invalid k or cambium_single for ", sample_id_run, " ", mode_tag)
    next
  }

  outdir_source <- file.path(
    grid_name,
    sprintf("kmeans_k%02d_raw_out", k_source)
  )

  if (!dir.exists(outdir_source)) {
    message("[skip] Cannot find source outdir: ", outdir_source)
    next
  }

  meta_source_file <- file.path(outdir_source, "bin_metadata_with_cluster_raw.tsv")
  if (!file.exists(meta_source_file)) {
    message("[skip] Source metadata not found: ", meta_source_file)
    next
  }

  dt_k10_bins <- read_meta_safely(meta_k10_file)[, .(x_center, y_center)]
  dt_source   <- read_meta_safely(meta_source_file)[, .(x_center, y_center, cluster_raw)]

  dt_proj <- merge(
    dt_k10_bins,
    dt_source,
    by = c("x_center", "y_center"),
    all.x = TRUE
  )

  dt_cambium <- dt_proj[cluster_raw == cambium_single_target]

  title1 <- paste0(
    sample_id_run, out_suffix,
    " | projected from k=", sprintf("%02d", k_source),
    " | cambium_single=", cambium_single_target
  )

  out_png1 <- file.path(
    outdir_source,
    paste0("projected_cambium_single_from_k10", out_suffix, ".png")
  )

  plot_bin_tiles(
    dt_bg = dt_k10_bins,
    dt_fg = dt_cambium,
    title_txt = title1,
    out_png = out_png1,
    bin_um = bin_um
  )

  title2 <- paste0(
    sample_id_run, out_suffix,
    " | only projected cambium bins",
    " | k=", sprintf("%02d", k_source),
    " | cluster=", cambium_single_target
  )

  out_png2 <- file.path(
    outdir_source,
    paste0("projected_cambium_single_only", out_suffix, ".png")
  )

  plot_only_tiles(
    dt_fg = dt_cambium,
    title_txt = title2,
    out_png = out_png2,
    bin_um = bin_um
  )

  out_tsv <- file.path(
    outdir_source,
    paste0("projected_cambium_single_bins", out_suffix, ".tsv")
  )

  out_dt <- dt_cambium[, .(x_center, y_center, cluster_raw)]
  out_dt[, sample_id := sample_id_run]
  out_dt[, mode := mode_tag]
  out_dt[, um := grid_name]
  out_dt[, source_k := sprintf("%02d", k_source)]
  out_dt[, cambium_single := cambium_single_target]

  setcolorder(out_dt, c(
    "sample_id", "mode", "um", "source_k", "cambium_single",
    "x_center", "y_center", "cluster_raw"
  ))

  fwrite(out_dt, out_tsv, sep = "\t")

  summary_list[[length(summary_list) + 1]] <- data.table(
    sample_id = sample_id_run,
    mode = mode_tag,
    um = grid_name,
    source_k = sprintf("%02d", k_source),
    cambium_single = cambium_single_target,
    n_k10_bins = nrow(dt_k10_bins),
    n_projected_bins = nrow(dt_cambium),
    out_png1 = out_png1,
    out_png2 = out_png2,
    out_tsv = out_tsv
  )

  message("Saved: ", out_png1)
  message("Saved: ", out_png2)
  message("Saved: ", out_tsv)
}

if (length(summary_list) > 0) {
  summary_dt <- rbindlist(summary_list, fill = TRUE)
  summary_file <- file.path(outdir_k10, "projected_cambium_single_summary.tsv")
  fwrite(summary_dt, summary_file, sep = "\t")
  message("Saved: ", summary_file)
}

cat("Done.\\n")
EOF

done



##############################
# 跑全部 output-* 資料夾
##############################
for d in /home/woodydrylab/FileShare/20260121_Xenium/output-*
do
  cd "$d"

  if [ ! -d "grid05um_out" ]; then
    echo "[skip] grid05um_out not found in $d"
    continue
  fi

  bash /home/woodydrylab/FileShare/20260121_Xenium/plot_projected_cambium_from_k10.sh
done