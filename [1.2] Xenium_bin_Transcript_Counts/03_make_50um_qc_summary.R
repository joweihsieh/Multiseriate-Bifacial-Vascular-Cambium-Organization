#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

this_file <- if (!is.null(sys.frames()[[1]]$ofile)) sys.frames()[[1]]$ofile else "."
script_dir <- dirname(normalizePath(this_file, mustWork = FALSE))
source(file.path(script_dir, "bin_qc_helpers.R"))

XENIUM_BASE <- "/home/woodydrylab/FileShare/20260121_Xenium"
GRID_NAME <- "grid05um_out"
ORIG_BIN_UM <- 5
TARGET_BIN_UM <- 50
OUTDIR_HIST <- file.path(XENIUM_BASE, paste0("hist_transcripts_per_", TARGET_BIN_UM, "um_bin"))
OUT_TSV <- file.path(XENIUM_BASE, paste0("Xenium_full_QC_summary_", TARGET_BIN_UM, "um.tsv"))
BOTTOM_COORD_FILE <- "Selection_Bottom_coordinates.csv"
UPPER_COORD_FILE <- "Selection_Upper_coordinates.csv"

dir.create(OUTDIR_HIST, showWarnings = FALSE, recursive = TRUE)

sample_dirs <- get_sample_dirs(XENIUM_BASE)
res_list <- list()

for (sdir in sample_dirs) {
  sample_name <- basename(sdir)
  message("Processing: ", sample_name)

  grid_dir <- file.path(sdir, GRID_NAME)
  meta_path <- file.path(grid_dir, "bin_metadata.tsv")
  sparse_rds <- file.path(grid_dir, "counts_bins_by_genes_sparse.rds")

  if (!file.exists(meta_path)) {
    message("  [skip] missing: ", meta_path)
    next
  }

  if (!file.exists(sparse_rds)) {
    message("  [skip] missing: ", sparse_rds)
    next
  }

  dt <- fread(meta_path)
  mat <- readRDS(sparse_rds)

  if (nrow(dt) != nrow(mat)) {
    stop(
      "Row mismatch in sample ", sample_name,
      ": nrow(bin_metadata.tsv) = ", nrow(dt),
      ", nrow(counts_bins_by_genes_sparse.rds) = ", nrow(mat)
    )
  }

  dt[, row_id := .I]
  regions <- split_bin_meta_by_selection(dt, sdir, BOTTOM_COORD_FILE, UPPER_COORD_FILE)

  for (region_name in names(regions)) {
    dt_region <- regions[[region_name]]

    if (nrow(dt_region) == 0) {
      next
    }

    dt_super <- aggregate_to_superbins(dt_region, target_bin_um = TARGET_BIN_UM)

    if (nrow(dt_super) == 0) {
      next
    }

    dt_super[, n_genes := calc_superbin_gene_counts(.SD, mat)]

    hist_file <- file.path(
      OUTDIR_HIST,
      paste0(sample_name, "_", region_name, "_hist_", TARGET_BIN_UM, "um.png")
    )

    plot_hist(
      x = dt_super$total_counts,
      outfile = hist_file,
      title = paste0(sample_name, " | ", region_name, " | ", TARGET_BIN_UM, " um"),
      xlab_txt = paste0("Transcripts per ", TARGET_BIN_UM, " um bin")
    )

    total_bins <- nrow(dt_super)
    dt_nonzero <- dt_super[total_counts > 0]
    nonzero_bins <- nrow(dt_nonzero)

    t_stat <- calc_distribution(dt_nonzero$total_counts)
    g_stat <- calc_distribution(dt_nonzero$n_genes)

    res_list[[paste0(sample_name, "_", region_name)]] <- data.table(
      sample = sample_name,
      region = region_name,
      orig_bin_um = ORIG_BIN_UM,
      target_bin_um = TARGET_BIN_UM,
      total_bins = total_bins,
      nonzero_bins = nonzero_bins,
      pct_nonzero_bins = ifelse(total_bins > 0, nonzero_bins / total_bins, NA_real_),
      t_min = t_stat$min,
      t_Q1 = t_stat$Q1,
      t_median = t_stat$median,
      t_mean = t_stat$mean,
      t_Q3 = t_stat$Q3,
      t_max = t_stat$max,
      g_min = g_stat$min,
      g_Q1 = g_stat$Q1,
      g_median = g_stat$median,
      g_mean = g_stat$mean,
      g_Q3 = g_stat$Q3,
      g_max = g_stat$max
    )
  }
}

if (length(res_list) == 0) {
  message("No valid samples found.")
} else {
  res <- rbindlist(res_list, fill = TRUE)
  fwrite(res, OUT_TSV, sep = "\t")

  message("Done.")
  message("Histogram dir: ", OUTDIR_HIST)
  message("Summary TSV: ", OUT_TSV)
}
