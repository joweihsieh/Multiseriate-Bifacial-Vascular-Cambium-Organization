#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

this_file <- if (!is.null(sys.frames()[[1]]$ofile)) sys.frames()[[1]]$ofile else "."
script_dir <- dirname(normalizePath(this_file, mustWork = FALSE))
source(file.path(script_dir, "bin_qc_helpers.R"))

XENIUM_BASE <- "/home/woodydrylab/FileShare/20260121_Xenium"
GRID_NAME <- "grid05um_out"
OUT_TSV <- file.path(XENIUM_BASE, "Xenium_full_QC_summary.tsv")
BOTTOM_COORD_FILE <- "Selection_Bottom_coordinates.csv"
UPPER_COORD_FILE <- "Selection_Upper_coordinates.csv"

sample_dirs <- get_sample_dirs(XENIUM_BASE)
res_list <- list()

for (sdir in sample_dirs) {
  sample_name <- basename(sdir)
  message("Processing: ", sample_name)

  grid_dir <- file.path(sdir, GRID_NAME)
  meta_path <- file.path(grid_dir, "bin_metadata.tsv")

  if (!file.exists(meta_path)) {
    message("  [skip] missing: ", meta_path)
    next
  }

  dt <- fread(meta_path)
  regions <- split_bin_meta_by_selection(dt, sdir, BOTTOM_COORD_FILE, UPPER_COORD_FILE)

  for (region_name in names(regions)) {
    dt_region <- regions[[region_name]]
    total_bins <- nrow(dt_region)
    dt_nonzero <- dt_region[total_counts > 0]
    nonzero_bins <- nrow(dt_nonzero)

    t_stat <- calc_distribution(dt_nonzero$total_counts)
    g_stat <- calc_distribution(dt_nonzero$n_genes)

    res_list[[paste0(sample_name, "_", region_name)]] <- data.table(
      sample = sample_name,
      region = region_name,
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

res <- rbindlist(res_list, fill = TRUE)
fwrite(res, OUT_TSV, sep = "\t")

message("Done.")
message("Saved to: ", OUT_TSV)
