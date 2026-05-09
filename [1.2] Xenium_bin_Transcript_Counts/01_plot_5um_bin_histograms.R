#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

this_file <- if (!is.null(sys.frames()[[1]]$ofile)) sys.frames()[[1]]$ofile else "."
script_dir <- dirname(normalizePath(this_file, mustWork = FALSE))
source(file.path(script_dir, "bin_qc_helpers.R"))

XENIUM_BASE <- "/home/woodydrylab/FileShare/20260121_Xenium"
GRID_NAME <- "grid05um_out"
OUTDIR <- file.path(XENIUM_BASE, "hist_transcripts_per_bin")
BOTTOM_COORD_FILE <- "Selection_Bottom_coordinates.csv"
UPPER_COORD_FILE <- "Selection_Upper_coordinates.csv"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

sample_dirs <- get_sample_dirs(XENIUM_BASE)

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

    hist_file <- file.path(
      OUTDIR,
      paste0(sample_name, "_", region_name, "_hist_5um.png")
    )

    plot_hist(
      x = dt_region$total_counts,
      outfile = hist_file,
      title = paste0(sample_name, " | ", region_name, " | 5 um"),
      xlab_txt = "Transcripts per 5 um bin"
    )
  }
}

message("Done.")
message("Histogram dir: ", OUTDIR)
