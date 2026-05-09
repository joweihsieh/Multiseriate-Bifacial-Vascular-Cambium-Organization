#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

############################################################
# user settings
############################################################
INFILE <- "/home/woodydrylab/FileShare/20260121_Xenium/xenium_sc_with_controls_nonzero_mean_multi_metric/observed_metrics.tsv"
OUTDIR <- "/home/woodydrylab/FileShare/20260121_Xenium/xenium_sc_with_controls_nonzero_mean_multi_metric/observed_only_pairwise_stats"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OUT_PLOT_PDF_RAW           <- file.path(OUTDIR, "observed_euclidean_boxplot_raw.pdf")
OUT_PLOT_PNG_RAW           <- file.path(OUTDIR, "observed_euclidean_boxplot_raw.png")
OUT_PLOT_PDF_NOOUT         <- file.path(OUTDIR, "observed_euclidean_boxplot_no_outlier.pdf")
OUT_PLOT_PNG_NOOUT         <- file.path(OUTDIR, "observed_euclidean_boxplot_no_outlier.png")

OUT_SUMMARY_RAW_TSV        <- file.path(OUTDIR, "observed_euclidean_summary_raw.tsv")
OUT_SUMMARY_NOOUT_TSV      <- file.path(OUTDIR, "observed_euclidean_summary_no_outlier.tsv")
OUT_OUTLIER_SUMMARY_TSV    <- file.path(OUTDIR, "observed_euclidean_outlier_summary.tsv")

OUT_TTEST_RAW_TSV          <- file.path(OUTDIR, "observed_euclidean_pairwise_ttest_raw.tsv")
OUT_WILCOX_RAW_TSV         <- file.path(OUTDIR, "observed_euclidean_pairwise_wilcox_raw.tsv")
OUT_TTEST_NOOUT_TSV        <- file.path(OUTDIR, "observed_euclidean_pairwise_ttest_no_outlier.tsv")
OUT_WILCOX_NOOUT_TSV       <- file.path(OUTDIR, "observed_euclidean_pairwise_wilcox_no_outlier.tsv")

############################################################
# helper
############################################################
stopf <- function(...) stop(sprintf(...), call. = FALSE)

save_metric_boxplot2 <- function(dt, metric_col, outfile_pdf, outfile_png,
                                 title_txt = "", ylab_txt = "Euclidean distance",
                                 width = 5, height = 7) {
  p <- ggplot(dt, aes(x = comparison, y = get(metric_col))) +
    geom_boxplot(width = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(
      x = NULL,
      y = ylab_txt,
      title = title_txt
    )

  ggsave(outfile_pdf, p, width = width, height = height)
  ggsave(outfile_png, p, width = width, height = height, dpi = 300)
}

pairwise_matrix_to_dt <- function(pmat) {
  out <- as.data.table(as.table(pmat))
  setnames(out, c("group1", "group2", "p_value"))
  out[!is.na(p_value)]
}

make_summary <- function(dt) {
  dt[, .(
    n = .N,
    mean = mean(euclidean_distance, na.rm = TRUE),
    median = median(euclidean_distance, na.rm = TRUE),
    sd = sd(euclidean_distance, na.rm = TRUE),
    min = min(euclidean_distance, na.rm = TRUE),
    max = max(euclidean_distance, na.rm = TRUE)
  ), by = comparison]
}

############################################################
# step 1. read data
############################################################
if (!file.exists(INFILE)) {
  stopf("Input file not found: %s", INFILE)
}

dt <- fread(INFILE)

req_cols <- c("xenium_group", "sc_group", "control_type", "euclidean_distance")
miss_cols <- setdiff(req_cols, names(dt))
if (length(miss_cols) > 0) {
  stopf("Missing required columns: %s", paste(miss_cols, collapse = ", "))
}

############################################################
# step 2. keep observed-only four groups
############################################################
dt <- dt[control_type == "observed"]
dt[, comparison := paste(xenium_group, "vs", sc_group)]

obs_levels4 <- c(
  "ray initials vs ray organizer",
  "ray initials vs fusiform organizer",
  "fusiform initials vs fusiform organizer",
  "fusiform initials vs ray organizer"
)

dt <- dt[comparison %in% obs_levels4]
dt[, comparison := factor(comparison, levels = obs_levels4)]

if (nrow(dt) == 0) {
  stopf("No rows remaining after filtering to observed-only four groups.")
}

############################################################
# step 3. raw summary and raw plot
############################################################
summary_raw <- make_summary(dt)
fwrite(summary_raw, OUT_SUMMARY_RAW_TSV, sep = "\t")

save_metric_boxplot2(
  dt = dt,
  metric_col = "euclidean_distance",
  outfile_pdf = OUT_PLOT_PDF_RAW,
  outfile_png = OUT_PLOT_PNG_RAW,
  title_txt = "",
  ylab_txt = "Euclidean distance"
)

############################################################
# step 4. remove outliers within each comparison
############################################################
dt_no_outlier <- dt[, {
  q1  <- quantile(euclidean_distance, 0.25, na.rm = TRUE)
  q3  <- quantile(euclidean_distance, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  low <- q1 - 1.5 * iqr
  high <- q3 + 1.5 * iqr

  .SD[euclidean_distance >= low & euclidean_distance <= high]
}, by = comparison]

dt_no_outlier[, comparison := factor(as.character(comparison), levels = obs_levels4)]

if (nrow(dt_no_outlier) == 0) {
  stopf("No rows remaining after outlier removal.")
}

############################################################
# step 5. outlier summary
############################################################
n_before <- dt[, .N, by = comparison]
setnames(n_before, "N", "n_before")

n_after <- dt_no_outlier[, .N, by = comparison]
setnames(n_after, "N", "n_after")

outlier_summary <- merge(n_before, n_after, by = "comparison", all = TRUE)
outlier_summary[is.na(n_after), n_after := 0L]
outlier_summary[, n_removed := n_before - n_after]

fwrite(outlier_summary, OUT_OUTLIER_SUMMARY_TSV, sep = "\t")


############################################################
# step 6. pairwise tests after outlier removal
############################################################
ttest_noout <- pairwise.t.test(
  x = dt_no_outlier$euclidean_distance,
  g = dt_no_outlier$comparison,
  p.adjust.method = "none",
  pool.sd = FALSE
)

wilcox_noout <- pairwise.wilcox.test(
  x = dt_no_outlier$euclidean_distance,
  g = dt_no_outlier$comparison,
  p.adjust.method = "none",
  exact = FALSE
)

fwrite(pairwise_matrix_to_dt(ttest_noout$p.value), OUT_TTEST_NOOUT_TSV, sep = "\t")
fwrite(pairwise_matrix_to_dt(wilcox_noout$p.value), OUT_WILCOX_NOOUT_TSV, sep = "\t")

############################################################
# done
############################################################
message("Done.")
message("Output directory: ", OUTDIR)