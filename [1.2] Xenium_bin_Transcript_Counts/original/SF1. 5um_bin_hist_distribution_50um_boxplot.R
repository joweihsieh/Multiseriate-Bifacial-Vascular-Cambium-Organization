#!/usr/bin/env Rscript

# This legacy script has been split into smaller task-specific scripts:
# 01_plot_5um_bin_histograms.R
# 02_make_5um_qc_summary.R
# 03_make_50um_qc_summary.R
# 04_plot_50um_gmean_boxplot.R

suppressPackageStartupMessages({
  library(data.table)
})

# =======================
# SETTINGS
# =======================
XENIUM_BASE <- "/home/woodydrylab/FileShare/20260121_Xenium"
GRID_NAME <- "grid05um_out"

OUTDIR <- "hist_transcripts_per_bin"
dir.create(OUTDIR, showWarnings = FALSE)

BOTTOM_COORD_FILE <- "Selection_Bottom_coordinates.csv"
UPPER_COORD_FILE  <- "Selection_Upper_coordinates.csv"

# =======================
# HELPERS
# =======================
guess_xy_cols <- function(dt) {
  nms <- names(dt)
  low <- tolower(nms)

  x_idx <- which(low %in% c("x", "x_center"))
  y_idx <- which(low %in% c("y", "y_center"))

  list(x = nms[x_idx[1]], y = nms[y_idx[1]])
}

read_polygon_csv <- function(path) {
  dt <- fread(path)
  xy <- guess_xy_cols(dt)

  poly <- dt[, .(
    x = as.numeric(get(xy$x)),
    y = as.numeric(get(xy$y))
  )]

  if (!(poly$x[1] == poly$x[nrow(poly)])) {
    poly <- rbind(poly, poly[1])
  }

  poly
}

point_in_polygon_vec <- function(px, py, vx, vy) {
  n <- length(vx)
  inside <- rep(FALSE, length(px))

  j <- n
  for (i in seq_len(n)) {
    xi <- vx[i]; yi <- vy[i]
    xj <- vx[j]; yj <- vy[j]

    intersect <- ((yi > py) != (yj > py)) &
      (px < (xj - xi) * (py - yi) / ((yj - yi) + 1e-30) + xi)

    inside <- xor(inside, intersect)
    j <- i
  }
  inside
}

split_bin_meta_by_selection <- function(bin_meta, sample_dir) {

  bottom_path <- file.path(sample_dir, BOTTOM_COORD_FILE)
  upper_path  <- file.path(sample_dir, UPPER_COORD_FILE)

  if (!file.exists(bottom_path) || !file.exists(upper_path)) {
    return(list(ALL = bin_meta))
  }

  poly_bottom <- read_polygon_csv(bottom_path)
  poly_upper  <- read_polygon_csv(upper_path)

  in_bottom <- point_in_polygon_vec(
    bin_meta$x_center, bin_meta$y_center,
    poly_bottom$x, poly_bottom$y
  )

  in_upper <- point_in_polygon_vec(
    bin_meta$x_center, bin_meta$y_center,
    poly_upper$x, poly_upper$y
  )

  list(
    Bottom = bin_meta[in_bottom],
    Upper  = bin_meta[in_upper]
  )
}

plot_hist <- function(x, outfile, title) {

  x <- x[x > 0]  # Remove zero-count bins

  if (length(x) < 10) return()

  # Limit extreme outliers to keep the plot readable
  xmax <- quantile(x, 0.99)

  png(outfile, width = 1200, height = 900)

  par(bg = "black",
      col = "white",
      col.axis = "white",
      col.lab = "white")

  hist(x[x <= xmax],
       breaks = 100,
       col = "white",
       border = "white",
       main = title,
       xlab = "Transcripts per bin",
       ylab = "Bin count")

  dev.off()
}

######################
summary
######################
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# =======================
# SETTINGS
# =======================
XENIUM_BASE <- "/home/woodydrylab/FileShare/20260121_Xenium"
GRID_NAME <- "grid05um_out"

OUT_TSV <- file.path(XENIUM_BASE, "Xenium_full_QC_summary.tsv")

BOTTOM_COORD_FILE <- "Selection_Bottom_coordinates.csv"
UPPER_COORD_FILE  <- "Selection_Upper_coordinates.csv"

# =======================
# HELPERS
# =======================
guess_xy_cols <- function(dt) {
  nms <- names(dt)
  low <- tolower(nms)

  x_idx <- which(low %in% c("x", "x_center"))
  y_idx <- which(low %in% c("y", "y_center"))

  list(x = nms[x_idx[1]], y = nms[y_idx[1]])
}

read_polygon_csv <- function(path) {
  dt <- fread(path)
  xy <- guess_xy_cols(dt)

  poly <- dt[, .(
    x = as.numeric(get(xy$x)),
    y = as.numeric(get(xy$y))
  )]

  if (!(poly$x[1] == poly$x[nrow(poly)])) {
    poly <- rbind(poly, poly[1])
  }

  poly
}

point_in_polygon_vec <- function(px, py, vx, vy) {
  n <- length(vx)
  inside <- rep(FALSE, length(px))

  j <- n
  for (i in seq_len(n)) {
    xi <- vx[i]; yi <- vy[i]
    xj <- vx[j]; yj <- vy[j]

    intersect <- ((yi > py) != (yj > py)) &
      (px < (xj - xi) * (py - yi) / ((yj - yi) + 1e-30) + xi)

    inside <- xor(inside, intersect)
    j <- i
  }
  inside
}

split_bin_meta_by_selection <- function(bin_meta, sample_dir) {

  bottom_path <- file.path(sample_dir, BOTTOM_COORD_FILE)
  upper_path  <- file.path(sample_dir, UPPER_COORD_FILE)

  if (!file.exists(bottom_path) || !file.exists(upper_path)) {
    return(list(ALL = bin_meta))
  }

  poly_bottom <- read_polygon_csv(bottom_path)
  poly_upper  <- read_polygon_csv(upper_path)

  in_bottom <- point_in_polygon_vec(
    bin_meta$x_center, bin_meta$y_center,
    poly_bottom$x, poly_bottom$y
  )

  in_upper <- point_in_polygon_vec(
    bin_meta$x_center, bin_meta$y_center,
    poly_upper$x, poly_upper$y
  )

  list(
    Bottom = bin_meta[in_bottom],
    Upper  = bin_meta[in_upper]
  )
}

calc_distribution <- function(x) {

  if (length(x) == 0) {
    return(list(min=NA, Q1=NA, median=NA, mean=NA, Q3=NA, max=NA))
  }

  q <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

  list(
    min = q[1],
    Q1 = q[2],
    median = q[3],
    mean = mean(x),
    Q3 = q[4],
    max = q[5]
  )
}

# =======================
# MAIN
# =======================
sample_dirs <- list.dirs(XENIUM_BASE, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^.*/output-", sample_dirs)]

res_list <- list()

for (sdir in sample_dirs) {

  sample_name <- basename(sdir)
  message("Processing: ", sample_name)

  grid_dir <- file.path(sdir, GRID_NAME)
  meta_path <- file.path(grid_dir, "bin_metadata.tsv")

  if (!file.exists(meta_path)) next

  dt <- fread(meta_path)

  regions <- split_bin_meta_by_selection(dt, sdir)

  for (region_name in names(regions)) {

    dt_region <- regions[[region_name]]

    # =======================
    # BIN COUNT
    # =======================
    total_bins <- nrow(dt_region)

    dt_nonzero <- dt_region[total_counts > 0]

    nonzero_bins <- nrow(dt_nonzero)

    # =======================
    # TRANSCRIPTS
    # =======================
    t_stat <- calc_distribution(dt_nonzero$total_counts)

    # =======================
    # GENES
    # =======================
    g_stat <- calc_distribution(dt_nonzero$n_genes)

    res_list[[paste0(sample_name, "_", region_name)]] <- data.table(
      sample = sample_name,
      region = region_name,

      # bin info
      total_bins = total_bins,
      nonzero_bins = nonzero_bins,
      pct_nonzero_bins = nonzero_bins / total_bins,

      # transcripts
      t_min = t_stat$min,
      t_Q1 = t_stat$Q1,
      t_median = t_stat$median,
      t_mean = t_stat$mean,
      t_Q3 = t_stat$Q3,
      t_max = t_stat$max,

      # genes
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

message("\nDone.")
message("Saved to: ", OUT_TSV)

######################
# check if those negative probes are removed..
######################

#/home/woodydrylab/FileShare/20260121_Xenium/Xenium_transcript_distribution_summary.tsv

#mat <- readRDS("/home/woodydrylab/FileShare/20260121_Xenium/output-XETG00360__0079031__BIO2_TISSUE_2__20260115__224443/grid05um_out/counts_bins_by_genes_sparse.rds")
#genes <- colnames(mat)
#non_potri <- genes[!grepl("^Potri", genes)]
#
#length(non_potri)
#head(non_potri)


#####################
# use 50um (same with our previous Cla stereo-seq)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

# =======================
# SETTINGS
# =======================
XENIUM_BASE <- "/home/woodydrylab/FileShare/20260121_Xenium"
GRID_NAME   <- "grid05um_out"

ORIG_BIN_UM   <- 5
TARGET_BIN_UM <- 50

OUTDIR_HIST <- file.path(XENIUM_BASE, paste0("hist_transcripts_per_", TARGET_BIN_UM, "um_bin"))
dir.create(OUTDIR_HIST, showWarnings = FALSE, recursive = TRUE)

OUT_TSV <- file.path(XENIUM_BASE, paste0("Xenium_full_QC_summary_", TARGET_BIN_UM, "um.tsv"))

BOTTOM_COORD_FILE <- "Selection_Bottom_coordinates.csv"
UPPER_COORD_FILE  <- "Selection_Upper_coordinates.csv"

# =======================
# HELPERS
# =======================
guess_xy_cols <- function(dt) {
  nms <- names(dt)
  low <- tolower(nms)

  x_idx <- which(low %in% c("x", "x_center"))
  y_idx <- which(low %in% c("y", "y_center"))

  if (length(x_idx) == 0 || length(y_idx) == 0) {
    stop("Cannot find x/y columns in polygon csv.")
  }

  list(x = nms[x_idx[1]], y = nms[y_idx[1]])
}

read_polygon_csv <- function(path) {
  dt <- fread(path)
  xy <- guess_xy_cols(dt)

  poly <- dt[, .(
    x = as.numeric(get(xy$x)),
    y = as.numeric(get(xy$y))
  )]

  if (nrow(poly) < 3) {
    stop("Polygon file has fewer than 3 points: ", path)
  }

  # close polygon if needed
  if (!(poly$x[1] == poly$x[nrow(poly)] && poly$y[1] == poly$y[nrow(poly)])) {
    poly <- rbind(poly, poly[1])
  }

  poly
}

point_in_polygon_vec <- function(px, py, vx, vy) {
  n <- length(vx)
  inside <- rep(FALSE, length(px))

  j <- n
  for (i in seq_len(n)) {
    xi <- vx[i]; yi <- vy[i]
    xj <- vx[j]; yj <- vy[j]

    intersect <- ((yi > py) != (yj > py)) &
      (px < (xj - xi) * (py - yi) / ((yj - yi) + 1e-30) + xi)

    inside <- xor(inside, intersect)
    j <- i
  }
  inside
}

split_bin_meta_by_selection <- function(bin_meta, sample_dir) {
  bottom_path <- file.path(sample_dir, BOTTOM_COORD_FILE)
  upper_path  <- file.path(sample_dir, UPPER_COORD_FILE)

  if (!file.exists(bottom_path) || !file.exists(upper_path)) {
    return(list(ALL = copy(bin_meta)))
  }

  poly_bottom <- read_polygon_csv(bottom_path)
  poly_upper  <- read_polygon_csv(upper_path)

  in_bottom <- point_in_polygon_vec(
    bin_meta$x_center, bin_meta$y_center,
    poly_bottom$x, poly_bottom$y
  )

  in_upper <- point_in_polygon_vec(
    bin_meta$x_center, bin_meta$y_center,
    poly_upper$x, poly_upper$y
  )

  list(
    Bottom = copy(bin_meta[in_bottom]),
    Upper  = copy(bin_meta[in_upper])
  )
}

calc_distribution <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(list(min=NA_real_, Q1=NA_real_, median=NA_real_,
                mean=NA_real_, Q3=NA_real_, max=NA_real_))
  }

  q <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

  list(
    min    = unname(q[1]),
    Q1     = unname(q[2]),
    median = unname(q[3]),
    mean   = mean(x, na.rm = TRUE),
    Q3     = unname(q[4]),
    max    = unname(q[5])
  )
}

plot_hist <- function(x, outfile, title, xlab_txt) {
  x <- x[is.finite(x)]
  x <- x[x > 0]

  if (length(x) < 10) {
    return(invisible(NULL))
  }

  xmax <- as.numeric(quantile(x, 0.99, na.rm = TRUE))

  png(outfile, width = 1200, height = 900)

  par(
    bg = "black",
    col = "white",
    col.axis = "white",
    col.lab = "white",
    col.main = "white"
  )

  hist(
    x[x <= xmax],
    breaks = 100,
    col = "white",
    border = "white",
    main = title,
    xlab = xlab_txt,
    ylab = "Bin count"
  )

  dev.off()
}

# Aggregate 5 um bins into 50 um super-bins
# Keep the original row_id here so the corresponding rows can be retrieved from the sparse matrix later
aggregate_to_superbins <- function(dt_region, target_bin_um = 50) {
  if (nrow(dt_region) == 0) {
    return(data.table(
      super_x = integer(),
      super_y = integer(),
      total_counts = numeric(),
      n_subbins = integer(),
      x_center = numeric(),
      y_center = numeric(),
      row_ids = list()
    ))
  }

  dt2 <- copy(dt_region)

  # Use the minimum coordinates in this region as the anchor
  x0 <- min(dt2$x_center, na.rm = TRUE)
  y0 <- min(dt2$y_center, na.rm = TRUE)

  dt2[, super_x := floor((x_center - x0) / target_bin_um)]
  dt2[, super_y := floor((y_center - y0) / target_bin_um)]

  agg <- dt2[, .(
    total_counts = sum(total_counts, na.rm = TRUE),
    n_subbins    = .N,
    x_center     = mean(x_center, na.rm = TRUE),
    y_center     = mean(y_center, na.rm = TRUE),
    row_ids      = list(row_id)
  ), by = .(super_x, super_y)]

  agg[]
}

# Correctly calculate the number of unique genes in each super-bin from the sparse matrix
# mat: rows = original 5 um bins, cols = genes
calc_superbin_gene_counts <- function(agg_dt, mat) {
  n <- nrow(agg_dt)

  if (n == 0) {
    return(numeric())
  }

  g_counts <- numeric(n)

  for (i in seq_len(n)) {
    rows <- agg_dt$row_ids[[i]]

    if (length(rows) == 0) {
      g_counts[i] <- 0
      next
    }

    # Submatrix colSums gives the total counts for each gene in this super-bin
    gene_sum <- Matrix::colSums(mat[rows, , drop = FALSE])

    # Count uniquely detected genes
    g_counts[i] <- sum(gene_sum > 0)
  }

  g_counts
}

# =======================
# MAIN
# =======================
sample_dirs <- list.dirs(XENIUM_BASE, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^.*/output-", sample_dirs)]

res_list <- list()

for (sdir in sample_dirs) {
  sample_name <- basename(sdir)
  message("Processing: ", sample_name)

  grid_dir   <- file.path(sdir, GRID_NAME)
  meta_path  <- file.path(grid_dir, "bin_metadata.tsv")
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

  # Check row alignment
  if (nrow(dt) != nrow(mat)) {
    stop(
      "Row mismatch in sample ", sample_name,
      ": nrow(bin_metadata.tsv) = ", nrow(dt),
      ", nrow(counts_bins_by_genes_sparse.rds) = ", nrow(mat)
    )
  }

  # Allow each row in dt to map back to the corresponding sparse matrix row
  dt[, row_id := .I]

  regions <- split_bin_meta_by_selection(dt, sdir)

  for (region_name in names(regions)) {
    dt_region <- regions[[region_name]]

    if (nrow(dt_region) == 0) {
      next
    }

    # Convert 5 um bins to 50 um bins
    dt_super <- aggregate_to_superbins(
      dt_region = dt_region,
      target_bin_um = TARGET_BIN_UM
    )

    if (nrow(dt_super) == 0) {
      next
    }

    # Correctly calculate unique genes per 50 um bin
    dt_super[, n_genes := calc_superbin_gene_counts(.SD, mat)]

    # =======================
    # HISTOGRAM
    # =======================
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

    # =======================
    # SUMMARY
    # =======================
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

      # transcript stats
      t_min = t_stat$min,
      t_Q1 = t_stat$Q1,
      t_median = t_stat$median,
      t_mean = t_stat$mean,
      t_Q3 = t_stat$Q3,
      t_max = t_stat$max,

      # gene stats
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

  message("\nDone.")
  message("Histogram dir: ", OUTDIR_HIST)
  message("Summary TSV: ", OUT_TSV)
}



##############

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

IN_TSV <- "/home/woodydrylab/FileShare/20260121_Xenium/Xenium_full_QC_summary_50um.tsv"
OUT_PNG <- "/home/woodydrylab/FileShare/20260121_Xenium/Xenium_50um_gmean_boxplot.png"

OLD_VALUE <- 25.91
USE_REGION <- NULL

dt <- fread(IN_TSV)
dt <- dt[is.finite(g_mean) & !is.na(g_mean)]

if (!is.null(USE_REGION)) {
  dt <- dt[region == USE_REGION]
}

if (nrow(dt) == 0) stop("No data to plot.")

# Key step: manually control the y-axis range
ymin <- min(dt$g_mean, na.rm = TRUE)
ymax <- max(dt$g_mean, na.rm = TRUE)

# Leave space for the arrow annotation
ymax <- max(ymax, OLD_VALUE + 10)

png(OUT_PNG, width = 800, height = 1200, res = 150)

par(mar = c(6, 5, 4, 2) + 0.1)

boxplot(
  dt$g_mean,
  outline = FALSE,
  names = "",
  ylim = c(0, ymax),   # Important for keeping the annotation visible
  ylab = " #Genes per 50 µm bin",
  main = ""
)

points(
  jitter(rep(1, nrow(dt)), amount = 0.05),
  dt$g_mean,
  pch = 16
)

# Reference line (always shown)
abline(h = OLD_VALUE, lty = 2, lwd = 2, col = "red")

# Arrow annotation (top to bottom)
arrows(
  x0 = 1,
  y0 = OLD_VALUE + 6,
  x1 = 1,
  y1 = OLD_VALUE,
  length = 0.12,
  lwd = 2,
  col = "red"
)

# Text annotation
text(
  x = 1,
  y = OLD_VALUE + 7,
  labels = paste0("Previous mean = ", OLD_VALUE),
  pos = 3,
  col = "red",
  cex = 1
)

dev.off()

cat("Saved to:", OUT_PNG, "\n")


png(OUT_PNG, width = 800, height = 1200, res = 150)
par(
  mar = c(6, 5, 4, 2) + 0.1,
  font.lab = 2,
  cex.lab = 1.8,    # y-axis label
  cex.axis = 1.5,   # axis tick labels
  cex.main = 1.8    # title
)

boxplot(
  dt$g_mean,
  outline = FALSE,
  names = "",
  ylim = c(0, ymax),
  ylab = "#Genes per 50 µm bin",
  main = "",
  lwd = 2           # Thicker box outline
)

# Increase point size
points(
  jitter(rep(1, nrow(dt)), amount = 0.05),
  dt$g_mean,
  pch = 16,
  cex = 1.8         # Point size
)

# Make the reference line thicker
abline(
  h = OLD_VALUE,
  lty = 2,
  lwd = 3,
  col = "red"
)

# Make the arrow thicker and larger
arrows(
  x0 = 1,
  y0 = OLD_VALUE + 6,
  x1 = 1,
  y1 = OLD_VALUE,
  length = 0.15,
  lwd = 3,
  col = "red"
)

# Increase annotation text size
text(
  x = 1,
  y = OLD_VALUE + 7,
  labels = paste0("Previous Stereo-seq = ", OLD_VALUE),
  pos = 3,
  col = "red",
  cex = 1.6        # Annotation text size
)

dev.off()

cat("Saved to:", OUT_PNG, "\n")
