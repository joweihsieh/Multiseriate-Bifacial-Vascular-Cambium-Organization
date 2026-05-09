#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

# ============================================================
# USER SETTINGS (EDIT HERE)
# ============================================================

# (1) Xenium outputs root folder (contains many "output-..." folders)
XENIUM_BASE <- "/home/woodydrylab/FileShare/20260121_Xenium"

# (2) LCM files
LCM_FILES <- c(
  cambium = "/home/chingweilu/Multi_Bilayer_Unifacial/Adjusting_UMAP_neighbor/WeiL_LCM/quantification/Ptr_Cambium_gene_abundances.csv",
  xylem   = "/home/chingweilu/Multi_Bilayer_Unifacial/Adjusting_UMAP_neighbor/WeiL_LCM/quantification/Ptr_Xylem_gene_abundances.csv",
  phloem  = "/home/chingweilu/Multi_Bilayer_Unifacial/Adjusting_UMAP_neighbor/WeiL_LCM/quantification/Ptr_Phloem_gene_abundances.csv"
)

# (3) correlation chunk size (smaller = safer, larger = faster)
CHUNK_SIZE <- 5000

# (4) plotting thresholds
PLOT_LOWER <- 0.10
PLOT_MAXI  <- 0.75

# ============================================================
# FUNCTIONS
# ============================================================

clean_gene_id <- function(x) {
  x <- sub("\\.v4\\.1$", "", x)
  x <- gsub("^Potri\\.", "Potri_", x)
  x
}

clean_filename <- function(x) {
  gsub("[^A-Za-z0-9_\\-]+", "_", x)
}

# chunked cor() per bin
cor_by_bin_cor <- function(mat_bins_genes, ref_vec, chunk_size = 5000, method = "pearson") {
  n <- nrow(mat_bins_genes)
  out <- numeric(n)

  for (s in seq(1, n, by = chunk_size)) {
    e <- min(n, s + chunk_size - 1)
    X <- as.matrix(mat_bins_genes[s:e, , drop = FALSE])  # chunk × genes
    out[s:e] <- apply(
      X, 1,
      function(x) cor(x, ref_vec, method = method, use = "pairwise.complete.obs")
    )
  }
  out
}

# build aligned reference vector (one replicate column)
make_ref <- function(dt, genes_x, col_name) {
  ref_dt <- dt[gene_clean %in% genes_x, .(gene_clean, ref = get(col_name))]
  ref_dt$ref[match(genes_x, ref_dt$gene_clean)]
}

# plot (base), with INSIDE color bar so it won't be cropped
ori_par <- par(no.readonly = TRUE)

output_png_figure <- function(plotting_function,
                              output_figure = FALSE,
                              output_path = "temp.png",
                              output_without_margin = FALSE,
                              ...) {
  if (output_figure) {
    png(output_path, pointsize = 10, res = 300, width = 20, height = 15, units = "cm")
  }
  if (output_without_margin) par(mai = c(0, 0, 0, 0)) else par(mai = ori_par$mai)

  plotting_function(output_without_margin = output_without_margin, ...)

  par(mai = ori_par$mai)
  if (output_figure) dev.off()
}

plot_xenium_with_correlation <- function(bin_meta,
                                         cor_colname,
                                         bin_um = 5,
                                         lower = 0.10,
                                         maxivalue = 0.75,
                                         reverse_y = TRUE,
                                         sorted_order = TRUE,
                                         output_without_margin = FALSE,
                                         ...) {
  stopifnot(all(c("x_center", "y_center") %in% names(bin_meta)))
  stopifnot(cor_colname %in% names(bin_meta))

  x <- bin_meta$x_center
  y <- bin_meta$y_center
  cor_vector <- bin_meta[[cor_colname]]

  keep <- !is.na(x) & !is.na(y) & !is.na(cor_vector)
  x <- x[keep]; y <- y[keep]; cor_vector <- cor_vector[keep]
  if (reverse_y) y <- -y

  # ---- thresholding for color
  cor_for_color <- cor_vector
  cor_for_color[cor_for_color < lower] <- 0

  denom <- (maxivalue - lower)
  if (denom <= 0) stop("maxivalue must be > lower")

  color_index <- (cor_for_color - lower) / denom
  color_index <- pmax(0, pmin(1, color_index))
  color_index <- round(color_index * 500) + 1

  # light -> red (keep your palette choice)
  color_tick <- c("#EEF2F9", "#C44233")
  color_pool <- colorRampPalette(color_tick)(501)

  if (sorted_order) {
    plot_order <- order(cor_for_color)   # low first, high last (high on top)
  } else {
    plot_order <- sample(length(x))
  }

  # ---- canvas: IMPORTANT (asp=1 + half-bin padding)
  half <- bin_um / 2

  par(bg = "black")
  plot(NA, NA,
       xlim = range(x) + c(-half, half),
       ylim = range(y) + c(-half, half),
       xlab = "", ylab = "",
       axes = !output_without_margin, las = 1,
       asp = 1,
       main = ifelse(output_without_margin, "", cor_colname)
  )

  # ---- draw true non-overlapping tiles in data units (um)
  xo <- x[plot_order]
  yo <- y[plot_order]
  ci <- color_index[plot_order]

  rect(xleft   = xo - half,
       ybottom = yo - half,
       xright  = xo + half,
       ytop    = yo + half,
       col     = color_pool[ci],
       border  = NA)

  # ---- color bar inside plot (same logic as before)
  current_xlim <- par()$usr[1:2]
  current_ylim <- par()$usr[3:4]
  xspan <- diff(current_xlim)
  yspan <- diff(current_ylim)

  legend_x_range <- c(current_xlim[1] + 0.55 * xspan, current_xlim[1] + 0.95 * xspan)
  legend_y <- current_ylim[1] + 0.06 * yspan

  legend_x_vector <- seq(legend_x_range[1], legend_x_range[2], length.out = length(color_pool) + 1)

  for (i in seq_along(color_pool)) {
    segments(legend_x_vector[i], legend_y, legend_x_vector[i + 1], legend_y,
             col = color_pool[i], lwd = 12, lend = "butt")
  }

  text(legend_x_range[1], legend_y + 0.03 * yspan, sprintf("%.2f", lower),
       cex = 0.9, col = "white", adj = c(0, 0))
  text(legend_x_range[2], legend_y + 0.03 * yspan, sprintf("%.2f", maxivalue),
       cex = 0.9, col = "white", adj = c(1, 0))
  text(mean(legend_x_range), legend_y + 0.06 * yspan, "correlation (r)",
       cex = 0.9, col = "white")
}

# ============================================================
# LOAD LCM TABLES ONCE
# ============================================================

dt_list <- lapply(LCM_FILES, function(f) {
  dt <- as.data.table(read.csv(f, check.names = FALSE))
  dt[, gene_clean := clean_gene_id(Gene.ID)]
  dt
})

# ============================================================
# DISCOVER ALL XENIUM SAMPLES + ALL GRID FOLDERS
# ============================================================

sample_dirs <- list.dirs(XENIUM_BASE, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^.*/output-", sample_dirs)]

if (length(sample_dirs) == 0) stop("No sample folders found under XENIUM_BASE: ", XENIUM_BASE)

message("Found ", length(sample_dirs), " sample folders.")

for (sdir in sample_dirs) {
  message("\n============================")
  message("Sample: ", sdir)
  message("============================")

  grid_dirs <- list.dirs(sdir, full.names = TRUE, recursive = FALSE)
  grid_dirs <- grid_dirs[grepl("/grid[0-9]+um_out$", grid_dirs)]

  if (length(grid_dirs) == 0) {
    message("  [skip] no grid*um_out folder found.")
    next
  }

  for (gdir in grid_dirs) {
    message("  Grid: ", gdir)

    mat_path <- file.path(gdir, "counts_bins_by_genes_sparse.rds")
    meta_path <- file.path(gdir, "bin_metadata.tsv")

    if (!file.exists(mat_path) || !file.exists(meta_path)) {
      message("    [skip] missing matrix or metadata.")
      next
    }

    # parse bin_um from folder name
    bin_um <- as.numeric(sub(".*grid([0-9]+)um_out$", "\\1", gdir))

    # load
    mat <- readRDS(mat_path)
    bin_meta <- fread(meta_path)
    if (nrow(bin_meta) != nrow(mat)) {
      message("    [skip] row mismatch: bin_meta vs mat")
      next
    }
    genes_x <- colnames(mat)

    # compute correlations and add columns
    for (tissue in names(dt_list)) {
      dt <- dt_list[[tissue]]
      lcm_cols <- setdiff(names(dt), c("Gene.ID", "gene_clean"))
      if (length(lcm_cols) < 1) next

      # each replicate
      for (col_name in lcm_cols) {
        ref <- make_ref(dt, genes_x, col_name)
        keep_g <- !is.na(ref)
        mat2 <- mat[, keep_g, drop = FALSE]
        ref2 <- ref[keep_g]

        if (sd(ref2, na.rm = TRUE) == 0) next

        safe_col <- gsub("[^A-Za-z0-9_]+", "_", col_name)

        rP <- cor_by_bin_cor(mat2, ref2, chunk_size = CHUNK_SIZE, method = "pearson")
        rS <- cor_by_bin_cor(mat2, ref2, chunk_size = CHUNK_SIZE, method = "spearman")

        colP <- paste0("corr_LCM_", tissue, "_", safe_col, "_pearson")
        colS <- paste0("corr_LCM_", tissue, "_", safe_col, "_spearman")
        bin_meta[, (colP) := rP]
        bin_meta[, (colS) := rS]
      }

      # mean reference
      ref_mean_dt <- dt[gene_clean %in% genes_x,
                        .(ref = rowMeans(.SD, na.rm = TRUE)),
                        by = gene_clean, .SDcols = lcm_cols]
      ref_mean <- ref_mean_dt$ref[match(genes_x, ref_mean_dt$gene_clean)]
      keep_g <- !is.na(ref_mean)
      mat2 <- mat[, keep_g, drop = FALSE]
      ref2 <- ref_mean[keep_g]

      if (sd(ref2, na.rm = TRUE) != 0) {
        rP_mean <- cor_by_bin_cor(mat2, ref2, chunk_size = CHUNK_SIZE, method = "pearson")
        rS_mean <- cor_by_bin_cor(mat2, ref2, chunk_size = CHUNK_SIZE, method = "spearman")

        meanP <- paste0("corr_LCM_", tissue, "_MEAN_pearson")
        meanS <- paste0("corr_LCM_", tissue, "_MEAN_spearman")

        bin_meta[, (meanP) := rP_mean]
        bin_meta[, (meanS) := rS_mean]
      }
    }

    # write bin_metadata_LCM.tsv
    out_meta <- file.path(gdir, "bin_metadata_LCM.tsv")
    fwrite(bin_meta, out_meta, sep = "\t")
    message("    wrote: ", out_meta)

    # plotting
    out_plot_dir <- file.path(gdir, "corr_maps_spearman_pearson")
    dir.create(out_plot_dir, showWarnings = FALSE, recursive = TRUE)

    cols_to_plot <- grep("^corr_LCM_(cambium|xylem|phloem)_", names(bin_meta), value = TRUE)
    if (length(cols_to_plot) == 0) {
      message("    [skip] no corr_LCM_* columns to plot.")
      next
    }

    for (cc in cols_to_plot) {
      outpng <- file.path(out_plot_dir, paste0(clean_filename(cc), ".png"))
      output_png_figure(
        plotting_function = plot_xenium_with_correlation,
        output_figure = TRUE,
        output_path = outpng,
        output_without_margin = TRUE,
        bin_meta = bin_meta,
        cor_colname = cc,
        bin_um = bin_um,
        lower = PLOT_LOWER,
        maxivalue = PLOT_MAXI,
        reverse_y = TRUE,
        sorted_order = TRUE
      )
    }

    message("    plots: ", out_plot_dir)
  }
}

message("\nAll done.")
