setwd("/Users/jo-wei/Library/CloudStorage/Dropbox/YCL/Spatial_Transcriptomes/Xenium/20260121_Data/simulation/20260406")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

set.seed(123)

############################################################
# settings
############################################################
OUTDIR <- "segregate_gradient_26x9_5um"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# each grid = 5 um
GRID_UM <- 5

# data size
NX <- 26   # width  = 26 grids
NY <- 9    # height = 9 grids

# each cell layer = 2 grids wide
LAYER_WIDTH <- 2
stopifnot(NX %% LAYER_WIDTH == 0)

N_LAYER <- NX / LAYER_WIDTH   # 13 layers total

# variation control
# within-layer should be smaller than between-layer
WITHIN_LAYER_SD  <- 0.02
WITHIN_COLUMN_SD <- 0.005

# colors: left -> right only
COL_LEFT  <- "#B3A250"
COL_RIGHT <- "#F2B700"

############################################################
# helper
############################################################
make_grid <- function(nx = NX, ny = NY, grid_um = GRID_UM) {
  dt <- CJ(ix = 1:nx, iy = 1:ny)

  dt[, x_idx := ix]
  dt[, y_idx := iy]

  # physical coordinates in um (center of each square)
  dt[, x_um := (ix - 0.5) * grid_um]
  dt[, y_um := (iy - 0.5) * grid_um]

  dt
}

make_linear_palette <- function(n = 256,
                                col_left = COL_LEFT,
                                col_right = COL_RIGHT) {
  scales::gradient_n_pal(
    colours = c(col_left, col_right),
    values = c(0, 1),
    space = "Lab"
  )(seq(0, 1, length.out = n))
}

############################################################
# build rightward layer-based gradient
############################################################
make_segregate_pattern_layer <- function(nx = NX,
                                   ny = NY,
                                   layer_width = LAYER_WIDTH,
                                   within_layer_sd = WITHIN_LAYER_SD,
                                   within_column_sd = WITHIN_COLUMN_SD) {

  dt <- make_grid(nx, ny)

  # assign 2-grid-wide layers
  dt[, layer_id := ceiling(ix / layer_width)]

  n_layer <- nx / layer_width
  layer_dt <- data.table(layer_id = 1:n_layer)

  # single-direction linear gradient from left to right
  # leftmost = 0, rightmost = 1
  layer_dt[, layer_mean := seq(from = 0, to = 1, length.out = n_layer)]
  layer_dt[, region := "rightward"]

  dt <- merge(dt, layer_dt, by = "layer_id", all.x = TRUE)

  # tiny column offset within each layer
  col_dt <- unique(dt[, .(ix, layer_id)])
  col_dt[, column_offset := rnorm(.N, mean = 0, sd = within_column_sd)]

  dt <- merge(dt, col_dt, by = c("ix", "layer_id"), all.x = TRUE)

  # final value:
  # main structure = layer_mean
  # small within-layer variation
  # very small within-column offset
  dt[, value := layer_mean + column_offset + rnorm(.N, mean = 0, sd = within_layer_sd)]

  # clamp to [0, 1]
  dt[, value := pmax(0, pmin(1, value))]

  dt[, region := factor(region, levels = c("rightward"))]
  dt[, layer_label := paste0("L", sprintf("%02d", layer_id))]

  dt[]
}

############################################################
# summary check
############################################################
calc_variation_summary <- function(dt) {

  layer_stats <- dt[, .(
    layer_mean_obs = mean(value),
    layer_sd_obs   = sd(value)
  ), by = layer_id]

  col_stats <- dt[, .(
    col_mean_obs = mean(value),
    col_sd_obs   = sd(value)
  ), by = ix]

  list(
    layer_stats = layer_stats,
    col_stats = col_stats,
    mean_within_layer_sd = mean(layer_stats$layer_sd_obs),
    sd_across_layer_means = sd(layer_stats$layer_mean_obs),
    mean_within_column_sd = mean(col_stats$col_sd_obs),
    sd_across_column_means = sd(col_stats$col_mean_obs)
  )
}

############################################################
# plots
############################################################
plot_pattern <- function(dt, outfile) {
  pal <- make_linear_palette(256)

  p <- ggplot(dt, aes(x = x_idx, y = y_idx, fill = value)) +
    geom_tile(color = NA) +
    geom_vline(xintercept = seq(1.5, NX, by = 1), color = "#D0D0D0", linewidth = 0.2) +
    geom_hline(yintercept = seq(1.5, NY, by = 1), color = "#D0D0D0", linewidth = 0.2) +

    coord_fixed() +
    scale_y_reverse() +
    scale_fill_gradientn(
      colours = pal,
      limits = c(0, 1),
      oob = squish
    ) +
    scale_x_continuous(
      breaks = 1:NX,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = 1:NY,
      expand = c(0, 0)
    ) +
    theme_bw(base_size = 14) +
    labs(
      title = "segregate pattern (rightward gradient)",
      subtitle = "26 × 9 grids; each layer is 2 grids wide; one-direction linear gradient",
      x = "Grid index (x)",
      y = "Grid index (y)",
      fill = "Value"
    )

  ggsave(outfile, p, width = 10, height = 4)
}

plot_pattern_um <- function(dt, outfile) {
  pal <- make_linear_palette(256)

  p <- ggplot(dt, aes(x = x_um, y = y_um, fill = value)) +
    geom_tile(width = GRID_UM, height = GRID_UM, color = NA) +
    coord_fixed() +
    scale_y_reverse() +
    scale_fill_gradientn(
      colours = pal,
      limits = c(0, 1),
      oob = squish
    ) +
    theme_bw(base_size = 14) +
    labs(
      title = "Rightward gradient in physical scale",
      subtitle = paste0("Each square = ", GRID_UM, " µm; total size = ", NX * GRID_UM, " × ", NY * GRID_UM, " µm"),
      x = "x (µm)",
      y = "y (µm)",
      fill = "Value"
    )

  ggsave(outfile, p, width = 10, height = 4)
}

plot_layer_means <- function(dt, outfile) {
  sum_dt <- dt[, .(
    mean_value = mean(value),
    sd_value   = sd(value)
  ), by = layer_id]

  p <- ggplot(sum_dt, aes(x = layer_id, y = mean_value)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    theme_bw(base_size = 14) +
    scale_x_continuous(breaks = 1:N_LAYER) +
    labs(
      title = "Mean value of each 2-grid layer",
      x = "Layer ID",
      y = "Observed layer mean"
    )

  ggsave(outfile, p, width = 8, height = 4)
}

plot_column_means <- function(dt, outfile) {
  sum_dt <- dt[, .(
    mean_value = mean(value),
    sd_value   = sd(value)
  ), by = ix]

  p <- ggplot(sum_dt, aes(x = ix, y = mean_value)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    theme_bw(base_size = 14) +
    scale_x_continuous(breaks = 1:NX) +
    labs(
      title = "Mean value of each column",
      x = "Column",
      y = "Observed column mean"
    )

  ggsave(outfile, p, width = 10, height = 4)
}

plot_variation_check <- function(dt, outfile) {
  layer_stats <- dt[, .(
    mean_value = mean(value),
    sd_value   = sd(value)
  ), by = layer_id]

  stat_dt <- data.table(
    metric = c(
      "Mean within-layer SD",
      "SD across layer means"
    ),
    value = c(
      mean(layer_stats$sd_value),
      sd(layer_stats$mean_value)
    )
  )

  p <- ggplot(stat_dt, aes(x = metric, y = value)) +
    geom_col(width = 0.65) +
    theme_bw(base_size = 14) +
    labs(
      title = "Variation check at layer level",
      x = NULL,
      y = "Value"
    )

  ggsave(outfile, p, width = 7, height = 4)
}

############################################################
# palette check
############################################################
check_palette <- function(n = N_LAYER,
                          outfile = file.path(OUTDIR, "palette_check_13_layers.png")) {
  pal <- make_linear_palette(n)

  dt <- data.table(
    idx = 1:n,
    color = pal
  )

  p <- ggplot(dt, aes(x = idx, y = 1, fill = color)) +
    geom_tile(height = 0.8) +
    geom_text(aes(label = color), angle = 90, size = 3, vjust = 0.5, hjust = 0) +
    scale_fill_identity() +
    scale_x_continuous(breaks = 1:n, expand = c(0.01, 0.01)) +
    theme_void(base_size = 12) +
    labs(title = paste0("Rightward color gradient check (", n, " layer steps)")) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(10, 10, 40, 10)
    )

  ggsave(outfile, p, width = 10, height = 3)
  print(p)

  invisible(dt)
}

check_palette_with_hex <- function(n = N_LAYER) {
  pal <- make_linear_palette(n)
  cat("Palette hex codes:\n")
  print(pal)
  invisible(pal)
}

############################################################
# main
############################################################
dt_segregate <- make_segregate_pattern_layer()

# save raw table
fwrite(dt_segregate, file.path(OUTDIR, "segregate_pattern_26x9_rightward.tsv"), sep = "\t")

# save plots
plot_pattern(
  dt_segregate,
  file.path(OUTDIR, "segregate_pattern_heatmap_grid_index.png")
)

plot_pattern_um(
  dt_segregate,
  file.path(OUTDIR, "segregate_pattern_heatmap_um.png")
)

plot_layer_means(
  dt_segregate,
  file.path(OUTDIR, "segregate_layer_means.png")
)

plot_column_means(
  dt_segregate,
  file.path(OUTDIR, "segregate_column_means.png")
)

plot_variation_check(
  dt_segregate,
  file.path(OUTDIR, "segregate_variation_check_layer.png")
)

# save variation summary
var_sum <- calc_variation_summary(dt_segregate)

summary_dt <- data.table(
  metric = c(
    "mean_within_layer_sd",
    "sd_across_layer_means",
    "mean_within_column_sd",
    "sd_across_column_means"
  ),
  value = c(
    var_sum$mean_within_layer_sd,
    var_sum$sd_across_layer_means,
    var_sum$mean_within_column_sd,
    var_sum$sd_across_column_means
  )
)

fwrite(summary_dt, file.path(OUTDIR, "segregate_variation_summary.tsv"), sep = "\t")
fwrite(var_sum$layer_stats, file.path(OUTDIR, "segregate_layer_stats.tsv"), sep = "\t")
fwrite(var_sum$col_stats, file.path(OUTDIR, "segregate_column_stats.tsv"), sep = "\t")

# palette check
palette_dt  <- check_palette()
palette_hex <- check_palette_with_hex()

cat("Done.\n")
cat("Output directory:", OUTDIR, "\n")
cat("Grid size (um):", GRID_UM, "\n")
cat("Total size (um):", NX * GRID_UM, "x", NY * GRID_UM, "\n")
cat("N_LAYER:", N_LAYER, "\n")
cat("Mean within-layer SD:", var_sum$mean_within_layer_sd, "\n")
cat("SD across layer means:", var_sum$sd_across_layer_means, "\n")
cat("Mean within-column SD:", var_sum$mean_within_column_sd, "\n")
cat("SD across column means:", var_sum$sd_across_column_means, "\n")



############################################################
# kmeans on single value per bin
############################################################
K_VALUES <- 2:6

run_kmeans_on_value <- function(dt, k, seed = 123) {
  set.seed(seed)

  km <- kmeans(dt$value, centers = k, nstart = 50)

  centers_dt <- data.table(
    cluster_raw = 1:k,
    center_value = as.numeric(km$centers[, 1])
  )[order(center_value)]

  # reorder clusters from low -> high value
  centers_dt[, cluster := seq_len(.N)]

  map_dt <- centers_dt[, .(cluster_raw, cluster)]

  out_dt <- copy(dt)
  out_dt[, cluster_raw := km$cluster]
  out_dt <- merge(out_dt, map_dt, by = "cluster_raw", all.x = TRUE, sort = FALSE)

  setorder(out_dt, y_idx, x_idx)

  list(
    dt = out_dt,
    centers = centers_dt
  )
}

plot_kmeans_heatmap <- function(dt, k, outfile) {
  p <- ggplot(dt, aes(x = x_idx, y = y_idx, fill = factor(cluster))) +
    geom_tile(color = "white", linewidth = 0.2) +
    coord_fixed() +
    scale_y_reverse() +
    scale_x_continuous(
      breaks = 1:NX,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = 1:NY,
      expand = c(0, 0)
    ) +
    theme_bw(base_size = 14) +
    labs(
      title = paste0("K-means on segregate pattern (K = ", k, ")"),
      subtitle = "Clusters reordered from low to high center value",
      x = "Grid index (x)",
      y = "Grid index (y)",
      fill = "Cluster"
    )

  ggsave(outfile, p, width = 10, height = 4)
}

plot_kmeans_heatmap_um <- function(dt, k, outfile) {
  p <- ggplot(dt, aes(x = x_um, y = y_um, fill = factor(cluster))) +
    geom_tile(width = GRID_UM, height = GRID_UM, color = "white", linewidth = 0.2) +
    coord_fixed() +
    scale_y_reverse() +
    theme_bw(base_size = 14) +
    labs(
      title = paste0("K-means on segregate pattern in physical scale (K = ", k, ")"),
      subtitle = paste0("Each square = ", GRID_UM, " µm"),
      x = "x (µm)",
      y = "y (µm)",
      fill = "Cluster"
    )

  ggsave(outfile, p, width = 10, height = 4)
}

plot_layer_cluster_fraction <- function(dt, k, outfile) {
  comp_dt <- dt[, .N, by = .(layer_id, cluster)]
  comp_dt[, frac := N / sum(N), by = layer_id]

  p <- ggplot(comp_dt, aes(x = layer_id, y = frac, fill = factor(cluster))) +
    geom_col(width = 0.8) +
    theme_bw(base_size = 14) +
    scale_x_continuous(breaks = 1:N_LAYER) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = paste0("Layer composition by cluster (K = ", k, ")"),
      x = "Layer ID",
      y = "Fraction",
      fill = "Cluster"
    )

  ggsave(outfile, p, width = 8, height = 4)
}

plot_cluster_centers <- function(centers_dt, k, outfile) {
  p <- ggplot(centers_dt, aes(x = cluster, y = center_value)) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    theme_bw(base_size = 14) +
    scale_x_continuous(breaks = centers_dt$cluster) +
    labs(
      title = paste0("Cluster centers (K = ", k, ")"),
      x = "Cluster (low -> high)",
      y = "Center value"
    )

  ggsave(outfile, p, width = 6, height = 4)
}

############################################################
# run all K
############################################################
KMEANS_OUTDIR <- file.path(OUTDIR, "kmeans_on_single_value")
dir.create(KMEANS_OUTDIR, showWarnings = FALSE, recursive = TRUE)

kmeans_summary_list <- list()

for (k in K_VALUES) {

  message("Running kmeans for K = ", k)

  kres <- run_kmeans_on_value(dt_segregate, k = k, seed = 123 + k)
  dt_k <- kres$dt
  centers_k <- kres$centers

  # save per-bin cluster table
  fwrite(
    dt_k,
    file.path(KMEANS_OUTDIR, paste0("segregate_pattern_26x9_k", k, ".tsv")),
    sep = "\t"
  )

  # save centers
  fwrite(
    centers_k,
    file.path(KMEANS_OUTDIR, paste0("kmeans_centers_k", k, ".tsv")),
    sep = "\t"
  )

  # layer x cluster counts
  layer_cluster_dt <- dt_k[, .N, by = .(layer_id, cluster)][order(layer_id, cluster)]
  layer_cluster_dt[, fraction := N / sum(N), by = layer_id]

  fwrite(
    layer_cluster_dt,
    file.path(KMEANS_OUTDIR, paste0("layer_cluster_composition_k", k, ".tsv")),
    sep = "\t"
  )

  # plots
  plot_kmeans_heatmap(
    dt_k,
    k,
    file.path(KMEANS_OUTDIR, paste0("kmeans_heatmap_grid_k", k, ".png"))
  )

  plot_kmeans_heatmap_um(
    dt_k,
    k,
    file.path(KMEANS_OUTDIR, paste0("kmeans_heatmap_um_k", k, ".png"))
  )

  plot_layer_cluster_fraction(
    dt_k,
    k,
    file.path(KMEANS_OUTDIR, paste0("layer_cluster_fraction_k", k, ".png"))
  )

  plot_cluster_centers(
    centers_k,
    k,
    file.path(KMEANS_OUTDIR, paste0("cluster_centers_k", k, ".png"))
  )

  # summary
  kmeans_summary_list[[as.character(k)]] <- data.table(
    K = k,
    n_cluster = nrow(centers_k),
    min_center = min(centers_k$center_value),
    max_center = max(centers_k$center_value)
  )
}

kmeans_summary_dt <- rbindlist(kmeans_summary_list)
fwrite(
  kmeans_summary_dt,
  file.path(KMEANS_OUTDIR, "kmeans_summary.tsv"),
  sep = "\t"
)

cat("K-means done.\n")
cat("K-means output directory:", KMEANS_OUTDIR, "\n")


############################################################
# Radial distance calculation & plots (Density + Boxplot)
############################################################
############################################################
# Radial distance calculation & plots (Histogram + Boxplot)
############################################################
plot_radial_distance_dist <- function(dt, k, center_x, center_y, outdir) {
  # Copy the data to avoid modifying the original data.table
  plot_dt <- copy(dt)
  
  # Calculate the distance r to the virtual circle center, in um
  plot_dt[, dx := x_um - center_x]
  plot_dt[, dy := y_um - center_y]
  plot_dt[, r := sqrt(dx^2 + dy^2)]
  #plot_dt[, r := abs(dx)]

  # ==========================================
  # 1. Plot a histogram with binwidth = 5
  # ==========================================
  p_hist <- ggplot(plot_dt, aes(x = r, fill = factor(cluster))) +
    # Use position = "identity" so bars from different clusters overlap; alpha makes the overlap easier to inspect
    geom_histogram(binwidth = 5, position = "identity", alpha = 0.6, color = "white", linewidth = 0.2) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0("Radial Distance Histogram (K = ", k, ")"),
      subtitle = paste0("Center at X = ", center_x, " µm, Y = ", center_y, " µm | Bin width = 5 µm"),
      x = "Distance to right-side center r (µm)",
      y = "Count (Number of grids)",
      fill = "Cluster"
    )
  
  # Use histogram in the output filename
  ggsave(file.path(outdir, paste0("radial_distance_histogram_k", k, ".png")), p_hist, width = 7, height = 4.5)
  
  # ==========================================
  # 2. Plot a boxplot, retained for comparison
  # ==========================================
  p_box <- ggplot(plot_dt, aes(x = factor(cluster), y = r, fill = factor(cluster))) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0("Radial Distance Boxplot (K = ", k, ")"),
      subtitle = paste0("Center at X = ", center_x, " µm, Y = ", center_y, " µm"),
      x = "Cluster",
      y = "Distance to right-side center r (µm)",
      fill = "Cluster"
    )
  ggsave(file.path(outdir, paste0("radial_distance_boxplot_k", k, ".png")), p_box, width = 7, height = 4.5)
  
  # Return the table with calculated r values for downstream export
  return(plot_dt)
}

############################################################
# run all K
############################################################
KMEANS_OUTDIR <- file.path(OUTDIR, "kmeans_on_single_value")
dir.create(KMEANS_OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Set the position of the virtual circle center on the right side
RIGHT_OFFSET_UM <- 50 
CENTER_X_UM <- (NX * GRID_UM) + RIGHT_OFFSET_UM
CENTER_Y_UM <- (NY * GRID_UM) / 2 

cat(sprintf("Setting right-side center for radial density at X=%.1f, Y=%.1f\n", CENTER_X_UM, CENTER_Y_UM))

kmeans_summary_list <- list()

for (k in K_VALUES) {

  message("Running kmeans for K = ", k)

  # Use dt_segregate
  kres <- run_kmeans_on_value(dt_segregate, k = k, seed = 123 + k)
  dt_k <- kres$dt
  centers_k <- kres$centers

  # Save the per-bin cluster table, using segregate in the filename
  fwrite(
    dt_k,
    file.path(KMEANS_OUTDIR, paste0("segregate_pattern_26x9_k", k, ".tsv")),
    sep = "\t"
  )

  # save centers
  fwrite(
    centers_k,
    file.path(KMEANS_OUTDIR, paste0("kmeans_centers_k", k, ".tsv")),
    sep = "\t"
  )

  # layer x cluster counts
  layer_cluster_dt <- dt_k[, .N, by = .(layer_id, cluster)][order(layer_id, cluster)]
  layer_cluster_dt[, fraction := N / sum(N), by = layer_id]

  fwrite(
    layer_cluster_dt,
    file.path(KMEANS_OUTDIR, paste0("layer_cluster_composition_k", k, ".tsv")),
    sep = "\t"
  )

  # plots
  plot_kmeans_heatmap(
    dt_k,
    k,
    file.path(KMEANS_OUTDIR, paste0("kmeans_heatmap_grid_k", k, ".png"))
  )

  plot_kmeans_heatmap_um(
    dt_k,
    k,
    file.path(KMEANS_OUTDIR, paste0("kmeans_heatmap_um_k", k, ".png"))
  )

  plot_layer_cluster_fraction(
    dt_k,
    k,
    file.path(KMEANS_OUTDIR, paste0("layer_cluster_fraction_k", k, ".png"))
  )

  plot_cluster_centers(
    centers_k,
    k,
    file.path(KMEANS_OUTDIR, paste0("cluster_centers_k", k, ".png"))
  )

  # --------------------------------------------------------
  # Calculate radial distance and generate the histogram and boxplot
  # --------------------------------------------------------
  dt_with_r <- plot_radial_distance_dist(
    dt = dt_k,
    k = k,
    center_x = CENTER_X_UM,
    center_y = CENTER_Y_UM,
    outdir = KMEANS_OUTDIR
  )
  
  # Save the data table containing r-distance values
  fwrite(
    dt_with_r[, .(x_idx, y_idx, x_um, y_um, value, cluster, dx, dy, r)],
    file.path(KMEANS_OUTDIR, paste0("radial_distance_data_k", k, ".tsv")),
    sep = "\t"
  )
  # --------------------------------------------------------

  # summary
  kmeans_summary_list[[as.character(k)]] <- data.table(
    K = k,
    n_cluster = nrow(centers_k),
    min_center = min(centers_k$center_value),
    max_center = max(centers_k$center_value)
  )
} 

kmeans_summary_dt <- rbindlist(kmeans_summary_list)
fwrite(
  kmeans_summary_dt,
  file.path(KMEANS_OUTDIR, "kmeans_summary.tsv"),
  sep = "\t"
)

cat("K-means done.\n")
cat("K-means output directory:", KMEANS_OUTDIR, "\n")