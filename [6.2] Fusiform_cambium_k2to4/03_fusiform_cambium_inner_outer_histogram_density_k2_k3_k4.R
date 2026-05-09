#!/usr/bin/env Rscript
# run at 20260406
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

############################################################
# user settings
############################################################
BASE_ROOT <- "/home/woodydrylab/FileShare/20260121_Xenium"

# Keep only directories matching the subset_from_k10_c01-style format
SUBSET_PATTERN <- "^subset_from_k[0-9]+_c[0-9]{2}$"

# Analyze only k = 2 and k = 3
K_VALUES <- c(2, 3, 4)
#K_VALUES <- c(4)


# angle bins
N_ANGLE_BINS <- 180

# Additional density plot using only r >= cutoff
R_CUTOFF_FOR_EXTRA_DENSITY <- 1000

PAIR_COLORS <- c(
  "1" = "#F2B700",
  "2" = "#377EB8",
  "3" = "#9467BD",
  "4" = "#F8766D"
)

OUT_SUMMARY_ALL <- file.path(BASE_ROOT, "all_pairwise_radial_ordering_summary.tsv")

############################################################
# helpers
############################################################
get_pair_info <- function(k_value) {
  if (k_value == 2) {
    list(
      target_clusters = c(1L, 2L),
      pairs = list(c(1L, 2L))
    )
  } else if (k_value == 3) {
    list(
      target_clusters = c(1L, 2L, 3L),
      pairs = list(c(1L, 2L), c(1L, 3L), c(2L, 3L))
    )
  } else if (k_value == 4) {
    list(
      target_clusters = c(1L, 2L, 3L, 4L),
      pairs = list(
        c(1L, 2L), c(1L, 3L), c(1L, 4L),
        c(2L, 3L), c(2L, 4L), c(3L, 4L)
      )
    )
  } else {
    stop("Only k = 2, 3, or 4 is supported.")
  }
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

  # rounded x/y fallback
  if (any(is.na(out$stem))) {
    message("  Some bins not matched in first pass; trying rounded x/y fallback...")

    sub_xy2   <- find_xy_cols(sub_dt)
    stem1_xy2 <- find_xy_cols(stem1_dt)
    stem2_xy2 <- find_xy_cols(stem2_dt)

    sub2 <- copy(sub_dt)
    s12  <- copy(stem1_dt)
    s22  <- copy(stem2_dt)

    sub2[, rx := round(get(sub_xy2$x), 4)]
    sub2[, ry := round(get(sub_xy2$y), 4)]

    s12[, rx := round(get(stem1_xy2$x), 4)]
    s12[, ry := round(get(stem1_xy2$y), 4)]
    s12 <- unique(s12[, .(rx, ry, stem = "stem1")])

    s22[, rx := round(get(stem2_xy2$x), 4)]
    s22[, ry := round(get(stem2_xy2$y), 4)]
    s22 <- unique(s22[, .(rx, ry, stem = "stem2")])

    stem_map2 <- unique(rbindlist(list(s12, s22)))

    out2 <- copy(sub_dt)
    out2[, rx := round(get(sub_xy2$x), 4)]
    out2[, ry := round(get(sub_xy2$y), 4)]
    out2 <- merge(out2, stem_map2, by = c("rx", "ry"), all.x = TRUE)

    miss_idx <- which(is.na(out$stem))
    if (length(miss_idx) > 0 && "stem" %in% names(out2)) {
      out$stem[miss_idx] <- out2$stem[miss_idx]
    }
  }

  out[]
}

safe_density_plot <- function(df, colors, title_text, outfile, cutoff = NULL) {
  if (nrow(df) == 0) return(invisible(NULL))

  p <- ggplot(df, aes(x = r, color = factor(cluster_raw), fill = factor(cluster_raw))) +
    geom_density(alpha = 0.25, linewidth = 1) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_classic() +
    labs(
      title = title_text,
      x = "r (distance to center)",
      color = "Cluster",
      fill = "Cluster"
    )

  if (!is.null(cutoff)) {
    p <- p + xlim(cutoff, NA)
  }

  ggsave(outfile, p, width = 6, height = 4.8, dpi = 600)
}


safe_density_plot_split_fixed <- function(df, colors, title_text, out_prefix, cutoff = NULL) {

  clusters <- sort(unique(df$cluster_raw))

  # Fix the x-axis range using all data
  x_min <- min(df$r, na.rm = TRUE)
  x_max <- max(df$r, na.rm = TRUE)

  # Compute the global density first to obtain the maximum height
  dens_all <- density(df$r)
  y_max <- max(dens_all$y, na.rm = TRUE)

  for (cl in clusters) {

    df_cl <- df[cluster_raw == cl]
    if (nrow(df_cl) == 0) next

    p <- ggplot(df_cl, aes(x = r)) +
      geom_density(
        color = colors[as.character(cl)],
        fill  = colors[as.character(cl)],
        alpha = 0.3,
        linewidth = 1.2
      ) +
      theme_classic() +
      labs(
        title = paste0(title_text, " | cluster ", cl),
        x = "r (distance to center)",
        y = "Density"
      ) +
      coord_cartesian(xlim = c(x_min, x_max), ylim = c(0, y_max))

    if (!is.null(cutoff)) {
      p <- p + coord_cartesian(xlim = c(cutoff, x_max), ylim = c(0, y_max))
    }

    outfile <- paste0(out_prefix, "_cluster", cl, ".png")

    ggsave(outfile, p, width = 6, height = 4.8, dpi = 600)

    message("Saved: ", outfile)
  }
}

safe_histogram_plot <- function(df, colors, title_text, outfile, cutoff = NULL) {
  if (nrow(df) == 0) return(invisible(NULL))

  # Use binwidth = 5 and overlay semi-transparent histograms
  p <- ggplot(df, aes(x = r, fill = factor(cluster_raw))) +
    geom_histogram(binwidth = 5, position = "identity", alpha = 0.6, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = colors) +
    theme_classic() +
    labs(
      title = title_text,
      x = "r (distance to center)",
      y = "Count (Number of bins/cells)",
      fill = "Cluster"
    )

  if (!is.null(cutoff)) {
    p <- p + xlim(cutoff, NA)
  }

  ggsave(outfile, p, width = 6, height = 4.8, dpi = 600)
}

extract_context_from_path <- function(k_dir) {
  list(
    sample_id  = basename(dirname(dirname(dirname(k_dir)))),
    grid_dir   = basename(dirname(dirname(k_dir))),
    subset_dir = basename(dirname(k_dir)),
    k_dir      = k_dir
  )
}

analyze_one_df <- function(df, k_dir, k_value, unit_name = "whole", out_prefix = NULL) {
  required_cols <- c("x_center", "y_center", "cluster_raw")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    cat("[skip] missing columns:", paste(missing_cols, collapse = ", "), "\n")
    return(NULL)
  }

  df <- copy(df)
  df[, cluster_raw := as.integer(cluster_raw)]
  df <- df[!is.na(cluster_raw)]

  pair_info <- get_pair_info(k_value)
  target_clusters <- pair_info$target_clusters
  pairs <- pair_info$pairs

  df <- df[cluster_raw %in% target_clusters]

  cat("Cluster counts for", unit_name, ":\n")
  print(df[, .N, by = cluster_raw][order(cluster_raw)])

  if (nrow(df) == 0) {
    cat("[skip] no target clusters found for", unit_name, "\n")
    return(NULL)
  }

  center_x <- mean(df$x_center)
  center_y <- mean(df$y_center)

  df[, dx := x_center - center_x]
  df[, dy := y_center - center_y]
  df[, theta := atan2(dy, dx)]
  df[, r := sqrt(dx^2 + dy^2)]

  breaks <- seq(-pi, pi, length.out = N_ANGLE_BINS + 1)
  df[, angle_bin := cut(theta, breaks = breaks, include.lowest = TRUE)]

  angle_stats <- df[, .(
    n = .N,
    med_r = median(r)
  ), by = .(angle_bin, cluster_raw)]

  vote_list <- rbindlist(lapply(pairs, function(p) {
    a <- p[1]
    b <- p[2]

    w <- angle_stats[cluster_raw %in% c(a, b)]

    bins_ok <- w[, .(k = uniqueN(cluster_raw)), by = angle_bin][k == 2, angle_bin]
    w <- w[angle_bin %in% bins_ok]

    if (nrow(w) == 0) return(data.table())

    wide <- dcast(w, angle_bin ~ cluster_raw, value.var = "med_r")
    setnames(wide, c("angle_bin", "ra", "rb"))

    wide[, `:=`(a = a, b = b)]
    wide[, winner := ifelse(ra < rb, a, b)]
    wide[, margin := abs(ra - rb)]
    wide
  }), fill = TRUE)

  if (nrow(vote_list) == 0) {
    cat("[skip] no usable angle bins for", unit_name, "\n")
    return(NULL)
  }

  pair_summary <- vote_list[, .(
    n_bins = .N,
    frac_a_inner = mean(winner == a),
    frac_b_inner = mean(winner == b),
    winner_majority = ifelse(mean(winner == a) > 0.5, a, b),
    median_margin = median(margin)
  ), by = .(a, b)]

  pair_summary[, pair := paste0(a, " vs ", b)]
  pair_summary[, k := k_value]

  ctx <- extract_context_from_path(k_dir)
  pair_summary[, sample_id := ctx$sample_id]
  pair_summary[, grid_dir := ctx$grid_dir]
  pair_summary[, subset_dir := ctx$subset_dir]
  pair_summary[, k_dir := ctx$k_dir]
  pair_summary[, unit := unit_name]

  print(pair_summary)

  prefix <- if (is.null(out_prefix)) "" else paste0(out_prefix, "_")

  fwrite(
    pair_summary,
    file.path(k_dir, paste0(prefix, "pairwise_radial_ordering_summary.tsv")),
    sep = "\t"
  )

  cat("\nMajority-inner call for", unit_name, ":\n")
  for (ii in seq_len(nrow(pair_summary))) {
    cat(
      "k =", k_value, "|",
      pair_summary$pair[ii], "-> cluster",
      pair_summary$winner_majority[ii],
      "is more inner (frac left inner =",
      round(pair_summary$frac_a_inner[ii], 3), ")\n"
    )
  }
  cat("\n")

  p1 <- ggplot(pair_summary, aes(x = pair, y = frac_a_inner)) +
    geom_col() +
    theme_classic() +
    labs(
      title = paste0("k = ", k_value, " | ", unit_name, " | pairwise radial ordering"),
      x = "Pair",
      y = "Fraction where left cluster is more inner (smaller median r)"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(
    file.path(k_dir, paste0(prefix, "pairwise_radial_ordering_bar.png")),
    p1, width = 7, height = 4.5, dpi = 600
  )

  cluster_ids_present <- sort(unique(df$cluster_raw))
  cluster_colors <- PAIR_COLORS[as.character(cluster_ids_present)]
  cluster_colors <- cluster_colors[!is.na(cluster_colors)]

  # ==========================================
  # 1. Export density plots and histograms for the full range
  # ==========================================
  safe_density_plot(
    df = df,
    colors = cluster_colors,
    title_text = paste0("k = ", k_value, " | ", unit_name, " | radial r density"),
    outfile = file.path(k_dir, paste0(prefix, "radial_r_density.png"))
  )

  safe_density_plot_split_fixed(
    df = df,
    colors = cluster_colors,
    title_text = paste0("k = ", k_value, " | ", unit_name, " | radial r density"),
    out_prefix = file.path(k_dir, paste0(prefix, "radial_r_density"))
  )
  safe_histogram_plot(
    df = df,
    colors = cluster_colors,
    title_text = paste0("k = ", k_value, " | ", unit_name, " | radial r histogram (bin=5)"),
    outfile = file.path(k_dir, paste0(prefix, "radial_r_histogram.png"))
  )

  # ==========================================
  # Export the underlying data used to draw the density plots
  # ==========================================
  fwrite(
    df[, .(x_center, y_center, cluster_raw, dx, dy, r, theta, angle_bin)], 
    file.path(k_dir, paste0(prefix, "radial_r_density_data.tsv")),
    sep = "\t"
  )
  # ==========================================

  # ==========================================
  # 2. Export density plots and histograms for the outer region only (r >= cutoff)
  # ==========================================
  df_plot <- df[r >= R_CUTOFF_FOR_EXTRA_DENSITY]
  if (nrow(df_plot) > 0) {
    safe_density_plot(
      df = df_plot,
      colors = cluster_colors,
      title_text = paste0(
        "k = ", k_value, " | ", unit_name,
        " | radial r density (r >= ", R_CUTOFF_FOR_EXTRA_DENSITY, ")"
      ),
      outfile = file.path(
        k_dir,
        paste0(prefix, "radial_r_density_rge", R_CUTOFF_FOR_EXTRA_DENSITY, ".png")
      ),
      cutoff = R_CUTOFF_FOR_EXTRA_DENSITY
    )

    safe_histogram_plot(
      df = df_plot,
      colors = cluster_colors,
      title_text = paste0(
        "k = ", k_value, " | ", unit_name,
        " | radial r histogram (r >= ", R_CUTOFF_FOR_EXTRA_DENSITY, ")"
      ),
      outfile = file.path(
        k_dir,
        paste0(prefix, "radial_r_histogram_rge", R_CUTOFF_FOR_EXTRA_DENSITY, ".png")
      ),
      cutoff = R_CUTOFF_FOR_EXTRA_DENSITY
    )
  }

  pair_summary[]
}

run_one_k_dir <- function(k_dir, k_value) {
  meta_file <- file.path(k_dir, "bin_metadata_with_cluster_raw.tsv")

  cat("====================================================\n")
  cat("Processing:", k_dir, "\n")
  cat("k =", k_value, "\n")

  if (!file.exists(meta_file)) {
    cat("[skip] missing file:", meta_file, "\n")
    return(NULL)
  }

  sub_dt <- fread(meta_file)
  ctx <- extract_context_from_path(k_dir)

  sample_root <- file.path(BASE_ROOT, ctx$sample_id)
  grid_root   <- file.path(sample_root, ctx$grid_dir)
  k10_root    <- file.path(grid_root, "kmeans_k10_raw_out")

  stem1_file <- file.path(k10_root, "bin_metadata_with_cluster_raw_stem1.tsv")
  stem2_file <- file.path(k10_root, "bin_metadata_with_cluster_raw_stem2.tsv")

  # dual-sample case
  if (file.exists(stem1_file) && file.exists(stem2_file)) {
    cat("Detected dual-sample folder. Splitting into stem1 / stem2 first...\n")
    cat("  stem1:", stem1_file, "\n")
    cat("  stem2:", stem2_file, "\n")

    stem1_dt <- fread(stem1_file)
    stem2_dt <- fread(stem2_file)

    split_dt <- split_subset_by_stem(sub_dt, stem1_dt, stem2_dt)

    n_miss <- sum(is.na(split_dt$stem))
    if (n_miss > 0) {
      cat("[warning]", n_miss, "bins could not be assigned to stem1/stem2\n")
    }

    res_list <- list()

    dt1 <- split_dt[stem == "stem1"]
    if (nrow(dt1) > 0) {
      res1 <- analyze_one_df(
        df = dt1,
        k_dir = k_dir,
        k_value = k_value,
        unit_name = "stem1",
        out_prefix = "stem1"
      )
      res_list[[length(res_list) + 1]] <- res1
    } else {
      cat("[skip] no bins assigned to stem1\n")
    }

    dt2 <- split_dt[stem == "stem2"]
    if (nrow(dt2) > 0) {
      res2 <- analyze_one_df(
        df = dt2,
        k_dir = k_dir,
        k_value = k_value,
        unit_name = "stem2",
        out_prefix = "stem2"
      )
      res_list[[length(res_list) + 1]] <- res2
    } else {
      cat("[skip] no bins assigned to stem2\n")
    }

    return(rbindlist(res_list, fill = TRUE))
  }

  # single-sample case
  cat("Single-sample folder detected. Analyzing whole subset directly...\n")
  analyze_one_df(
    df = sub_dt,
    k_dir = k_dir,
    k_value = k_value,
    unit_name = "whole",
    out_prefix = NULL
  )
}

############################################################
# main
############################################################
sample_dirs <- list.dirs(BASE_ROOT, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("^.*/output-", sample_dirs)]

all_k_dirs <- character()

for (sdir in sample_dirs) {
  subset_dirs <- list.dirs(sdir, recursive = TRUE, full.names = TRUE)
  subset_dirs <- subset_dirs[grepl(SUBSET_PATTERN, basename(subset_dirs))]

  for (subdir in subset_dirs) {
    for (k_value in K_VALUES) {
      k_dir <- file.path(subdir, paste0("kmeans_subset_k", sprintf("%02d", k_value), "_raw_out"))
      if (dir.exists(k_dir)) {
        all_k_dirs <- c(all_k_dirs, k_dir)
      }
    }
  }
}

all_k_dirs <- unique(all_k_dirs)

cat("====================================================\n")
cat("Found", length(all_k_dirs), "k directories to process.\n")
cat("====================================================\n")

all_summary <- rbindlist(lapply(all_k_dirs, function(k_dir) {
  if (grepl("kmeans_subset_k02_raw_out$", k_dir)) {
    run_one_k_dir(k_dir, 2)
  } else if (grepl("kmeans_subset_k03_raw_out$", k_dir)) {
    run_one_k_dir(k_dir, 3)
  } else if (grepl("kmeans_subset_k04_raw_out$", k_dir)) {
    run_one_k_dir(k_dir, 4)
  } else {
    NULL
  }
}), fill = TRUE)

if (!is.null(all_summary) && nrow(all_summary) > 0) {
  fwrite(all_summary, OUT_SUMMARY_ALL, sep = "\t")
  cat("Combined summary saved to:\n")
  cat(OUT_SUMMARY_ALL, "\n")
} else {
  cat("No valid results found.\n")
}

cat("Done.\n")