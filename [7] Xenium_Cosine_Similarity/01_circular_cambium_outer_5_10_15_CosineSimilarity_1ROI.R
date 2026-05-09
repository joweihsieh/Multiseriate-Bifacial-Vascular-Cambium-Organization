#20260409
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(RANN)
  library(Matrix)
  library(readxl)
})

# ============================================================
# USER SETTINGS
# ============================================================
BASE_DIR <- "/home/woodydrylab/FileShare/20260121_Xenium"

SUBSET_XLSX <- file.path(BASE_DIR, "subset_cambium_5um_k10_k3.xlsx")
K10_MAP_XLSX <- file.path(BASE_DIR, "k10_5domain.xlsx")

UM_DIR <- "grid05um_out"
SUB_K <- 3
K10_K <- 10

N_ANGLE_BINS <- 180
TRIM_Q <- 0.95
PUSH_STEPS_UM <- c(5, 10, 15)
NN_MAX_DIST <- 6

GLOBAL_TRIM_METHOD <- "quantile"   # "quantile" or "mad"
GLOBAL_TRIM_Q <- 0.95
GLOBAL_K_MAD <- 3

GRID_UM <- 5
SET_SEED <- 1
set.seed(SET_SEED)

# ============================================================
# Helpers
# ============================================================
parse_pairs_to_clusters <- function(x) {
  x <- gsub("\\s+", "", as.character(x))
  vals <- as.integer(unlist(strsplit(x, ",", fixed = TRUE)))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) stop("Cannot parse valid clusters from pairs: ", x)
  vals
}

get_sclerenchyma_cluster <- function(sample_id_input, k10_map_dt, k_use = 10, um_use = "grid05um_out") {
  sub <- k10_map_dt[
    sample_id == sample_id_input &
      k == k_use &
      um == um_use &
      tolower(domain) == "sclerenchyma"
  ]

  scl <- sort(unique(as.integer(sub$cluster)))
  scl <- scl[!is.na(scl)]

  if (length(scl) == 0) {
    stop("No sclerenchyma cluster found in k10_5domain.xlsx for sample_id_input: ", sample_id_input)
  }
  if (length(scl) > 1) {
    stop(
      "More than one sclerenchyma cluster found for sample_id_input: ", sample_id_input,
      " | clusters = ", paste(scl, collapse = ", ")
    )
  }

  scl
}

detect_shared_bin_id_col <- function(df_subset, df_full) {
  candidate_id_cols <- c(
    "bin_id","binID","barcode","bin_barcode","id","ID","BinID","BinId","BIN_ID"
  )

  bin_col <- candidate_id_cols[
    candidate_id_cols %in% names(df_subset) &
      candidate_id_cols %in% names(df_full)
  ][1]

  if (is.na(bin_col) || is.null(bin_col) || length(bin_col) == 0) {
    shared_cols <- intersect(names(df_subset), names(df_full))
    shared_char <- shared_cols[sapply(shared_cols, function(cc) {
      is.character(df_subset[[cc]]) || is.factor(df_subset[[cc]])
    })]

    if (length(shared_char) == 0) {
      stop("Could not auto-detect a shared bin id column. Please set BIN_ID_COL manually.")
    }

    uprop <- sapply(shared_char, function(cc) uniqueN(df_full[[cc]]) / nrow(df_full))
    bin_col <- shared_char[which.max(uprop)]
    cat(sprintf("Auto-picked BIN_ID_COL = '%s'\n", bin_col))
  } else {
    cat(sprintf("Detected shared BIN_ID_COL = '%s'\n", bin_col))
  }

  bin_col
}

safe_mkdir <- function(x) {
  dir.create(x, showWarnings = FALSE, recursive = TRUE)
}

# ============================================================
# Load Excel configs
# ============================================================
subset_cfg <- as.data.table(read_excel(SUBSET_XLSX))
k10_map_dt <- as.data.table(read_excel(K10_MAP_XLSX))

# normalize column names
if ("Sample_ID" %in% names(subset_cfg) && !("sample_id" %in% names(subset_cfg))) {
  setnames(subset_cfg, "Sample_ID", "sample_id")
}
if ("subset" %in% names(subset_cfg) && !("subset_dir" %in% names(subset_cfg))) {
  setnames(subset_cfg, "subset", "subset_dir")
}
if ("Subset" %in% names(subset_cfg) && !("subset_dir" %in% names(subset_cfg))) {
  setnames(subset_cfg, "Subset", "subset_dir")
}
if ("Pairs" %in% names(subset_cfg) && !("pairs" %in% names(subset_cfg))) {
  setnames(subset_cfg, "Pairs", "pairs")
}
if ("K" %in% names(subset_cfg) && !("k" %in% names(subset_cfg))) {
  setnames(subset_cfg, "K", "k")
}

if ("Sample_ID" %in% names(k10_map_dt) && !("sample_id" %in% names(k10_map_dt))) {
  setnames(k10_map_dt, "Sample_ID", "sample_id")
}
if ("Cluster" %in% names(k10_map_dt) && !("cluster" %in% names(k10_map_dt))) {
  setnames(k10_map_dt, "Cluster", "cluster")
}
if ("Domain" %in% names(k10_map_dt) && !("domain" %in% names(k10_map_dt))) {
  setnames(k10_map_dt, "Domain", "domain")
}
if ("UM" %in% names(k10_map_dt) && !("um" %in% names(k10_map_dt))) {
  setnames(k10_map_dt, "UM", "um")
}
if ("K" %in% names(k10_map_dt) && !("k" %in% names(k10_map_dt))) {
  setnames(k10_map_dt, "K", "k")
}

stopifnot(all(c("sample_id", "subset_dir", "k", "pairs") %in% names(subset_cfg)))
stopifnot(all(c("sample_id", "um", "k", "cluster", "domain") %in% names(k10_map_dt)))

subset_cfg[, sample_id := as.character(sample_id)]
subset_cfg[, subset_dir := as.character(subset_dir)]
subset_cfg[, k := as.integer(k)]
subset_cfg[, pairs := as.character(pairs)]

k10_map_dt[, sample_id := as.character(sample_id)]
k10_map_dt[, um := as.character(um)]
k10_map_dt[, k := as.integer(k)]
k10_map_dt[, cluster := as.integer(cluster)]
k10_map_dt[, domain := as.character(domain)]

subset_cfg <- subset_cfg[k == SUB_K]

cat("Samples to process:\n")
print(subset_cfg)

# ============================================================
# Main per-sample function
# ============================================================
run_one_sample <- function(sample_row) {
  sample_id_input <- sample_row$sample_id
  subset_dir <- sample_row$subset_dir
  keep_clusters <- parse_pairs_to_clusters(sample_row$pairs)


  sclerenchyma_cluster <- get_sclerenchyma_cluster(
    sample_id_input = sample_id_input,
    k10_map_dt = k10_map_dt,
    k_use = K10_K,
    um_use = UM_DIR
  )

  sample_label <- unique(
    k10_map_dt[
      sample_id == sample_id_input &
      um == UM_DIR &
      k == K10_K,
      label
    ]
  )

  sample_label <- sample_label[!is.na(sample_label)][1]
  cat("\n============================================================\n")
  cat("Processing sample_id: ", sample_id_input, "\n", sep = "")
  cat("subset_dir: ", subset_dir, "\n", sep = "")
  cat("pairs: ", sample_row$pairs, "\n", sep = "")
  cat("KEEP_CLUSTERS: c(", paste0(keep_clusters, "L", collapse = ", "), ")\n", sep = "")
  cat("k10 sclerenchyma cluster: ", sclerenchyma_cluster, "\n", sep = "")

  sample_base <- file.path(BASE_DIR, sample_id_input, UM_DIR)
  workdir <- file.path(sample_base, subset_dir, sprintf("kmeans_subset_k%02d_raw_out", SUB_K))
  subset_meta_file <- file.path(workdir, "bin_metadata_with_cluster_raw.tsv")
  full_meta_file <- file.path(sample_base, "bin_metadata.tsv")
  counts_rds <- file.path(sample_base, "counts_bins_by_genes_sparse.rds")
  k10_meta_file <- file.path(sample_base, "kmeans_k10_raw_out", "bin_metadata_with_cluster_raw.tsv")

  if (!file.exists(subset_meta_file)) stop("Missing subset metadata: ", subset_meta_file)
  if (!file.exists(full_meta_file)) stop("Missing full metadata: ", full_meta_file)
  if (!file.exists(counts_rds)) stop("Missing counts RDS: ", counts_rds)
  if (!file.exists(k10_meta_file)) stop("Missing k10 metadata: ", k10_meta_file)

  out_dir <- file.path(workdir, "radial_push_with_sclerenchyma_auto_20260409")
  safe_mkdir(out_dir)

  # ------------------------------------------------------------
  # 0) Load subset + full metadata
  # ------------------------------------------------------------
  df_subset0 <- fread(subset_meta_file)
  stopifnot(all(c("x_center", "y_center", "cluster_raw") %in% names(df_subset0)))
  df_subset0[, cluster_raw := as.integer(cluster_raw)]
  df_subset <- df_subset0[cluster_raw %in% keep_clusters]

  cat("Subset counts per cluster_raw:\n")
  print(df_subset[, .N, by = cluster_raw][order(cluster_raw)])

  if (nrow(df_subset) == 0) {
    stop("No bins left after filtering KEEP_CLUSTERS in sample: ", sample_id_input)
  }

  df_full <- fread(full_meta_file)
  stopifnot(all(c("x_center", "y_center") %in% names(df_full)))

  # ------------------------------------------------------------
  # 0.1) Detect BIN_ID_COL shared by subset/full
  # ------------------------------------------------------------
  bin_id_col <- detect_shared_bin_id_col(df_subset, df_full)
  df_subset[, (bin_id_col) := as.character(get(bin_id_col))]
  df_full[, (bin_id_col) := as.character(get(bin_id_col))]

  # ------------------------------------------------------------
  # 1) Compute polar coords + GLOBAL r trimming
  # ------------------------------------------------------------
  center_x <- mean(df_subset$x_center)
  center_y <- mean(df_subset$y_center)

  df_subset[, dx := x_center - center_x]
  df_subset[, dy := y_center - center_y]
  df_subset[, theta := atan2(dy, dx)]
  df_subset[, r := sqrt(dx^2 + dy^2)]

  if (GLOBAL_TRIM_METHOD == "quantile") {
    r_cut_global <- as.numeric(quantile(df_subset$r, probs = GLOBAL_TRIM_Q, na.rm = TRUE))
    df_subset_clean <- df_subset[r <= r_cut_global]
    cat(sprintf("GLOBAL r trim (quantile): keep r <= %.3f (q=%.3f)\n", r_cut_global, GLOBAL_TRIM_Q))
  } else if (GLOBAL_TRIM_METHOD == "mad") {
    r_med <- median(df_subset$r, na.rm = TRUE)
    r_mad <- mad(df_subset$r, constant = 1.4826, na.rm = TRUE)
    if (is.na(r_mad) || r_mad == 0) r_mad <- 1
    df_subset_clean <- df_subset[abs(r - r_med) <= GLOBAL_K_MAD * r_mad]
    cat(sprintf("GLOBAL r trim (MAD): keep |r-med| <= %.3f (K=%g)\n", GLOBAL_K_MAD * r_mad, GLOBAL_K_MAD))
  } else {
    stop("GLOBAL_TRIM_METHOD must be 'quantile' or 'mad'")
  }

  cat(sprintf("Subset bins: before=%d, after_global_trim=%d\n", nrow(df_subset), nrow(df_subset_clean)))

  if (nrow(df_subset_clean) == 0) {
    stop("No bins left after global trimming in sample: ", sample_id_input)
  }

  breaks <- seq(-pi, pi, length.out = N_ANGLE_BINS + 1)
  df_subset_clean[, angle_bin := cut(theta, breaks = breaks, include.lowest = TRUE)]

  # ------------------------------------------------------------
  # 2) Anchor selection
  # ------------------------------------------------------------
  anchor_dt <- df_subset_clean[, {
    rr <- r
    cutoff <- as.numeric(quantile(rr, probs = TRIM_Q, na.rm = TRUE))
    dd <- .SD[r <= cutoff]
    if (nrow(dd) == 0) dd <- .SD
    dd[which.max(r), .(
      anchor_bin_id = get(bin_id_col),
      anchor_x = x_center,
      anchor_y = y_center,
      anchor_r = r,
      anchor_theta = theta,
      n_in_anglebin = .N,
      trim_cutoff_r = cutoff
    )]
  }, by = angle_bin]

  fwrite(anchor_dt, file.path(out_dir, "anglebin_outer_anchor_trimmed.tsv"), sep = "\t")
  cat(sprintf("Saved anchors (%d)\n", nrow(anchor_dt)))

  # ------------------------------------------------------------
  # 3) Push outward + match bins on FULL metadata
  # ------------------------------------------------------------
  targets <- anchor_dt[, {
    ux <- cos(anchor_theta)
    uy <- sin(anchor_theta)
    rbindlist(lapply(PUSH_STEPS_UM, function(s) {
      data.table(
        step_um = s,
        target_x = anchor_x + s * ux,
        target_y = anchor_y + s * uy
      )
    }))
  }, by = .(angle_bin, anchor_bin_id, anchor_x, anchor_y, anchor_r, anchor_theta)]

  cand_xy <- as.matrix(df_full[, .(x_center, y_center)])
  query_xy <- as.matrix(targets[, .(target_x, target_y)])

  nn <- RANN::nn2(data = cand_xy, query = query_xy, k = 1)
  targets[, nn_row := nn$nn.idx[, 1]]
  targets[, nn_dist_to_target := nn$nn.dists[, 1]]

  targets[, matched_bin_id := df_full[[bin_id_col]][nn_row]]
  targets[, matched_x := df_full$x_center[nn_row]]
  targets[, matched_y := df_full$y_center[nn_row]]

  targets[, dist_anchor_to_matched := sqrt((matched_x - anchor_x)^2 + (matched_y - anchor_y)^2)]

  targets[, is_valid := nn_dist_to_target <= NN_MAX_DIST]
  targets[is_valid == FALSE, `:=`(
    matched_bin_id = NA_character_,
    matched_x = NA_real_,
    matched_y = NA_real_,
    dist_anchor_to_matched = NA_real_
  )]

  fwrite(targets, file.path(out_dir, "radial_push3steps_matched_bins_fullcand.tsv"), sep = "\t")
  cat(sprintf("Saved: %s\n", file.path(out_dir, "radial_push3steps_matched_bins_fullcand.tsv")))
  cat(sprintf("Valid target rate: %.2f%% (NN_MAX_DIST=%g)\n", 100 * mean(targets$is_valid), NN_MAX_DIST))

  # ------------------------------------------------------------
  # 4) Visualization background + anchors + pushes + sclerenchyma
  # ------------------------------------------------------------
  anc_nn <- RANN::nn2(
    data = as.matrix(df_full[, .(x_center, y_center)]),
    query = as.matrix(anchor_dt[, .(anchor_x, anchor_y)]),
    k = 1
  )

  anchor_plot <- copy(anchor_dt)
  anchor_plot[, nn_row_full := anc_nn$nn.idx[, 1]]
  anchor_plot[, anchor_xg := df_full$x_center[nn_row_full]]
  anchor_plot[, anchor_yg := df_full$y_center[nn_row_full]]
  anchor_plot[, nn_dist_anchor_to_full := anc_nn$nn.dists[, 1]]

  cat("Anchor->fullgrid nn dist summary:\n")
  print(summary(anchor_plot$nn_dist_anchor_to_full))

  push_plot <- unique(targets[is_valid == TRUE, .(step_um, xg = matched_x, yg = matched_y)])

  df_k10 <- fread(k10_meta_file)
  df_k10[, cluster_raw := as.integer(cluster_raw)]
  df_k10_scl <- df_k10[cluster_raw == sclerenchyma_cluster]

  cat(sprintf("k10 sclerenchyma bins: %d\n", nrow(df_k10_scl)))

  bg_plot <- df_subset_clean[, .(
    x = x_center,
    y = y_center,
    category = "subset"
  )]

  scl_plot <- df_k10_scl[, .(
    x = x_center,
    y = y_center,
    category = "sclerenchyma"
  )]

  anchor_tile_plot <- anchor_plot[, .(
    x = anchor_xg,
    y = anchor_yg,
    category = "anchor"
  )]

  push_tile_plot <- push_plot[, .(
    x = xg,
    y = yg,
    category = paste0("step_", step_um)
  )]

  all_tiles <- rbindlist(
    list(bg_plot, scl_plot, anchor_tile_plot, push_tile_plot),
    fill = TRUE
  )

  all_tiles[, category := factor(
    category,
    levels = c("subset", "sclerenchyma", "anchor", "step_5", "step_10", "step_15")
  )]

  tile_colors <- c(
    subset = "#F8766D",
    sclerenchyma = "#CE9332",
    anchor = "black",
    step_5 = "#F8766D",
    step_10 = "#FBA9A5",
    step_15 = "#FEECEC"
  )

  p_grid_alltiles <- ggplot(all_tiles) +
    geom_tile(
      aes(x = x, y = y, fill = category),
      width = GRID_UM,
      height = GRID_UM
    ) +
    scale_fill_manual(values = tile_colors, drop = FALSE) +
    theme_classic() +
    theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank()
    )+
    labs(
      fill = "",
      title = paste0("\n", sample_label)
    ) +
    coord_equal() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0))

  ggsave(
    file.path(out_dir, "viz_grid_tiles_all_clean_plus_sclerenchyma.png"),
    p_grid_alltiles,
    width = 7.8,
    height = 6.8,
    dpi = 600
  )

  tile_fill_colors_blackwhite <- c(
    subset = "#F8766D",
    sclerenchyma = "#CE9332",
    anchor = "black",
    step_5 = "white",
    step_10 = "white",
    step_15 = "white"
  )

  # 
  tile_border_colors <- c(
    subset = "transparent",
    sclerenchyma = "transparent",
    anchor = "transparent",
    step_5 = "black",
    step_10 = "black",
    step_15 = "black"
  )
  p_grid_alltiles_black <- ggplot(all_tiles) +
    geom_tile(
      aes(x = x, y = y, fill = category, color = category),
      width = GRID_UM,
      height = GRID_UM,
      linewidth = 0.05 
    ) +
    scale_fill_manual(values = tile_fill_colors_blackwhite, drop = FALSE) +
    scale_color_manual(values = tile_border_colors, drop = FALSE, guide = "none") + 
    theme_classic() +
    theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank()
    )+
    labs(
      fill = "",
      title = paste0("\n", sample_label)
    ) +
    coord_equal() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0))
  
  ggsave(
    file.path(out_dir, "viz_grid_tiles_all_clean_plus_sclerenchyma_black.png"),
    p_grid_alltiles_black,
    width = 7.8,
    height = 6.8,
    dpi = 1200
  )

  cat(sprintf("Saved: %s\n", file.path(out_dir, "viz_grid_tiles_all_clean_plus_sclerenchyma.png")))

  # ------------------------------------------------------------
  # 5) Extract transcriptomes from sparse counts
  # ------------------------------------------------------------
  counts <- readRDS(counts_rds)

  if (inherits(counts, "dgCMatrix") || inherits(counts, "dgTMatrix") || inherits(counts, "Matrix")) {
    # OK
  } else {
    stop("counts_rds is not a sparse Matrix object: ", counts_rds)
  }

  if (is.null(rownames(counts))) {
    stop("counts matrix has no rownames; expected rownames = bin IDs")
  }

  bin_ids_in_counts <- rownames(counts)

  anchor_dt2 <- anchor_plot[anchor_bin_id %in% bin_ids_in_counts]
  targets2 <- targets[
    is_valid == TRUE &
      !is.na(matched_bin_id) &
      matched_bin_id %in% bin_ids_in_counts
  ]

  cat(sprintf("Anchor bins in counts: %d / %d\n", nrow(anchor_dt2), nrow(anchor_plot)))
  cat(sprintf(
    "Valid pushed bins in counts: %d / %d\n",
    nrow(targets2),
    nrow(targets[is_valid == TRUE & !is.na(matched_bin_id)])
  ))

  anchor_mat <- counts[anchor_dt2$anchor_bin_id, , drop = FALSE]
  pushed_mat <- counts[targets2$matched_bin_id, , drop = FALSE]

  saveRDS(anchor_mat, file.path(out_dir, "anchor_transcriptomes_binsXgenes.rds"))
  saveRDS(pushed_mat, file.path(out_dir, "pushed3steps_transcriptomes_binsXgenes.rds"))
  cat(sprintf("Saved: %s\n", file.path(out_dir, "anchor_transcriptomes_binsXgenes.rds")))
  cat(sprintf("Saved: %s\n", file.path(out_dir, "pushed3steps_transcriptomes_binsXgenes.rds")))

  # ------------------------------------------------------------
  # 6) Transcriptome distance on RAW counts
  # ------------------------------------------------------------
  anchor_lookup <- anchor_dt2[, .(angle_bin, anchor_bin_id_used = anchor_bin_id)]

  targets_dist <- merge(
    targets2,
    anchor_lookup,
    by = "angle_bin",
    all.x = TRUE,
    sort = FALSE
  )

  targets_dist <- targets_dist[
    !is.na(anchor_bin_id_used) &
      anchor_bin_id_used %in% bin_ids_in_counts &
      matched_bin_id %in% bin_ids_in_counts
  ]

  cat(sprintf("Pairs for transcriptome distance (raw): %d\n", nrow(targets_dist)))

  X <- counts[targets_dist$anchor_bin_id_used, , drop = FALSE]
  Y <- counts[targets_dist$matched_bin_id, , drop = FALSE]

  dot_xy <- Matrix::rowSums(X * Y)
  norm_x <- sqrt(Matrix::rowSums(X * X))
  norm_y <- sqrt(Matrix::rowSums(Y * Y))
  den <- norm_x * norm_y
  den[den == 0] <- NA_real_

  cosine_dist <- 1 - as.numeric(dot_xy / den)
  euclid_dist <- sqrt(Matrix::rowSums((X - Y) * (X - Y)))

  targets_dist[, transcriptome_cosine_dist_raw := cosine_dist]
  targets_dist[, transcriptome_euclid_dist_raw := as.numeric(euclid_dist)]

  dist_tbl <- targets_dist[, .(
    angle_bin,
    step_um,
    anchor_bin_id = anchor_bin_id_used,
    matched_bin_id,
    nn_dist_to_target,
    dist_anchor_to_matched,
    transcriptome_cosine_dist_raw,
    transcriptome_euclid_dist_raw
  )]

  fwrite(dist_tbl, file.path(out_dir, "transcriptome_distance_raw_anchor_vs_push.tsv"), sep = "\t")
  cat(sprintf("Saved: %s\n", file.path(out_dir, "transcriptome_distance_raw_anchor_vs_push.tsv")))

  p_dist <- ggplot(dist_tbl, aes(x = factor(step_um), y = transcriptome_cosine_dist_raw)) +
    geom_boxplot(outlier.size = 0.4) +
    theme_classic() +
    labs(
      x = "step (um)",
      y = "Cosine distance",
      title = paste0("\n", sample_label)
    )

  ggsave(
    file.path(out_dir, "transcriptome_cosine_distance_raw_by_step.png"),
    p_dist,
    width = 6,
    height = 6,
    dpi = 600
  )

  p_dist_color <- ggplot(dist_tbl, aes(
    x = factor(step_um, levels = c("15", "10", "5")),
    y = transcriptome_cosine_dist_raw, 
    fill = factor(step_um, levels = c("15", "10", "5"))
    )) +
    geom_boxplot(outlier.size = 0.4) +
    scale_fill_manual(values = c("5" = "#F8766D", "10" = "#FBA9A5", "15" = "#FEECEC")) +    
    theme_classic() +
    theme(legend.position = "none") +
    labs(
      x = "step (um)",
      y = "Cosine distance",
      title = paste0("\n", sample_label)
    )

  ggsave(
    file.path(out_dir, "transcriptome_cosine_distance_raw_by_step_color.png"),
    p_dist_color,
    width = 6,
    height = 6,
    dpi = 600
  )

  cat(sprintf("Saved: %s\n", file.path(out_dir, "transcriptome_cosine_distance_raw_by_step.png")))

  invisible(data.table(
    sample_id = sample_id_input,
    subset_dir = subset_dir,
    keep_clusters = paste(keep_clusters, collapse = ","),
    sclerenchyma_cluster = sclerenchyma_cluster,
    out_dir = out_dir,
    n_subset = nrow(df_subset),
    n_subset_clean = nrow(df_subset_clean),
    n_anchor = nrow(anchor_dt),
    valid_push_rate = mean(targets$is_valid),
    n_dist_pairs = nrow(dist_tbl)
  ))
}

# ============================================================
# Run all samples
# ============================================================
results <- vector("list", nrow(subset_cfg))

for (i in seq_len(nrow(subset_cfg))) {
  sample_row <- subset_cfg[i]
  results[[i]] <- run_one_sample(sample_row)
}

summary_dt <- rbindlist(results, use.names = TRUE, fill = TRUE)
summary_file <- file.path(BASE_DIR, paste0("radial_push_with_sclerenchyma_auto_summary_k", SUB_K, ".tsv"))
fwrite(summary_dt, summary_file, sep = "\t")

cat("\n============================================================\n")
cat("DONE.\n")
cat("Summary saved to:\n")
cat(summary_file, "\n")

#### draw boxplot using transcriptome_distance_raw_anchor_vs_push.tsv
