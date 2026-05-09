#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(readxl)
})

# ============================================================
# FILES / SETTINGS
# ============================================================
BASE_DIR <- "/home/woodydrylab/FileShare/20260121_Xenium"

UMAP_FILE <- "/home/chingweilu/Multi_Bilayer_Unifacial/Adjusting_UMAP_neighbor/UMAP_csv/WeiLi_projection_UMAP.csv"
TENX_FILE <- "/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/all_UMI_tables/geneUMI_TenX_PtrWT2forWOX2_v701.csv"

SUBSET_XLSX_SINGLE <- file.path(BASE_DIR, "subset_cambium_5um_k10_k3.xlsx")
SUBSET_XLSX_TWO    <- file.path(BASE_DIR, "subset_cambium_5um_k10_k3_two.xlsx")
K10_MAP_XLSX       <- file.path(BASE_DIR, "k10_5domain.xlsx")

UM_DIR <- "grid05um_out"
COUNTS_FILE_NAME <- "counts_subset.rds"
OUT_BASENAME <- "WeiLi_UMAP_XeniumCorr_20260409"

# ============================================================
# Gene ID cleaning
# ============================================================
clean_potri <- function(x) {
  x <- sub("\\.v[0-9]+\\.[0-9]+$", "", x)
  x <- gsub("^Potri\\.", "Potri_", x)
  x
}

# ============================================================
# Helpers
# ============================================================
normalize_subset_dir <- function(x) {
  x <- as.character(x)

  is_target <- grepl("^subset_from_k10_", x)

  x[is_target] <- vapply(x[is_target], function(s) {
    prefix <- "subset_from_k10_"
    rest <- sub("^subset_from_k10_", "", s)

    parts <- strsplit(rest, "_", fixed = TRUE)[[1]]

    parts2 <- vapply(parts, function(p) {
      if (grepl("^c[0-9]+$", p)) {
        num <- as.integer(sub("^c", "", p))
        sprintf("c%02d", num)
      } else {
        p
      }
    }, character(1))

    paste0(prefix, paste(parts2, collapse = "_"))
  }, character(1))

  x
}

is_two_cluster_subset <- function(x) {
  grepl("^subset_from_k10_c[0-9]{2}_c[0-9]{2}$", x)
}

detect_shared_bin_id_col <- function(df1, df2) {
  candidate_id_cols <- c(
    "bin_id","binID","barcode","bin_barcode","id","ID","BinID","BinId","BIN_ID"
  )

  bin_col <- candidate_id_cols[
    candidate_id_cols %in% names(df1) &
      candidate_id_cols %in% names(df2)
  ][1]

  if (is.na(bin_col) || is.null(bin_col) || length(bin_col) == 0) {
    shared_cols <- intersect(names(df1), names(df2))
    shared_char <- shared_cols[sapply(shared_cols, function(cc) {
      is.character(df1[[cc]]) || is.factor(df1[[cc]])
    })]

    if (length(shared_char) == 0) {
      stop("Could not auto-detect a shared bin id column.")
    }

    uprop <- sapply(shared_char, function(cc) uniqueN(df2[[cc]]) / nrow(df2))
    bin_col <- shared_char[which.max(uprop)]
    cat(sprintf("Auto-picked shared BIN_ID_COL = '%s'\n", bin_col))
  } else {
    cat(sprintf("Detected shared BIN_ID_COL = '%s'\n", bin_col))
  }

  bin_col
}

detect_row_id_col <- function(df) {
  candidate_id_cols <- c(
    "bin_id","binID","barcode","bin_barcode","id","ID","BinID","BinId","BIN_ID"
  )
  hit <- candidate_id_cols[candidate_id_cols %in% names(df)]
  if (length(hit) >= 1) return(hit[1])

  char_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
  if (length(char_cols) == 0) {
    stop("Could not detect row id column from metadata.")
  }
  char_cols[1]
}

split_combined_labels <- function(x) {
  labs <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE)))
  labs[nzchar(labs)]
}

read_subset_cfg <- function(fp, sample_type_label) {
  dt <- as.data.table(read_excel(fp))

  if ("Sample_ID" %in% names(dt) && !("sample_id" %in% names(dt))) setnames(dt, "Sample_ID", "sample_id")
  if ("subset" %in% names(dt) && !("subset_dir" %in% names(dt))) setnames(dt, "subset", "subset_dir")
  if ("Subset" %in% names(dt) && !("subset_dir" %in% names(dt))) setnames(dt, "Subset", "subset_dir")
  if ("K" %in% names(dt) && !("k" %in% names(dt))) setnames(dt, "K", "k")
  if ("UM" %in% names(dt) && !("um" %in% names(dt))) setnames(dt, "UM", "um")

  stopifnot(all(c("sample_id", "subset_dir") %in% names(dt)))

  dt[, sample_id := as.character(sample_id)]
  dt[, subset_dir := as.character(subset_dir)]
  dt[, subset_dir := normalize_subset_dir(subset_dir)]

  if ("k" %in% names(dt)) dt[, k := as.integer(k)]
  if ("um" %in% names(dt)) dt[, um := as.character(um)]

  dt[, sample_type := sample_type_label]
  dt
}

# ============================================================
# Robust UMAP loader
# ============================================================
load_umap <- function(fp) {
  umap <- fread(fp)

  if (!("barcode" %in% names(umap))) {
    bc_cand <- names(umap)[grepl("barcode|cell|spot", names(umap), ignore.case = TRUE)]
    if (length(bc_cand) >= 1) {
      setnames(umap, bc_cand[1], "barcode")
    } else {
      setnames(umap, names(umap)[1], "barcode")
    }
  }

  nms <- names(umap)
  u1_cand <- nms[grepl("^umap[_\\-\\.]?1$|umap.*1$", nms, ignore.case = TRUE)]
  u2_cand <- nms[grepl("^umap[_\\-\\.]?2$|umap.*2$", nms, ignore.case = TRUE)]

  pick_best <- function(cands, target_num) {
    if (length(cands) == 0) return(character(0))
    exact <- cands[grepl(paste0("^umap[_\\-\\.]?", target_num, "$"), cands, ignore.case = TRUE)]
    if (length(exact) >= 1) return(exact[1])
    cands[which.min(nchar(cands))]
  }

  u1 <- pick_best(u1_cand, 1)
  u2 <- pick_best(u2_cand, 2)

  if (length(u1) == 0 || length(u2) == 0) {
    stop("Cannot detect UMAP columns. Please check names(UMAP csv).")
  }

  setnames(umap, c(u1, u2), c("UMAP_1", "UMAP_2"))
  umap[, .(barcode, UMAP_1, UMAP_2)]
}

# ============================================================
# Load TenX: row=barcode, col=genes -> gene x barcode
# ============================================================
load_tenx_gene_by_cell <- function(fp) {
  tenx <- fread(fp)

  if (!("Barcode" %in% names(tenx))) setnames(tenx, names(tenx)[1], "Barcode")
  setnames(tenx, "Barcode", "barcode")

  gene_cols_raw <- setdiff(names(tenx), "barcode")
  gene_cols_clean <- clean_potri(gene_cols_raw)
  setnames(tenx, gene_cols_raw, gene_cols_clean)

  barcodes <- tenx$barcode
  mat_bg <- as.matrix(tenx[, ..gene_cols_clean])
  rownames(mat_bg) <- barcodes

  mat_gb <- t(mat_bg)
  log1p(mat_gb)
}

# ============================================================
# Load Xenium subset: row = bin, col = gene -> bin x gene
# ============================================================
load_xenium_bin_by_gene <- function(fp) {
  x <- readRDS(fp)
  x <- as.matrix(x)

  if (is.null(colnames(x))) stop("Xenium matrix has no gene colnames: ", fp)
  if (is.null(rownames(x))) stop("Xenium matrix has no rownames (bin IDs): ", fp)

  colnames(x) <- clean_potri(colnames(x))
  log1p(x)
}

# ============================================================
# Plot function
# ============================================================
plot_umap_correlation_base <- function(
    dt,
    value_col,
    main_text = "",
    output_png,
    sorted_order = TRUE,
    point_cex = 0.75,
    point_pch = 20,
    color_tick = c("#EEF2F9", "#C44233")
) {
  x <- dt$UMAP_1
  y <- dt$UMAP_2
  z <- dt[[value_col]]

  ok <- !is.na(x) & !is.na(y) & !is.na(z)
  x <- x[ok]
  y <- y[ok]
  z <- z[ok]

  if (length(z) == 0) {
    warning("No non-NA points for ", value_col)
    return(invisible(NULL))
  }

  max_bound <- as.numeric(quantile(z, 0.99, na.rm = TRUE))
  upper_bound <- max_bound * 0.99
  lower_bound <- max_bound * 0.75

  if (!is.finite(lower_bound) || !is.finite(upper_bound) || upper_bound <= lower_bound) {
    lower_bound <- min(z, na.rm = TRUE)
    upper_bound <- max(z, na.rm = TRUE)
    if (upper_bound <= lower_bound) upper_bound <- lower_bound + 1e-8
  }

  color_index <- round(pmin(pmax((z - lower_bound) / (upper_bound - lower_bound), 0), 1) * 500) + 1
  color_pool <- colorRampPalette(color_tick)(501)

  if (sorted_order) {
    plot_order <- order(z)
  } else {
    set.seed(1)
    plot_order <- sample(seq_along(z))
  }

  png(output_png, pointsize = 10, res = 300, width = 15, height = 15, units = "cm")
  par(mai = c(0.7, 0.7, 0.9, 0.5))

  plot(
    x[plot_order], y[plot_order],
    col = color_pool[color_index[plot_order]],
    pch = point_pch,
    cex = point_cex,
    xlab = "UMAP.1",
    ylab = "UMAP.2",
    main = main_text,
    las = 1,
    asp = 1
  )

  par(xpd = TRUE)
  current_xlim <- par()$usr[1:2]
  current_ylim <- par()$usr[3:4]

  legend_x_range <- c(sum(current_xlim * c(0.1, 0.9)), current_xlim[2])
  legend_x_vector <- seq(
    from = legend_x_range[1],
    to = legend_x_range[2],
    length.out = length(color_pool) + 1
  )
  legend_y_range <- c(
    sum(current_ylim * c(-0.1, 1.1)),
    sum(current_ylim * c(-0.15, 1.15))
  )

  for (i in seq_along(color_pool)) {
    lines(
      legend_x_vector[c(i, i + 1)],
      rep(legend_y_range[1], 2),
      col = color_pool[i],
      lwd = 20,
      lend = "square"
    )
  }

  text(legend_x_range[1], legend_y_range[2], sprintf("%.2f", lower_bound), cex = 1)
  text(legend_x_range[2], legend_y_range[2], sprintf("%.2f", upper_bound), cex = 1)
  par(xpd = FALSE)

  dev.off()
}

# ============================================================
# Split Xenium subset into stem1/stem2 if needed
# ============================================================
split_subset_matrix_by_stem <- function(sample_base, subset_dir, xenium_bg, force_split = FALSE) {
  k10_dir <- file.path(sample_base, "kmeans_k10_raw_out")
  stem1_file <- file.path(k10_dir, "bin_metadata_with_cluster_raw_stem1.tsv")
  stem2_file <- file.path(k10_dir, "bin_metadata_with_cluster_raw_stem2.tsv")

  subset_meta_stem1 <- file.path(sample_base, subset_dir, "bin_metadata_with_cluster_raw_stem1.tsv")
  subset_meta_stem2 <- file.path(sample_base, subset_dir, "bin_metadata_with_cluster_raw_stem2.tsv")

  subset_bin_ids <- rownames(xenium_bg)
  if (is.null(subset_bin_ids)) {
    stop("xenium_bg has no rownames, cannot split by stem.")
  }
  subset_bin_ids <- as.character(subset_bin_ids)

  out <- list()
  out_ids <- list()

  # case 1: subset has stem1/stem2 metadata
  if (file.exists(subset_meta_stem1) || file.exists(subset_meta_stem2)) {
    if (file.exists(subset_meta_stem1)) {
      meta1 <- fread(subset_meta_stem1)
      idcol1 <- detect_row_id_col(meta1)
      ids1 <- intersect(subset_bin_ids, as.character(meta1[[idcol1]]))
      cat(sprintf("subset stem1 direct split: %d bins\n", length(ids1)))
      if (length(ids1) > 0) {
        out[["stem1"]] <- xenium_bg[ids1, , drop = FALSE]
        out_ids[["stem1"]] <- ids1
      }
    }

    if (file.exists(subset_meta_stem2)) {
      meta2 <- fread(subset_meta_stem2)
      idcol2 <- detect_row_id_col(meta2)
      ids2 <- intersect(subset_bin_ids, as.character(meta2[[idcol2]]))
      cat(sprintf("subset stem2 direct split: %d bins\n", length(ids2)))
      if (length(ids2) > 0) {
        out[["stem2"]] <- xenium_bg[ids2, , drop = FALSE]
        out_ids[["stem2"]] <- ids2
      }
    }
  } else {
    # case 2: 用 parent k10 stem1/stem2 metadata 直接和 subset rownames 取交集
    if (file.exists(stem1_file)) {
      meta1 <- fread(stem1_file)
      idcol1 <- detect_row_id_col(meta1)
      ids1 <- intersect(subset_bin_ids, as.character(meta1[[idcol1]]))
      cat(sprintf("parent stem1 split: %d bins\n", length(ids1)))
      if (length(ids1) > 0) {
        out[["stem1"]] <- xenium_bg[ids1, , drop = FALSE]
        out_ids[["stem1"]] <- ids1
      }
    }

    if (file.exists(stem2_file)) {
      meta2 <- fread(stem2_file)
      idcol2 <- detect_row_id_col(meta2)
      ids2 <- intersect(subset_bin_ids, as.character(meta2[[idcol2]]))
      cat(sprintf("parent stem2 split: %d bins\n", length(ids2)))
      if (length(ids2) > 0) {
        out[["stem2"]] <- xenium_bg[ids2, , drop = FALSE]
        out_ids[["stem2"]] <- ids2
      }
    }
  }

  if (length(out) == 0) {
    if (force_split) {
      stop("Force split requested, but no usable stem1/stem2 bins found under: ", k10_dir)
    } else {
      return(list(mats = list(combined = xenium_bg), ids = list(combined = subset_bin_ids)))
    }
  }

  if ("stem1" %in% names(out_ids) && "stem2" %in% names(out_ids)) {
    ids1s <- sort(unique(out_ids[["stem1"]]))
    ids2s <- sort(unique(out_ids[["stem2"]]))
    overlap_n <- length(intersect(ids1s, ids2s))
    cat(sprintf("stem1/stem2 overlap bins: %d\n", overlap_n))

    if (identical(ids1s, ids2s)) {
      stop("stem1 and stem2 got identical bin sets after split. This means the split failed logically.")
    }
  }

  list(mats = out, ids = out_ids)
}

# ============================================================
# Load k10 map for labels
# ============================================================
k10_map_dt <- as.data.table(read_excel(K10_MAP_XLSX))

if ("Sample_ID" %in% names(k10_map_dt) && !("sample_id" %in% names(k10_map_dt))) {
  setnames(k10_map_dt, "Sample_ID", "sample_id")
}
if ("Label" %in% names(k10_map_dt) && !("label" %in% names(k10_map_dt))) {
  setnames(k10_map_dt, "Label", "label")
}

k10_map_dt[, sample_id := as.character(sample_id)]
if ("label" %in% names(k10_map_dt)) {
  k10_map_dt[, label := as.character(label)]
} else {
  stop("k10_5domain.xlsx does not contain label/Label column.")
}

# ============================================================
# Load config Excel
# ============================================================
cfg_single <- read_subset_cfg(SUBSET_XLSX_SINGLE, "single")
cfg_two    <- read_subset_cfg(SUBSET_XLSX_TWO, "two")

cfg <- rbindlist(list(cfg_single, cfg_two), fill = TRUE, use.names = TRUE)

if ("um" %in% names(cfg)) {
  cfg <- cfg[um == UM_DIR]
}

cfg <- cfg[is_two_cluster_subset(subset_dir)]
cfg <- unique(cfg[, .(sample_id, subset_dir, sample_type)])

cat("Jobs to process:\n")
print(cfg)

# ============================================================
# Shared data
# ============================================================
cat("Loading UMAP...\n")
umap <- load_umap(UMAP_FILE)

cat("Loading TenX...\n")
mat_tenx_gb <- load_tenx_gene_by_cell(TENX_FILE)

# ============================================================
# Per-subset runner
# ============================================================
run_one_subset <- function(sample_id_input, subset_dir_input, sample_type_input) {
  subset_dir_input <- normalize_subset_dir(subset_dir_input)

  if (!is_two_cluster_subset(subset_dir_input)) {
    warning("Skip non-two-cluster subset: ", subset_dir_input)
    return(data.table(
      sample_id = sample_id_input,
      subset_dir = subset_dir_input,
      sample_type = sample_type_input,
      stem_id = NA_character_,
      status = "skip_non_two_cluster_subset"
    ))
  }

  cat("\n============================================================\n")
  cat("Processing sample_id :", sample_id_input, "\n")
  cat("Processing subset_dir:", subset_dir_input, "\n")
  cat("sample_type         :", sample_type_input, "\n")

  sample_base <- file.path(BASE_DIR, sample_id_input, UM_DIR)
  subset_base <- file.path(sample_base, subset_dir_input)
  xenium_file <- file.path(subset_base, COUNTS_FILE_NAME)

  if (!file.exists(xenium_file)) {
    warning("Missing Xenium file: ", xenium_file)
    return(data.table(
      sample_id = sample_id_input,
      subset_dir = subset_dir_input,
      sample_type = sample_type_input,
      stem_id = NA_character_,
      status = "missing_counts_subset"
    ))
  }

  # label
  sample_label_all <- unique(k10_map_dt[sample_id == sample_id_input, label])
  sample_label_all <- sample_label_all[!is.na(sample_label_all)][1]
  if (is.na(sample_label_all) || length(sample_label_all) == 0) {
    sample_label_all <- sample_id_input
  }
  split_labels <- split_combined_labels(sample_label_all)

  cat("Loading Xenium...\n")
  xenium_bg <- load_xenium_bin_by_gene(xenium_file)

  if (sample_type_input == "two") {
    cat("Two-sample entry: force split into stem1/stem2...\n")
    split_res <- split_subset_matrix_by_stem(
      sample_base = sample_base,
      subset_dir = subset_dir_input,
      xenium_bg = xenium_bg,
      force_split = TRUE
    )
  } else {
    cat("Single-sample entry: run as combined only...\n")
    split_res <- list(
      mats = list(combined = xenium_bg),
      ids  = list(combined = rownames(xenium_bg))
    )
  }

  split_mats <- split_res$mats
  split_ids  <- split_res$ids

  res_one <- vector("list", length(split_mats))
  idx <- 1L

  for (stem_id in names(split_mats)) {
    cat("Running stem:", stem_id, "\n")

    xenium_part_bg <- split_mats[[stem_id]]
    stem_bin_ids   <- split_ids[[stem_id]]

    cat(sprintf("Stem %s uses %d bins\n", stem_id, length(stem_bin_ids)))

    if (nrow(xenium_part_bg) == 0) {
      res_one[[idx]] <- data.table(
        sample_id = sample_id_input,
        subset_dir = subset_dir_input,
        sample_type = sample_type_input,
        stem_id = stem_id,
        status = "empty_after_split"
      )
      idx <- idx + 1L
      next
    }

    out_prefix <- if (stem_id == "combined") {
      file.path(subset_base, OUT_BASENAME)
    } else {
      file.path(subset_base, paste0(OUT_BASENAME, "_", stem_id))
    }

    fwrite(
      data.table(bin_id = stem_bin_ids),
      paste0(out_prefix, "_bin_ids.tsv"),
      sep = "\t"
    )

    common_genes <- intersect(rownames(mat_tenx_gb), colnames(xenium_part_bg))
    cat("Common genes:", length(common_genes), "\n")

    if (length(common_genes) < 50) {
      warning("Very few common genes in ", sample_id_input, " / ", subset_dir_input, " / ", stem_id)
    }

    if (length(common_genes) < 2) {
      res_one[[idx]] <- data.table(
        sample_id = sample_id_input,
        subset_dir = subset_dir_input,
        sample_type = sample_type_input,
        stem_id = stem_id,
        status = "too_few_common_genes",
        n_common_genes = length(common_genes)
      )
      idx <- idx + 1L
      next
    }

    mat_tenx2 <- mat_tenx_gb[common_genes, , drop = FALSE]
    xenium2_bg <- xenium_part_bg[, common_genes, drop = FALSE]
    xenium2_gb <- t(xenium2_bg)
    ref_vec <- rowMeans(xenium2_gb)

    cat(sprintf(
      "Stem %s ref summary: min=%.6f median=%.6f max=%.6f\n",
      stem_id, min(ref_vec), median(ref_vec), max(ref_vec)
    ))

    cors_pearson <- apply(mat_tenx2, 2, function(v) {
      if (sd(v) == 0 || sd(ref_vec) == 0) return(NA_real_)
      cor(v, ref_vec, method = "pearson")
    })

    cors_spearman <- apply(mat_tenx2, 2, function(v) {
      if (sd(v) == 0 || sd(ref_vec) == 0) return(NA_real_)
      suppressWarnings(cor(v, ref_vec, method = "spearman"))
    })

    cor_dt <- data.table(
      barcode = colnames(mat_tenx2),
      cor_pearson = as.numeric(cors_pearson),
      cor_spearman = as.numeric(cors_spearman)
    )

    plot_dt <- merge(umap, cor_dt, by = "barcode", all.x = TRUE)

    out_tsv <- paste0(out_prefix, ".tsv")
    pearson_png <- paste0(out_prefix, "_pearson.png")
    spearman_png <- paste0(out_prefix, "_spearman.png")

    fwrite(plot_dt, out_tsv, sep = "\t")
    cat("Wrote:", out_tsv, "\n")

    if (stem_id == "stem1" && length(split_labels) >= 1) {
      plot_label <- split_labels[1]
    } else if (stem_id == "stem2" && length(split_labels) >= 2) {
      plot_label <- split_labels[2]
    } else {
      plot_label <- sample_label_all
    }

    plot_umap_correlation_base(
      dt = plot_dt,
      value_col = "cor_pearson",
      main_text = paste0(plot_label, "\nPearson"),
      output_png = pearson_png
    )
    cat("Wrote:", pearson_png, "\n")

    plot_umap_correlation_base(
      dt = plot_dt,
      value_col = "cor_spearman",
      main_text = paste0(plot_label, "\nSpearman"),
      output_png = spearman_png
    )
    cat("Wrote:", spearman_png, "\n")

    res_one[[idx]] <- data.table(
      sample_id = sample_id_input,
      subset_dir = subset_dir_input,
      sample_type = sample_type_input,
      stem_id = stem_id,
      plot_label = plot_label,
      status = "ok",
      n_bins = nrow(xenium_part_bg),
      n_common_genes = length(common_genes),
      out_tsv = out_tsv,
      pearson_png = pearson_png,
      spearman_png = spearman_png
    )

    idx <- idx + 1L
  }

  rbindlist(res_one, fill = TRUE, use.names = TRUE)
}

# ============================================================
# Run all jobs
# ============================================================
results <- vector("list", nrow(cfg))

for (i in seq_len(nrow(cfg))) {
  results[[i]] <- tryCatch(
    run_one_subset(
      sample_id_input = cfg$sample_id[i],
      subset_dir_input = cfg$subset_dir[i],
      sample_type_input = cfg$sample_type[i]
    ),
    error = function(e) {
      warning("FAILED: ", cfg$sample_id[i], " / ", cfg$subset_dir[i], " | ", conditionMessage(e))
      data.table(
        sample_id = cfg$sample_id[i],
        subset_dir = cfg$subset_dir[i],
        sample_type = cfg$sample_type[i],
        stem_id = NA_character_,
        status = paste0("error: ", conditionMessage(e))
      )
    }
  )
}

summary_dt <- rbindlist(results, fill = TRUE, use.names = TRUE)
summary_file <- file.path(BASE_DIR, paste0(OUT_BASENAME, "_all_cambium_summary.tsv"))
fwrite(summary_dt, summary_file, sep = "\t")

cat("\n============================================================\n")
cat("DONE.\n")
cat("Summary saved to:\n")
cat(summary_file, "\n")