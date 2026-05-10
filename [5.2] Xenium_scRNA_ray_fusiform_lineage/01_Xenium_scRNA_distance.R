#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(Matrix)
  library(ggplot2)
})

set.seed(123)

############################################################
# user settings
############################################################
XENIUM_BASE <- "/home/woodydrylab/FileShare/20260121_Xenium"
DOMAIN_XLSX <- "/home/woodydrylab/FileShare/20260121_Xenium/k10_5domain.xlsx"

K10_DIRNAME <- "kmeans_k10_raw_out"

PLOTTING_FILE <- "/home/f06b22037/DiskArray_f06b22037/SSD2/RK/1136project_SingleCell/results/Single_species_analysis/all_plotting_tables/plotting_TenX_Ptr.csv"
UMI_FILE <- "/home/f06b22037/DiskArray_f06b22037/SSD2/RK/1136project_SingleCell/results/Single_species_analysis/all_UMI_tables/geneUMI_TenX_Ptr.csv"

RAY_CLUSTER_SC <- 3
FUS_CLUSTER_SC <- 6
N_PERM <- 20

# filtering thresholds
MIN_DET_FRAC_XENIUM <- 0.05
MIN_DET_FRAC_SC <- 0.05
VAR_QUANTILE_XENIUM <- 0
VAR_QUANTILE_SC <- 0.5

OUTDIR <- file.path(XENIUM_BASE, "xenium_sc_with_controls_nonzero_mean_euclidean")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OUT_XENIUM_ASSIGN_TSV      <- file.path(OUTDIR, "xenium_initials_assignment.tsv")
OUT_XENIUM_BINCOUNT_TSV    <- file.path(OUTDIR, "xenium_group_bin_count.tsv")
OUT_SC_EXPR_TSV            <- file.path(OUTDIR, "singlecell_ray_fusiform_organizer_gene_nonzero_mean.tsv")
OUT_SC_FILTER_TSV          <- file.path(OUTDIR, "singlecell_filter_stats.tsv")
OUT_FILTER_SUMMARY_TSV     <- file.path(OUTDIR, "filter_summary.tsv")

OUT_OBSERVED_TSV           <- file.path(OUTDIR, "observed_metrics.tsv")
OUT_PERM_TSV               <- file.path(OUTDIR, "permuted_metrics.tsv")
OUT_COMBINED_TSV           <- file.path(OUTDIR, "combined_metrics.tsv")

OUT_PLOT_EUCLIDEAN_PDF     <- file.path(OUTDIR, "observed_euclidean_boxplot.pdf")
OUT_PLOT_EUCLIDEAN_PNG     <- file.path(OUTDIR, "observed_euclidean_boxplot.png")

OUT_PLOT_MATCH_EUC_PDF     <- file.path(OUTDIR, "matched_vs_controls_euclidean_boxplot.pdf")
OUT_PLOT_MATCH_EUC_PNG     <- file.path(OUTDIR, "matched_vs_controls_euclidean_boxplot.png")

OUT_RESULT_RDS             <- file.path(OUTDIR, "xenium_sc_with_controls_nonzero_mean_euclidean_result.rds")

############################################################
# helper
############################################################
stopf <- function(...) stop(sprintf(...), call. = FALSE)

to_char <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x
}

safe_fread <- function(path, ...) {
  if (!file.exists(path)) stopf("File not found: %s", path)
  fread(path, ...)
}

safe_readRDS <- function(path) {
  if (!file.exists(path)) stopf("File not found: %s", path)
  readRDS(path)
}

find_first_col <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

clean_gene_id <- function(x) {
  x <- sub("\\.v4\\.1$", "", x)
  x <- gsub("^Potri\\.", "Potri_", x)
  x
}

get_stem_suffix <- function(sample_label) {
  if (grepl("stem1$", sample_label)) return("stem1")
  if (grepl("stem2$", sample_label)) return("stem2")
  NA_character_
}

get_meta_file <- function(k10_dir, sample_label) {
  stem_suffix <- get_stem_suffix(sample_label)

  if (!is.na(stem_suffix)) {
    f <- file.path(k10_dir, paste0("bin_metadata_with_cluster_raw_", stem_suffix, ".tsv"))
    if (!file.exists(f)) {
      stopf("Expected stem-specific metadata not found for %s: %s", sample_label, f)
    }
    return(f)
  }

  f <- file.path(k10_dir, "bin_metadata_with_cluster_raw.tsv")
  if (!file.exists(f)) {
    stopf("Metadata file not found for %s: %s", sample_label, f)
  }
  f
}

get_count_file <- function(grid_dir) {
  f <- file.path(grid_dir, "counts_bins_by_genes_sparse.rds")
  if (!file.exists(f)) {
    stopf("Shared count matrix not found in %s", grid_dir)
  }
  f
}

orient_counts_matrix <- function(mat, meta_bin_ids) {
  if (!inherits(mat, "Matrix")) {
    mat <- as(mat, "dgCMatrix")
  }

  rn <- rownames(mat)
  cn <- colnames(mat)

  row_hit <- if (!is.null(rn)) sum(rn %in% meta_bin_ids) else 0L
  col_hit <- if (!is.null(cn)) sum(cn %in% meta_bin_ids) else 0L

  if (col_hit >= row_hit) {
    return(mat)
  } else {
    message("  counts matrix appears to be bins x genes, transposing...")
    return(as(t(mat), "dgCMatrix"))
  }
}

# non-zero mean
calc_row_nonzero_mean <- function(mat_sub) {
  if (!inherits(mat_sub, "Matrix")) {
    mat_sub <- as(mat_sub, "dgCMatrix")
  }

  rs <- Matrix::rowSums(mat_sub)
  nnz <- Matrix::rowSums(mat_sub != 0)

  out <- ifelse(nnz > 0, rs / nnz, 0)
  names(out) <- rownames(mat_sub)
  out
}

calc_row_detect_frac <- function(mat_sub) {
  if (!inherits(mat_sub, "Matrix")) {
    mat_sub <- as(mat_sub, "dgCMatrix")
  }

  out <- Matrix::rowMeans(mat_sub != 0)
  names(out) <- rownames(mat_sub)
  out
}

calc_row_log1p_var <- function(mat_sub) {
  if (!inherits(mat_sub, "Matrix")) {
    mat_sub <- as(mat_sub, "dgCMatrix")
  }

  n <- ncol(mat_sub)
  out <- numeric(nrow(mat_sub))

  if (n <= 1) {
    names(out) <- rownames(mat_sub)
    return(out)
  }

  m1 <- mat_sub
  m1@x <- log1p(m1@x)

  mu <- Matrix::rowMeans(m1)

  m2 <- m1
  m2@x <- m2@x ^ 2
  mu2 <- Matrix::rowMeans(m2)

  out <- (mu2 - mu^2) * n / (n - 1)
  out[out < 0] <- 0
  names(out) <- rownames(mat_sub)
  out
}

calc_euclidean <- function(x, y) {
  sqrt(sum((x - y)^2))
}

is_barcode_like <- function(x) {
  x <- x[!is.na(x)]
  x <- x[x != ""]
  if (length(x) == 0) return(FALSE)
  mean(grepl("^[A-Za-z0-9_-]+$", x)) > 0.8 &&
    mean(nchar(x) >= 8) > 0.8
}

is_gene_like <- function(x) {
  x <- x[!is.na(x)]
  x <- x[x != ""]
  if (length(x) == 0) return(FALSE)
  mean(
    grepl("^Potri[._]", x) |
    grepl("^[A-Za-z]{2,}[0-9]{3,}", x) |
    grepl("^[A-Za-z0-9_.-]+\\.v[0-9]", x)
  ) > 0.3
}

aggregate_gene_expr_table <- function(dt, gene_col = "gene") {
  if (!anyDuplicated(dt[[gene_col]])) return(dt)
  num_cols <- setdiff(names(dt), gene_col)
  dt[, lapply(.SD, mean), by = gene, .SDcols = num_cols]
}

aggregate_gene_stats_table <- function(dt, gene_col = "gene") {
  if (!anyDuplicated(dt[[gene_col]])) return(dt)

  dt[, .(
    detect_frac = max(detect_frac, na.rm = TRUE),
    var_log1p = max(var_log1p, na.rm = TRUE)
  ), by = gene]
}

build_gene_stats <- function(mat_sub) {
  dt <- data.table(
    gene = clean_gene_id(rownames(mat_sub)),
    detect_frac = as.numeric(calc_row_detect_frac(mat_sub)),
    var_log1p = as.numeric(calc_row_log1p_var(mat_sub))
  )
  aggregate_gene_stats_table(dt)
}

get_keep_genes_from_stats <- function(stats_dt, min_det_frac, var_quantile) {
  dt <- copy(stats_dt)

  dt <- dt[!is.na(gene) & gene != ""]
  dt <- dt[!is.na(detect_frac) & !is.na(var_log1p)]

  if (nrow(dt) == 0) return(character(0))

  det_keep <- dt[detect_frac >= min_det_frac]
  if (nrow(det_keep) == 0) return(character(0))

  var_thr <- as.numeric(quantile(det_keep$var_log1p, probs = var_quantile, na.rm = TRUE, names = FALSE))

  keep <- det_keep[var_log1p >= var_thr, gene]
  unique(keep)
}

make_common_vectors <- function(xen_dt, sc_ray_vec, sc_fus_vec, xen_col, keep_genes) {
  common_genes <- Reduce(intersect, list(
    xen_dt$gene,
    names(sc_ray_vec),
    names(sc_fus_vec),
    keep_genes
  ))

  if (length(common_genes) == 0) return(NULL)

  xen_vec <- xen_dt[match(common_genes, gene), get(xen_col)]
  ray_vec <- sc_ray_vec[common_genes]
  fus_vec <- sc_fus_vec[common_genes]

  names(xen_vec) <- common_genes
  names(ray_vec) <- common_genes
  names(fus_vec) <- common_genes

  keep <- (xen_vec + ray_vec + fus_vec) > 0
  xen_vec <- xen_vec[keep]
  ray_vec <- ray_vec[keep]
  fus_vec <- fus_vec[keep]

  if (length(xen_vec) == 0) return(NULL)

  list(
    xen = log1p(xen_vec),
    ray = log1p(ray_vec),
    fus = log1p(fus_vec),
    genes = names(xen_vec)
  )
}

expand_combined_labels <- function(domain_dt) {
  out_list <- vector("list", nrow(domain_dt))
  idx <- 1L

  for (i in seq_len(nrow(domain_dt))) {
    one <- copy(domain_dt[i])
    labs <- trimws(unlist(strsplit(one$label, ",", fixed = TRUE)))
    labs <- labs[labs != ""]

    if (length(labs) == 0) {
      out_list[[idx]] <- one
      idx <- idx + 1L
    } else {
      for (lab in labs) {
        tmp <- copy(one)
        tmp[, label := lab]
        out_list[[idx]] <- tmp
        idx <- idx + 1L
      }
    }
  }

  rbindlist(out_list[seq_len(idx - 1L)], fill = TRUE)
}

save_metric_boxplot <- function(dt, metric_col, outfile_pdf, outfile_png, title_txt, ylab_txt, width = 13, height = 5) {
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

save_match_control_boxplot <- function(dt, metric_col, outfile_pdf, outfile_png, title_txt, ylab_txt, width = 9, height = 5) {
  p <- ggplot(dt, aes(x = box_group, y = get(metric_col))) +
    geom_boxplot(width = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.4) +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank()) +
    labs(
      x = NULL,
      y = ylab_txt,
      title = title_txt
    )

  ggsave(outfile_pdf, p, width = width, height = height)
  ggsave(outfile_png, p, width = width, height = height, dpi = 300)
}

############################################################
# Step 1. read and expand domain xlsx
############################################################
domain_dt <- as.data.table(read_excel(DOMAIN_XLSX))
setnames(domain_dt, names(domain_dt), make.names(names(domain_dt), unique = TRUE))

req_cols <- c("sample_id", "label", "um", "k", "cluster", "domain")
miss_cols <- setdiff(req_cols, names(domain_dt))
if (length(miss_cols) > 0) {
  stopf("Missing columns in domain xlsx: %s", paste(miss_cols, collapse = ", "))
}

domain_dt[, sample_id := to_char(sample_id)]
domain_dt[, label    := to_char(label)]
domain_dt[, um       := to_char(um)]
domain_dt[, cluster  := as.character(cluster)]
domain_dt[, domain   := to_char(domain)]
domain_dt[, k        := as.integer(k)]

domain_dt <- domain_dt[
  !is.na(sample_id) &
  !is.na(label) &
  !is.na(um) &
  !is.na(cluster) &
  !is.na(domain)
]

domain_dt <- domain_dt[k == 10]
domain_dt <- expand_combined_labels(domain_dt)

if (nrow(domain_dt) == 0) {
  stopf("No K=10 entries found in %s", DOMAIN_XLSX)
}

############################################################
# Step 2. read single-cell plotting table
############################################################
plot_dt <- safe_fread(PLOTTING_FILE)

barcode_col <- find_first_col(
  plot_dt,
  c("barcode", "Barcode", "cell", "cell_id", "Cell", "X")
)
cluster_col <- find_first_col(
  plot_dt,
  c("cluster", "Cluster", "seurat_clusters", "celltype", "group")
)

if (is.na(barcode_col)) stopf("Cannot find barcode column in plotting file.")
if (is.na(cluster_col)) stopf("Cannot find cluster column in plotting file.")

setnames(plot_dt, barcode_col, "barcode")
setnames(plot_dt, cluster_col, "cluster")

plot_dt[, barcode := as.character(barcode)]
plot_dt[, cluster := as.character(cluster)]

ray_bc_sc <- plot_dt[cluster == as.character(RAY_CLUSTER_SC), barcode]
fus_bc_sc <- plot_dt[cluster == as.character(FUS_CLUSTER_SC), barcode]

if (length(ray_bc_sc) == 0) stopf("No cells found for single-cell cluster %d", RAY_CLUSTER_SC)
if (length(fus_bc_sc) == 0) stopf("No cells found for single-cell cluster %d", FUS_CLUSTER_SC)

n_ray_sc <- length(ray_bc_sc)
n_fus_sc <- length(fus_bc_sc)

############################################################
# Step 3. read single-cell UMI and compute organizer non-zero means
############################################################
umi_dt_raw <- safe_fread(UMI_FILE)

first_col <- names(umi_dt_raw)[1]
first_vals <- as.character(umi_dt_raw[[first_col]])
first_vals <- first_vals[!is.na(first_vals) & first_vals != ""]
other_cols <- names(umi_dt_raw)[-1]

firstcol_barcode_like <- is_barcode_like(head(first_vals, 200))
firstcol_gene_like    <- is_gene_like(head(first_vals, 200))
colnames_barcode_like <- is_barcode_like(head(other_cols, 200))
colnames_gene_like    <- is_gene_like(head(other_cols, 200))

if (firstcol_barcode_like && colnames_gene_like) {
  setnames(umi_dt_raw, first_col, "barcode")
  umi_dt_raw[, barcode := as.character(barcode)]

  keep_bc <- intersect(plot_dt$barcode, umi_dt_raw$barcode)
  if (length(keep_bc) == 0) {
    stopf("No overlap between plotting barcodes and UMI table barcode rows.")
  }

  umi_dt_sub <- umi_dt_raw[barcode %in% keep_bc]

  barcode_vec <- umi_dt_sub$barcode
  gene_names <- clean_gene_id(names(umi_dt_sub)[names(umi_dt_sub) != "barcode"])

  expr_mat <- as.matrix(umi_dt_sub[, !"barcode"])
  rownames(expr_mat) <- barcode_vec
  colnames(expr_mat) <- gene_names

  umi_mat <- t(expr_mat)

  if (anyDuplicated(rownames(umi_mat))) {
    tmp_dt <- as.data.table(umi_mat, keep.rownames = "gene")
    tmp_dt <- tmp_dt[, lapply(.SD, mean), by = gene]
    gene_names2 <- tmp_dt$gene
    umi_mat <- as.matrix(tmp_dt[, !"gene"])
    rownames(umi_mat) <- gene_names2
  }

} else if (firstcol_gene_like && colnames_barcode_like) {
  setnames(umi_dt_raw, first_col, "gene")
  umi_dt_raw[, gene := as.character(gene)]
  umi_dt_raw[, gene := clean_gene_id(gene)]
  umi_dt_raw <- umi_dt_raw[!is.na(gene) & gene != ""]

  if (anyDuplicated(umi_dt_raw$gene)) {
    umi_dt_raw <- umi_dt_raw[, lapply(.SD, mean), by = gene]
  }

  keep_bc <- intersect(plot_dt$barcode, names(umi_dt_raw))
  if (length(keep_bc) == 0) {
    stopf("No overlap between plotting barcodes and UMI table columns.")
  }

  umi_mat <- as.matrix(umi_dt_raw[, c("gene", keep_bc), with = FALSE][, !"gene"])
  rownames(umi_mat) <- umi_dt_raw$gene

} else {
  stopf("Cannot confidently determine UMI table orientation.")
}

keep_ray_bc_sc <- intersect(ray_bc_sc, colnames(umi_mat))
keep_fus_bc_sc <- intersect(fus_bc_sc, colnames(umi_mat))

if (length(keep_ray_bc_sc) == 0) stopf("No overlap for single-cell ray organizer barcodes.")
if (length(keep_fus_bc_sc) == 0) stopf("No overlap for single-cell fusiform organizer barcodes.")

ray_org_expr <- calc_row_nonzero_mean(umi_mat[, keep_ray_bc_sc, drop = FALSE])
fus_org_expr <- calc_row_nonzero_mean(umi_mat[, keep_fus_bc_sc, drop = FALSE])

sc_expr_dt <- data.table(
  gene = names(ray_org_expr),
  ray_organizer_nonzero_mean = as.numeric(ray_org_expr),
  fusiform_organizer_nonzero_mean = as.numeric(fus_org_expr)
)
fwrite(sc_expr_dt, OUT_SC_EXPR_TSV, sep = "\t")

sc_stats_dt <- build_gene_stats(umi_mat)
sc_keep_genes <- get_keep_genes_from_stats(
  sc_stats_dt,
  min_det_frac = MIN_DET_FRAC_SC,
  var_quantile = VAR_QUANTILE_SC
)

sc_stats_out <- copy(sc_stats_dt)
sc_stats_out[, keep_sc := gene %in% sc_keep_genes]
fwrite(sc_stats_out, OUT_SC_FILTER_TSV, sep = "\t")

all_sc_barcodes <- colnames(umi_mat)

############################################################
# Step 4. process each Xenium sample
############################################################
sample_labels <- unique(domain_dt$label)

assign_list <- list()
bincount_list <- list()
observed_list <- list()
perm_list <- list()
per_sample_xenium_tables <- list()
filter_summary_list <- list()

for (lab in sample_labels) {
  message("Processing Xenium sample: ", lab)

  one_domain <- domain_dt[label == lab]
  if (nrow(one_domain) == 0) next

  sid <- unique(one_domain$sample_id)
  um_dir <- unique(one_domain$um)

  if (length(sid) != 1) stopf("Label %s maps to multiple sample_id values.", lab)
  if (length(um_dir) != 1) stopf("Label %s maps to multiple um values.", lab)

  cambium_clusters   <- unique(one_domain[domain == "cambium",   as.character(cluster)])
  phloem_clusters    <- unique(one_domain[domain == "phloem",    as.character(cluster)])
  epidermis_clusters <- unique(one_domain[domain == "epidermis", as.character(cluster)])

  if (length(cambium_clusters) < 2) next
  if (length(phloem_clusters) < 1) next
  if (length(epidermis_clusters) < 1) next

  grid_dir <- file.path(XENIUM_BASE, sid, um_dir)
  k10_dir  <- file.path(grid_dir, K10_DIRNAME)

  meta_path  <- get_meta_file(k10_dir, lab)
  count_path <- get_count_file(grid_dir)

  meta <- safe_fread(meta_path)

  bin_id_col <- find_first_col(meta, c("bin_id", "barcode", "id"))
  if (is.na(bin_id_col)) stopf("No bin ID column found in %s", meta_path)
  setnames(meta, bin_id_col, "bin_id")

  cluster_raw_col <- find_first_col(meta, c("cluster_raw", "cluster", "seurat_clusters"))
  if (is.na(cluster_raw_col)) stopf("No cluster column found in %s", meta_path)
  setnames(meta, cluster_raw_col, "cluster_raw")

  meta[, bin_id := as.character(bin_id)]
  meta[, cluster_raw := as.character(cluster_raw)]

  cambium_meta   <- meta[cluster_raw %in% cambium_clusters]
  phloem_meta    <- meta[cluster_raw %in% phloem_clusters]
  epidermis_meta <- meta[cluster_raw %in% epidermis_clusters]

  if (nrow(cambium_meta) == 0 || nrow(phloem_meta) == 0 || nrow(epidermis_meta) == 0) next

  count_dt <- cambium_meta[, .N, by = .(cluster_raw)]
  setorder(count_dt, -N, cluster_raw)
  count_dt <- count_dt[1:min(.N, 2)]
  if (nrow(count_dt) < 2) next

  count_dt[, xenium_group := c("fusiform initials", "ray initials")]
  count_dt[, label := lab]
  count_dt[, sample_id := sid]
  count_dt[, um := um_dir]
  assign_list[[lab]] <- copy(count_dt)

  cambium_meta <- merge(
    cambium_meta,
    count_dt[, .(cluster_raw, xenium_group)],
    by = "cluster_raw",
    all.x = FALSE
  )

  bincount_dt <- rbindlist(list(
    cambium_meta[, .N, by = .(xenium_group)],
    data.table(xenium_group = "phloem", N = nrow(phloem_meta)),
    data.table(xenium_group = "epidermis", N = nrow(epidermis_meta))
  ), fill = TRUE)
  bincount_dt[, label := lab]
  bincount_dt[, sample_id := sid]
  bincount_dt[, um := um_dir]
  bincount_list[[lab]] <- copy(bincount_dt)

  ray_bins <- cambium_meta[xenium_group == "ray initials", bin_id]
  fus_bins <- cambium_meta[xenium_group == "fusiform initials", bin_id]
  phl_bins <- phloem_meta$bin_id
  epi_bins <- epidermis_meta$bin_id

  mat <- safe_readRDS(count_path)
  mat <- orient_counts_matrix(mat, meta$bin_id)

  if (is.null(colnames(mat))) stopf("Counts matrix has no colnames: %s", count_path)
  if (is.null(rownames(mat))) stopf("Counts matrix has no rownames: %s", count_path)

  this_stem_bins <- unique(meta$bin_id)
  keep_stem_bins <- intersect(this_stem_bins, colnames(mat))
  if (length(keep_stem_bins) == 0) {
    stopf("No overlapping bin IDs between stem-specific metadata and shared counts matrix for %s", lab)
  }

  mat <- mat[, keep_stem_bins, drop = FALSE]

  xen_stats_dt <- build_gene_stats(mat)
  xen_keep_genes <- get_keep_genes_from_stats(
    xen_stats_dt,
    min_det_frac = MIN_DET_FRAC_XENIUM,
    var_quantile = VAR_QUANTILE_XENIUM
  )

  keep_genes_sample <- intersect(xen_keep_genes, sc_keep_genes)

  filter_summary_list[[lab]] <- data.table(
    label = lab,
    sample_id = sid,
    n_genes_total_sc = nrow(sc_stats_dt),
    n_genes_total_xen = nrow(xen_stats_dt),
    n_sc_keep_genes = length(sc_keep_genes),
    n_xen_keep_genes = length(xen_keep_genes),
    n_keep_genes_intersection = length(keep_genes_sample),
    pct_sc_kept = length(sc_keep_genes) / nrow(sc_stats_dt),
    pct_xen_kept = length(xen_keep_genes) / nrow(xen_stats_dt),
    pct_intersection_vs_sc = length(keep_genes_sample) / nrow(sc_stats_dt),
    pct_intersection_vs_xen = length(keep_genes_sample) / nrow(xen_stats_dt)
  )

  xen_stats_out <- copy(xen_stats_dt)
  xen_stats_out[, keep_xen := gene %in% xen_keep_genes]
  xen_stats_out[, keep_both := gene %in% keep_genes_sample]
  fwrite(
    xen_stats_out,
    file.path(OUTDIR, paste0(lab, "__xenium_filter_stats.tsv")),
    sep = "\t"
  )

  keep_ray_bins <- intersect(ray_bins, colnames(mat))
  keep_fus_bins <- intersect(fus_bins, colnames(mat))
  keep_phl_bins <- intersect(phl_bins, colnames(mat))
  keep_epi_bins <- intersect(epi_bins, colnames(mat))

  if (length(keep_ray_bins) == 0 ||
      length(keep_fus_bins) == 0 ||
      length(keep_phl_bins) == 0 ||
      length(keep_epi_bins) == 0) next

  ray_init_expr  <- calc_row_nonzero_mean(mat[, keep_ray_bins, drop = FALSE])
  fus_init_expr  <- calc_row_nonzero_mean(mat[, keep_fus_bins, drop = FALSE])
  phloem_expr    <- calc_row_nonzero_mean(mat[, keep_phl_bins, drop = FALSE])
  epidermis_expr <- calc_row_nonzero_mean(mat[, keep_epi_bins, drop = FALSE])

  xen_dt <- data.table(
    gene = clean_gene_id(names(ray_init_expr)),
    ray_initials_nonzero_mean = as.numeric(ray_init_expr),
    fusiform_initials_nonzero_mean = as.numeric(fus_init_expr),
    phloem_nonzero_mean = as.numeric(phloem_expr),
    epidermis_nonzero_mean = as.numeric(epidermis_expr)
  )

  xen_dt <- aggregate_gene_expr_table(xen_dt)
  per_sample_xenium_tables[[lab]] <- copy(xen_dt)

  fwrite(
    xen_dt,
    file.path(OUTDIR, paste0(lab, "__xenium_groups_gene_nonzero_mean.tsv")),
    sep = "\t"
  )

  group_specs <- list(
    list(xcol = "ray_initials_nonzero_mean",      xname = "ray initials"),
    list(xcol = "fusiform_initials_nonzero_mean", xname = "fusiform initials"),
    list(xcol = "phloem_nonzero_mean",            xname = "phloem"),
    list(xcol = "epidermis_nonzero_mean",         xname = "epidermis")
  )

  obs_rows <- list()
  for (gs in group_specs) {
    vv <- make_common_vectors(
      xen_dt = xen_dt,
      sc_ray_vec = ray_org_expr,
      sc_fus_vec = fus_org_expr,
      xen_col = gs$xcol,
      keep_genes = keep_genes_sample
    )
    if (is.null(vv)) next

    euc_ray <- calc_euclidean(vv$xen, vv$ray)
    euc_fus <- calc_euclidean(vv$xen, vv$fus)

    obs_rows[[length(obs_rows) + 1]] <- data.table(
      label = lab,
      sample_id = sid,
      xenium_group = gs$xname,
      sc_group = "ray organizer",
      control_type = fifelse(
        gs$xname == "phloem", "phloem_control",
        fifelse(gs$xname == "epidermis", "epidermis_control", "observed")
      ),
      euclidean_distance = euc_ray,
      n_common_genes = length(vv$xen)
    )
    obs_rows[[length(obs_rows) + 1]] <- data.table(
      label = lab,
      sample_id = sid,
      xenium_group = gs$xname,
      sc_group = "fusiform organizer",
      control_type = fifelse(
        gs$xname == "phloem", "phloem_control",
        fifelse(gs$xname == "epidermis", "epidermis_control", "observed")
      ),
      euclidean_distance = euc_fus,
      n_common_genes = length(vv$xen)
    )
  }
  if (length(obs_rows) > 0) {
    observed_list[[lab]] <- rbindlist(obs_rows, fill = TRUE)
  }

  perm_rows <- list()
  for (perm_i in seq_len(N_PERM)) {
    sampled_ray <- sample(all_sc_barcodes, size = n_ray_sc, replace = FALSE)
    remaining <- setdiff(all_sc_barcodes, sampled_ray)
    sampled_fus <- sample(remaining, size = n_fus_sc, replace = FALSE)

    perm_ray_expr <- calc_row_nonzero_mean(umi_mat[, sampled_ray, drop = FALSE])
    perm_fus_expr <- calc_row_nonzero_mean(umi_mat[, sampled_fus, drop = FALSE])

    for (gs in group_specs) {
      vv <- make_common_vectors(
        xen_dt = xen_dt,
        sc_ray_vec = perm_ray_expr,
        sc_fus_vec = perm_fus_expr,
        xen_col = gs$xcol,
        keep_genes = keep_genes_sample
      )
      if (is.null(vv)) next

      euc_ray <- calc_euclidean(vv$xen, vv$ray)
      euc_fus <- calc_euclidean(vv$xen, vv$fus)

      perm_rows[[length(perm_rows) + 1]] <- data.table(
        label = lab,
        sample_id = sid,
        perm_id = perm_i,
        xenium_group = gs$xname,
        sc_group = "shuffled ray organizer",
        control_type = "shuffled_sc_control",
        euclidean_distance = euc_ray,
        n_common_genes = length(vv$xen)
      )
      perm_rows[[length(perm_rows) + 1]] <- data.table(
        label = lab,
        sample_id = sid,
        perm_id = perm_i,
        xenium_group = gs$xname,
        sc_group = "shuffled fusiform organizer",
        control_type = "shuffled_sc_control",
        euclidean_distance = euc_fus,
        n_common_genes = length(vv$xen)
      )
    }
  }
  if (length(perm_rows) > 0) {
    perm_list[[lab]] <- rbindlist(perm_rows, fill = TRUE)
  }
}

############################################################
# Step 5. save summary tables
############################################################
assign_dt   <- rbindlist(assign_list, fill = TRUE)
bincount_dt <- rbindlist(bincount_list, fill = TRUE)
observed_dt <- rbindlist(observed_list, fill = TRUE)
perm_dt     <- rbindlist(perm_list, fill = TRUE)
filter_summary_dt <- rbindlist(filter_summary_list, fill = TRUE)

if (nrow(assign_dt) == 0) stopf("No valid Xenium samples processed.")
if (nrow(observed_dt) == 0) stopf("No observed metric results generated.")
if (nrow(perm_dt) == 0) stopf("No permuted metric results generated.")

combined_dt <- rbindlist(list(
  observed_dt[, !"perm_id"],
  perm_dt
), fill = TRUE)

fwrite(assign_dt,         OUT_XENIUM_ASSIGN_TSV,   sep = "\t")
fwrite(bincount_dt,       OUT_XENIUM_BINCOUNT_TSV, sep = "\t")
fwrite(filter_summary_dt, OUT_FILTER_SUMMARY_TSV,  sep = "\t")
fwrite(observed_dt,       OUT_OBSERVED_TSV,        sep = "\t")
fwrite(perm_dt,           OUT_PERM_TSV,            sep = "\t")
fwrite(combined_dt,       OUT_COMBINED_TSV,        sep = "\t")

############################################################
# Step 6. plots
############################################################
obs_plot_dt <- copy(observed_dt)
obs_plot_dt[, comparison := paste(xenium_group, "vs", sc_group)]
obs_levels <- c(
  "ray initials vs ray organizer",
  "ray initials vs fusiform organizer",
  "fusiform initials vs ray organizer",
  "fusiform initials vs fusiform organizer",
  "phloem vs ray organizer",
  "phloem vs fusiform organizer",
  "epidermis vs ray organizer",
  "epidermis vs fusiform organizer"
)
obs_plot_dt[, comparison := factor(comparison, levels = obs_levels)]

save_metric_boxplot(
  obs_plot_dt, "euclidean_distance",
  OUT_PLOT_EUCLIDEAN_PDF, OUT_PLOT_EUCLIDEAN_PNG,
  "Observed Xenium vs single-cell comparisons", "Euclidean distance"
)

match_dt <- rbindlist(list(
  observed_dt[
    xenium_group == "ray initials" & sc_group == "ray organizer",
    .(label, sample_id, box_group = "matched", euclidean_distance)
  ],
  observed_dt[
    xenium_group == "fusiform initials" & sc_group == "fusiform organizer",
    .(label, sample_id, box_group = "matched", euclidean_distance)
  ],
  observed_dt[
    xenium_group == "phloem",
    .(label, sample_id, box_group = "phloem control", euclidean_distance)
  ],
  observed_dt[
    xenium_group == "epidermis",
    .(label, sample_id, box_group = "epidermis control", euclidean_distance)
  ],
  perm_dt[
    xenium_group %in% c("ray initials", "fusiform initials"),
    .(label, sample_id, perm_id, box_group = "shuffled sc control", euclidean_distance)
  ]
), fill = TRUE)

match_dt[, box_group := factor(
  box_group,
  levels = c("matched", "phloem control", "epidermis control", "shuffled sc control")
)]

save_match_control_boxplot(
  match_dt, "euclidean_distance",
  OUT_PLOT_MATCH_EUC_PDF, OUT_PLOT_MATCH_EUC_PNG,
  "Matched comparisons versus control groups", "Euclidean distance"
)

############################################################
# Step 7. save RDS
############################################################
saveRDS(
  list(
    xenium_assignment = assign_dt,
    xenium_bin_count = bincount_dt,
    sc_nonzero_mean = sc_expr_dt,
    sc_filter_stats = sc_stats_dt,
    filter_summary = filter_summary_dt,
    observed = observed_dt,
    permuted = perm_dt,
    combined = combined_dt,
    per_sample_xenium_tables = per_sample_xenium_tables
  ),
  OUT_RESULT_RDS
)

message("Done.")
message("Output directory: ", OUTDIR)


############################################################
# Step 6A. statistics for observed euclidean boxplot
############################################################
obs_plot_dt <- copy(observed_dt)
obs_plot_dt[, comparison := paste(xenium_group, "vs", sc_group)]

obs_levels <- c(
  "ray initials vs ray organizer",
  "ray initials vs fusiform organizer",
  "fusiform initials vs fusiform organizer",
  "fusiform initials vs ray organizer"
)
obs_plot_dt[, comparison := factor(comparison, levels = obs_levels)]

# summary
obs_euc_summary <- obs_plot_dt[, .(
  n = .N,
  mean = mean(euclidean_distance, na.rm = TRUE),
  median = median(euclidean_distance, na.rm = TRUE),
  sd = sd(euclidean_distance, na.rm = TRUE)
), by = comparison]

fwrite(
  obs_euc_summary,
  file.path(OUTDIR, "observed_euclidean_group_summary.tsv"),
  sep = "\t"
)

# overall test
kruskal_res <- kruskal.test(euclidean_distance ~ comparison, data = obs_plot_dt)

overall_dt <- data.table(
  test = "kruskal_wallis",
  statistic = as.numeric(kruskal_res$statistic),
  df = as.numeric(kruskal_res$parameter),
  p_value = kruskal_res$p.value
)

fwrite(
  overall_dt,
  file.path(OUTDIR, "observed_euclidean_kruskal.tsv"),
  sep = "\t"
)

# pairwise test
pairwise_res <- pairwise.wilcox.test(
  x = obs_plot_dt$euclidean_distance,
  g = obs_plot_dt$comparison,
  p.adjust.method = "BH",
  exact = FALSE
)

pairwise_mat <- as.data.table(as.table(pairwise_res$p.value))
setnames(pairwise_mat, c("group1", "group2", "p_adj"))
pairwise_mat <- pairwise_mat[!is.na(p_adj)]

fwrite(
  pairwise_mat,
  file.path(OUTDIR, "observed_euclidean_pairwise_wilcox.tsv"),
  sep = "\t"
)

OUT_PLOT_EUCLIDEAN_PDF_SIMPLE     <- file.path(OUTDIR, "observed_euclidean_boxplot_simple.pdf")
OUT_PLOT_EUCLIDEAN_PNG_SIMPLE     <- file.path(OUTDIR, "observed_euclidean_boxplot_simple.png")


save_metric_boxplot2 <- function(dt, metric_col, outfile_pdf, outfile_png, title_txt, ylab_txt, width = 5, height = 7) {
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


obs_plot_dt[, comparison := factor(comparison, levels = obs_levels)]

save_metric_boxplot2(
  obs_plot_dt[obs_plot_dt$comparison%in%obs_levels,], "euclidean_distance",
  OUT_PLOT_EUCLIDEAN_PDF_SIMPLE, OUT_PLOT_EUCLIDEAN_PNG_SIMPLE,
  "", "Euclidean distance"
)
