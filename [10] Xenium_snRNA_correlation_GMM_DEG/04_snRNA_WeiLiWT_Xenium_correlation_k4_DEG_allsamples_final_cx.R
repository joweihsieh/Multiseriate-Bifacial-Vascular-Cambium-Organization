#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})



setwd("/home/woodydrylab/FileShare/20260117_MultiLayer_Unifacial/two_cambium")

# ============================================================
# Input
# ============================================================
UMI_CSV <- "/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/all_UMI_tables/geneUMI_TenX_PtrWT2forWOX2_v701.csv"

CELL_SUMMARY_TSV <- "GMM_on_mean_median_correlation/cell_summary_with_clusters.tsv"

SYNONYM_TXT <- "/home/woodydrylab/FileShare/20260121_Xenium/Ptrichocarpa_533_v4.1.synonym_20241227.txt"
FUNC_PATH   <- "/home/woodydrylab/FileShare/20260121_Xenium/Ptrichocarpa_v4.1_complete table_v1.1.csv"

OUTDIR_BASE <- "GMM_on_mean_median_correlation/snRNA_DEG_from_mean_median_GMM_codex"
dir.create(OUTDIR_BASE, recursive = TRUE, showWarnings = FALSE)

SAMPLE_ID <- "PtrWT2forWOX2_v701"

LOGFC_CUTOFF <- 1
PADJ_CUTOFF <- 0.05
MIN_PCT <- 0.1
TEST_USE <- "wilcox"
ASSAY_USE <- "RNA"

# ============================================================
# Mapping
# ============================================================
MEAN_CPHLOEM <- "2"
MEAN_CXYLEM  <- "4"

MEDIAN_CPHLOEM <- "2"
MEDIAN_CXYLEM  <- "3"

# ============================================================
# Helpers
# ============================================================
clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- toupper(x)
  x <- sub("\\.V4\\.1$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x <- gsub("^POTRI[\\.-]", "POTRI_", x)
  x
}

dedupe_bridge_by_gene <- function(dt, value_cols = character()) {
  dt <- as.data.table(dt)
  if (!"gene_clean" %in% names(dt)) {
    stop("dedupe_bridge_by_gene() requires gene_clean column")
  }

  value_cols <- intersect(value_cols, names(dt))
  keep_cols <- unique(c("gene_clean", value_cols))
  dt <- dt[, ..keep_cols]

  if (length(value_cols) == 0) {
    return(unique(dt, by = "gene_clean"))
  }

  dt[, lapply(.SD, function(col) {
    col <- as.character(col)
    col <- col[!is.na(col) & nzchar(trimws(col))]
    if (length(col) == 0) NA_character_ else col[1]
  }), by = gene_clean, .SDcols = value_cols]
}

cell_key1 <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\.", "-", x)
  x
}

cell_key2 <- function(x) {
  x <- cell_key1(x)
  x <- sub("-1$", "", x)
  x
}

read_umi_csv_as_sparse <- function(csv_file) {
  message("[1/10] Reading UMI csv: ", csv_file)
  dt <- fread(csv_file)

  if (ncol(dt) < 2) stop("UMI csv 至少要有 2 欄")

  first_col <- names(dt)[1]
  feat_ids <- dt[[first_col]]
  expr_dt <- dt[, -1, with = FALSE]

  mat <- as.matrix(expr_dt)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- feat_ids
  colnames(mat) <- names(expr_dt)

  gene_like_rows <- mean(grepl("^Potri[._]|^AT[1-5MC]G|^LOC|^evm", rownames(mat), ignore.case = TRUE), na.rm = TRUE)
  gene_like_cols <- mean(grepl("^Potri[._]|^AT[1-5MC]G|^LOC|^evm", colnames(mat), ignore.case = TRUE), na.rm = TRUE)

  if (is.finite(gene_like_cols) && is.finite(gene_like_rows) && gene_like_cols > gene_like_rows) {
    mat <- t(mat)
  }

  raw_gene_ids <- rownames(mat)
  bridge <- data.table(
    gene_v4.1_UMI_original = raw_gene_ids,
    gene_clean = clean_gene_id(raw_gene_ids)
  )

  if (anyDuplicated(colnames(mat))) colnames(mat) <- make.unique(colnames(mat))

  rownames(mat) <- clean_gene_id(rownames(mat))
  if (anyDuplicated(rownames(mat))) {
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  }

  list(
    mat = Matrix(mat, sparse = TRUE),
    bridge = unique(bridge, by = "gene_clean")
  )
}

run_findmarkers <- function(object, ident.1, ident.2 = NULL,
                            assay = ASSAY_USE, logfc.threshold = 0,
                            min.pct = MIN_PCT, test.use = TEST_USE, only.pos = FALSE) {
  tryCatch(
    {
      FindMarkers(
        object = object,
        ident.1 = ident.1,
        ident.2 = ident.2,
        only.pos = only.pos,
        assay = assay,
        slot = "data",
        logfc.threshold = logfc.threshold,
        min.pct = min.pct,
        test.use = test.use,
        verbose = FALSE
      )
    },
    error = function(e1) {
      msg <- conditionMessage(e1)
      slot_layer_issue <- grepl("unused argument.*slot|unused argument.*layer|slot|layer", msg, ignore.case = TRUE)
      if (!slot_layer_issue) stop(e1)

      FindMarkers(
        object = object,
        ident.1 = ident.1,
        ident.2 = ident.2,
        only.pos = only.pos,
        assay = assay,
        layer = "data",
        logfc.threshold = logfc.threshold,
        min.pct = min.pct,
        test.use = test.use,
        verbose = FALSE
      )
    }
  )
}

get_norm_data <- function(object, assay = ASSAY_USE) {
  x <- tryCatch(GetAssayData(object, assay = assay, layer = "data"), error = function(e) NULL)
  if (is.null(x)) x <- GetAssayData(object, assay = assay, slot = "data")
  x
}

add_deg_metadata <- function(deg_dt, sample_id, comparison, master_bridge) {
  if (!"gene" %in% names(deg_dt)) deg_dt <- as.data.table(deg_dt, keep.rownames = "gene")
  deg_dt <- as.data.table(deg_dt)

  deg_dt[, gene_clean := clean_gene_id(gene)]
  deg_dt[, sample := sample_id]
  deg_dt[, comparison := comparison]

  deg_dt <- merge(deg_dt, master_bridge, by = "gene_clean", all.x = TRUE)

  front_cols <- c(
    "sample", "comparison", "gene", "gene_clean",
    "gene_v4.1_UMI_original", "Gene.ID", "Common.name",
    "Best.hit.Arabidopsis", "Annotated.function", "gene_v3"
  )
  front_cols <- intersect(front_cols, names(deg_dt))
  setcolorder(deg_dt, c(front_cols, setdiff(names(deg_dt), front_cols)))

  deg_dt
}

write_deg_outputs <- function(deg_dt, out_prefix, logfc_cutoff = LOGFC_CUTOFF, padj_cutoff = PADJ_CUTOFF) {
  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), names(deg_dt))[1]
  if (length(fc_col) == 0 || is.na(fc_col)) {
    deg_dt[, avg_log2FC := NA_real_]
    fc_col <- "avg_log2FC"
  }

  if (!"p_val_adj" %in% names(deg_dt)) deg_dt[, p_val_adj := NA_real_]
  if (!"p_val" %in% names(deg_dt)) deg_dt[, p_val := NA_real_]
  if (!"pct.1" %in% names(deg_dt)) deg_dt[, `pct.1` := NA_real_]
  if (!"pct.2" %in% names(deg_dt)) deg_dt[, `pct.2` := NA_real_]

  deg_all  <- copy(deg_dt)
  deg_sig  <- deg_all[abs(get(fc_col)) >= logfc_cutoff & p_val_adj < padj_cutoff]
  deg_up   <- deg_all[get(fc_col) >= logfc_cutoff & p_val_adj < padj_cutoff]
  deg_down <- deg_all[get(fc_col) <= -logfc_cutoff & p_val_adj < padj_cutoff]

  fwrite(deg_all,  paste0(out_prefix, ".tsv"),      sep = "\t")
  fwrite(deg_sig,  paste0(out_prefix, "_sig.tsv"),  sep = "\t")
  fwrite(deg_up,   paste0(out_prefix, "_up.tsv"),   sep = "\t")
  fwrite(deg_down, paste0(out_prefix, "_down.tsv"), sep = "\t")

  if ("gene_v3" %in% names(deg_sig)) {
    writeLines(na.omit(unique(deg_sig$gene_v3)),  paste0(out_prefix, "_ShinyGO_v3_sig.txt"))
    writeLines(na.omit(unique(deg_up$gene_v3)),   paste0(out_prefix, "_ShinyGO_v3_up.txt"))
    writeLines(na.omit(unique(deg_down$gene_v3)), paste0(out_prefix, "_ShinyGO_v3_down.txt"))
  }

  invisible(list(all = deg_all, sig = deg_sig, up = deg_up, down = deg_down))
}

write_opposite_regulation <- function(sig_file, out_file) {
  dt_sig <- fread(sig_file)

  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), names(dt_sig))[1]
  if (length(fc_col) == 0 || is.na(fc_col)) {
    stop("Cannot find avg_log2FC / avg_logFC in: ", sig_file)
  }

  dt_sig[, abs_fc_tmp := abs(get(fc_col))]
  setorder(dt_sig, comparison, gene_clean, -abs_fc_tmp)

  dt_sig_uniq <- dt_sig[, .SD[1], by = .(comparison, gene_clean)]

  dt_wide <- dcast(
    dt_sig_uniq,
    gene_clean ~ comparison,
    value.var = fc_col,
    fill = NA_real_
  )

  needed_cols <- c("Cphloem_vs_rest", "Cxylem_vs_rest")
  if (!all(needed_cols %in% names(dt_wide))) {
    warning("Missing required comparisons in ", sig_file, ". Skip opposite regulation.")
    return(NULL)
  }

  dt_wide_opp <- dt_wide[
    !is.na(Cphloem_vs_rest) &
    !is.na(Cxylem_vs_rest) &
    (Cphloem_vs_rest * Cxylem_vs_rest < 0)
  ]

  target_genes <- dt_wide_opp$gene_clean

  dt_sig_opp <- dt_sig[gene_clean %in% target_genes]
  dt_sig_opp[, abs_fc_tmp := NULL]
  setorder(dt_sig_opp, gene_clean, comparison)

  fwrite(dt_sig_opp, out_file, sep = "\t")

  if ("gene_v3" %in% names(dt_sig_opp)) {
    writeLines(
      na.omit(unique(dt_sig_opp$gene_v3)),
      sub("\\.tsv$", "_v3_IDs.txt", out_file)
    )
  }

  invisible(dt_sig_opp)
}

run_one_pipeline <- function(seu, meta_dt, cluster_col, cphloem_cluster, cxylem_cluster, outdir, label_name) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message("============================================================")
  message("[RUN] ", label_name)
  message("  cluster_col = ", cluster_col)
  message("  Cphloem = ", cphloem_cluster)
  message("  Cxylem  = ", cxylem_cluster)

  obj_cells <- colnames(seu)
  meta_dt2 <- copy(meta_dt)

  meta_dt2[, cell_k1 := cell_key1(cell)]
  meta_dt2[, cell_k2 := cell_key2(cell)]

  cluster_map1 <- setNames(as.character(meta_dt2[[cluster_col]]), meta_dt2$cell_k1)
  cluster_map2 <- setNames(as.character(meta_dt2[[cluster_col]]), meta_dt2$cell_k2)

  match1 <- unname(cluster_map1[cell_key1(obj_cells)])
  match2 <- unname(cluster_map2[cell_key2(obj_cells)])

  seu_tmp <- seu

  if (sum(!is.na(match2)) > sum(!is.na(match1))) {
    seu_tmp[[cluster_col]] <- match2
  } else {
    seu_tmp[[cluster_col]] <- match1
  }

  seu_tmp <- subset(seu_tmp, cells = colnames(seu_tmp)[!is.na(seu_tmp[[cluster_col]][,1])])

  raw_cluster <- as.character(seu_tmp[[cluster_col]][,1])

  new_group <- ifelse(
    raw_cluster == cphloem_cluster, "Cphloem",
    ifelse(raw_cluster == cxylem_cluster, "Cxylem", "Other")
  )

  seu_tmp$new_group <- new_group
  Idents(seu_tmp) <- "new_group"

  group_count_dt <- data.table(group = names(table(Idents(seu_tmp))), n = as.integer(table(Idents(seu_tmp))))
  fwrite(group_count_dt, file.path(outdir, "group_counts.tsv"), sep = "\t")

  saveRDS(seu_tmp, file.path(outdir, paste0("seurat_", label_name, "_grouped.rds")))

  deg_list <- list()
  target_groups <- c("Cphloem", "Cxylem")

  for (grp in target_groups) {
    if (!grp %in% levels(Idents(seu_tmp))) next

    message("  - DEG: ", grp, " vs rest")
    deg <- run_findmarkers(object = seu_tmp, ident.1 = grp, ident.2 = NULL)
    deg <- add_deg_metadata(
      deg,
      sample_id = SAMPLE_ID,
      comparison = paste0(grp, "_vs_rest"),
      master_bridge = master_bridge
    )
    write_deg_outputs(
      deg_dt = deg,
      out_prefix = file.path(outdir, paste0("DEG_", grp, "_vs_rest"))
    )
    deg_list[[grp]] <- deg
  }

  if (length(deg_list) > 0) {
    deg_all_ovr <- rbindlist(deg_list, use.names = TRUE, fill = TRUE)
    fwrite(deg_all_ovr, file.path(outdir, "DEG_Cphloem_Cxylem_one_vs_rest_all.tsv"), sep = "\t")

    fc_col_all <- intersect(c("avg_log2FC", "avg_logFC"), names(deg_all_ovr))[1]
    if (length(fc_col_all) == 0 || is.na(fc_col_all)) {
      deg_all_ovr[, avg_log2FC := NA_real_]
      fc_col_all <- "avg_log2FC"
    }
    if (!"p_val_adj" %in% names(deg_all_ovr)) {
      deg_all_ovr[, p_val_adj := NA_real_]
    }

    deg_all_ovr_sig <- deg_all_ovr[abs(get(fc_col_all)) >= LOGFC_CUTOFF & p_val_adj < PADJ_CUTOFF]
    fwrite(deg_all_ovr_sig, file.path(outdir, "DEG_Cphloem_Cxylem_one_vs_rest_sig.tsv"), sep = "\t")

    write_opposite_regulation(
      sig_file = file.path(outdir, "DEG_Cphloem_Cxylem_one_vs_rest_sig.tsv"),
      out_file = file.path(outdir, "DEG_Cphloem_Cxylem_Opposite_Regulation.tsv")
    )
  }

  if (all(target_groups %in% levels(Idents(seu_tmp)))) {
    message("  - DEG: Cphloem vs Cxylem")
    deg_direct <- run_findmarkers(object = seu_tmp, ident.1 = "Cphloem", ident.2 = "Cxylem")
    deg_direct <- add_deg_metadata(
      deg_direct,
      sample_id = SAMPLE_ID,
      comparison = "Cphloem_vs_Cxylem",
      master_bridge = master_bridge
    )
    write_deg_outputs(
      deg_dt = deg_direct,
      out_prefix = file.path(outdir, "DEG_Cphloem_vs_Cxylem")
    )
  }

  expr_data <- get_norm_data(seu_tmp, assay = ASSAY_USE)

  for (grp in target_groups) {
    if (!grp %in% levels(Idents(seu_tmp))) next

    cells_use <- WhichCells(seu_tmp, idents = grp)
    expr_sub <- expr_data[, cells_use, drop = FALSE]
    avg_expr <- Matrix::rowMeans(expr_sub)

    saveRDS(expr_sub, file.path(outdir, paste0("expr_", grp, "_logNorm.rds")))

    avg_dt <- data.table(
      gene = names(avg_expr),
      gene_clean = clean_gene_id(names(avg_expr)),
      avg_expr = as.numeric(avg_expr)
    )
    avg_dt <- merge(avg_dt, master_bridge, by = "gene_clean", all.x = TRUE)
    fwrite(avg_dt, file.path(outdir, paste0("AverageExpression_", grp, ".tsv")), sep = "\t")
  }

  if (all(target_groups %in% levels(Idents(seu_tmp)))) {
    cells_cphloem <- WhichCells(seu_tmp, idents = "Cphloem")
    cells_cxylem  <- WhichCells(seu_tmp, idents = "Cxylem")

    avg_cphloem <- Matrix::rowMeans(expr_data[, cells_cphloem, drop = FALSE])
    avg_cxylem  <- Matrix::rowMeans(expr_data[, cells_cxylem,  drop = FALSE])

    avg_dt <- data.table(
      gene = names(avg_cphloem),
      gene_clean = clean_gene_id(names(avg_cphloem)),
      avg_expr_Cphloem = as.numeric(avg_cphloem),
      avg_expr_Cxylem  = as.numeric(avg_cxylem)
    )
    avg_dt <- merge(avg_dt, master_bridge, by = "gene_clean", all.x = TRUE)
    fwrite(avg_dt, file.path(outdir, "AverageExpression_Cphloem_Cxylem.tsv"), sep = "\t")
  }

  invisible(seu_tmp)
}

# ============================================================
# 1) Read snRNA UMI table and create Seurat object
# ============================================================
umi_data <- read_umi_csv_as_sparse(UMI_CSV)
mat <- umi_data$mat
master_bridge <- umi_data$bridge

message("UMI matrix: ", nrow(mat), " genes x ", ncol(mat), " cells")

seu <- CreateSeuratObject(
  counts = mat,
  project = SAMPLE_ID,
  assay = ASSAY_USE,
  min.cells = 0,
  min.features = 0
)

DefaultAssay(seu) <- ASSAY_USE
seu <- NormalizeData(
  object = seu,
  assay = ASSAY_USE,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ============================================================
# 2) Build master bridge
# ============================================================
message("[2/10] Building ID Translation & Annotation Bridge")

func_dt <- fread(FUNC_PATH)
if (!"Gene.ID" %in% names(func_dt)) {
  stop("FUNC_PATH must contain column: Gene.ID")
}
func_dt[, gene_clean := clean_gene_id(`Gene.ID`)]
func_map <- dedupe_bridge_by_gene(
  func_dt,
  value_cols = c("Gene.ID", "Common.name", "Best.hit.Arabidopsis", "Annotated.function")
)
master_bridge <- merge(master_bridge, func_map, by = "gene_clean", all.x = TRUE)

syn_dt <- fread(SYNONYM_TXT, header = "auto", fill = TRUE)
syn_orig_col <- names(syn_dt)[1]
syn_v3_col <- names(syn_dt)[2]
syn_dt[, gene_clean := clean_gene_id(get(syn_orig_col))]
setnames(syn_dt, syn_v3_col, "gene_v3")
syn_map <- dedupe_bridge_by_gene(syn_dt, value_cols = "gene_v3")
master_bridge <- merge(master_bridge, syn_map, by = "gene_clean", all.x = TRUE)
master_bridge <- dedupe_bridge_by_gene(
  master_bridge,
  value_cols = c(
    "gene_v4.1_UMI_original", "Gene.ID", "Common.name",
    "Best.hit.Arabidopsis", "Annotated.function", "gene_v3"
  )
)

fwrite(master_bridge, file.path(OUTDIR_BASE, "Master_Gene_Bridge_and_Annotation.tsv"), sep = "\t")

# ============================================================
# 3) Read cell summary with clusters
# ============================================================
message("[3/10] Reading cell summary: ", CELL_SUMMARY_TSV)
meta_dt <- fread(CELL_SUMMARY_TSV)

needed_cols <- c("cell", "cluster_mean", "cluster_median")
missing_cols <- setdiff(needed_cols, names(meta_dt))
if (length(missing_cols) > 0) {
  stop("CELL_SUMMARY_TSV 缺少欄位: ", paste(missing_cols, collapse = ", "))
}

meta_dt[, cell := as.character(cell)]
meta_dt[, cluster_mean := as.character(cluster_mean)]
meta_dt[, cluster_median := as.character(cluster_median)]

fwrite(meta_dt, file.path(OUTDIR_BASE, "input_cell_summary_used.tsv"), sep = "\t")

# ============================================================
# 4) Export mapping tables
# ============================================================
mean_map_dt <- unique(meta_dt[, .(
  cell,
  cluster_mean,
  new_group_mean = ifelse(
    cluster_mean == MEAN_CPHLOEM, "Cphloem",
    ifelse(cluster_mean == MEAN_CXYLEM, "Cxylem", "Other")
  )
)])

median_map_dt <- unique(meta_dt[, .(
  cell,
  cluster_median,
  new_group_median = ifelse(
    cluster_median == MEDIAN_CPHLOEM, "Cphloem",
    ifelse(cluster_median == MEDIAN_CXYLEM, "Cxylem", "Other")
  )
)])

fwrite(mean_map_dt, file.path(OUTDIR_BASE, "cell_mapping_mean.tsv"), sep = "\t")
fwrite(median_map_dt, file.path(OUTDIR_BASE, "cell_mapping_median.tsv"), sep = "\t")

# ============================================================
# 5) Run mean pipeline
# ============================================================
mean_outdir <- file.path(OUTDIR_BASE, "mean_based")
dir.create(mean_outdir, recursive = TRUE, showWarnings = FALSE)

seu_mean <- run_one_pipeline(
  seu = seu,
  meta_dt = meta_dt,
  cluster_col = "cluster_mean",
  cphloem_cluster = MEAN_CPHLOEM,
  cxylem_cluster = MEAN_CXYLEM,
  outdir = mean_outdir,
  label_name = "mean"
)

# ============================================================
# 6) Run median pipeline
# ============================================================
median_outdir <- file.path(OUTDIR_BASE, "median_based")
dir.create(median_outdir, recursive = TRUE, showWarnings = FALSE)

seu_median <- run_one_pipeline(
  seu = seu,
  meta_dt = meta_dt,
  cluster_col = "cluster_median",
  cphloem_cluster = MEDIAN_CPHLOEM,
  cxylem_cluster = MEDIAN_CXYLEM,
  outdir = median_outdir,
  label_name = "median"
)

message("============================================================")
message("All done.")
message("Base output directory: ", normalizePath(OUTDIR_BASE))
