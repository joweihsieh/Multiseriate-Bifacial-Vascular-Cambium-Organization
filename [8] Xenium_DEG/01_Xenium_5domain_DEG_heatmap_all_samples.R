#!/usr/bin/env Rscript
# date: 20260409

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
  library(pheatmap)
  library(readxl)
})

# ============================================================
# INPUTS
# ============================================================
BASE_DIR <- "/home/woodydrylab/FileShare/20260121_Xenium"
MAP_XLSX <- file.path(BASE_DIR, "k10_5domain.xlsx")

OUT_DIR  <- file.path(BASE_DIR, "Merged_DEG_k10_5domains_all")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

UM_USE <- "grid05um_out"
K_USE  <- 10

#DOMAIN_ORDER <- c("cambium", "phloem", "epidermis", "parenchyma", "sclerenchyma")
DOMAIN_ORDER <- c("epidermis", "phloem", "cambium", "sclerenchyma", "parenchyma")

TOP_N <- 300

# ============================================================
# Heatmap colors
# ============================================================
heat_colors <- colorRampPalette(
  c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B")
)(100)
heat_breaks <- seq(-3, 3, length.out = 101)

domain_colors <- c(
  cambium      = "#F8766D",
  phloem       = "#00B0F6",
  epidermis    = "#9590FF",
  parenchyma   = "#39B600",
  sclerenchyma = "#D89000"


)

# ============================================================
# Helpers
# ============================================================

expand_mapping_table <- function(map_dt) {
  map_dt <- as.data.table(copy(map_dt))

  if (!("part" %in% colnames(map_dt))) {
    map_dt[, part := NA_character_]
  }

  map_dt[, sample_id := as.character(sample_id)]
  map_dt[, label     := trimws(as.character(label))]
  map_dt[, um        := as.character(um)]
  map_dt[, k         := as.integer(k)]
  map_dt[, cluster   := as.integer(cluster)]
  map_dt[, domain    := as.character(domain)]
  map_dt[, part      := trimws(as.character(part))]

  map_dt[part %in% c("", "NA", "NaN"), part := NA_character_]
  map_dt[label %in% c("", "NA", "NaN"), label := NA_character_]

  out_list <- vector("list", nrow(map_dt) * 2L)
  idx <- 1L

  for (i in seq_len(nrow(map_dt))) {
    row <- copy(map_dt[i])

    # If the Excel file already contains a part value, keep it as is
    if (!is.na(row$part) && row$part != "") {
      out_list[[idx]] <- row
      idx <- idx + 1L
      next
    }

    # If the label contains a comma, split it into stem1 / stem2
    if (!is.na(row$label) && grepl(",", row$label, fixed = TRUE)) {
      labs <- trimws(unlist(strsplit(row$label, ",", fixed = TRUE)))

      for (lab in labs) {
        new_row <- copy(row)

        if (grepl("stem1$", lab, ignore.case = TRUE)) {
          new_row[, part := "stem1"]
        } else if (grepl("stem2$", lab, ignore.case = TRUE)) {
          new_row[, part := "stem2"]
        } else {
          stop(
            "Cannot infer part from label: ", lab,
            " | sample_id = ", row$sample_id
          )
        }

        new_row[, label := lab]
        out_list[[idx]] <- new_row
        idx <- idx + 1L
      }

    } else {
      # For a single label, treat it as allbins by default
      row[, part := "allbins"]
      out_list[[idx]] <- row
      idx <- idx + 1L
    }
  }

  out <- rbindlist(out_list[seq_len(idx - 1L)], use.names = TRUE, fill = TRUE)

  # Check whether each sample_id + part maps to a unique label
  chk <- out[, .(n_label = uniqueN(label), labels = paste(unique(label), collapse = " | ")),
             by = .(sample_id, part)]
  bad_chk <- chk[n_label != 1]
  if (nrow(bad_chk) > 0) {
    stop(
      "Some sample_id + part combinations do not map to a unique label:\n",
      paste(
        apply(bad_chk, 1, function(x) {
          paste0("sample_id=", x[["sample_id"]], ", part=", x[["part"]], ", labels=", x[["labels"]])
        }),
        collapse = "\n"
      )
    )
  }

  out
}

find_meta_tsvs <- function(workdir, map_sub) {
  stem1 <- file.path(workdir, "bin_metadata_with_cluster_raw_stem1.tsv")
  stem2 <- file.path(workdir, "bin_metadata_with_cluster_raw_stem2.tsv")
  mainf <- file.path(workdir, "bin_metadata_with_cluster_raw.tsv")

  expected_parts <- sort(unique(na.omit(as.character(map_sub$part))))
  file_parts <- character(0)

  if (file.exists(stem1)) file_parts <- c(file_parts, "stem1")
  if (file.exists(stem2)) file_parts <- c(file_parts, "stem2")
  if (length(file_parts) == 0 && file.exists(mainf)) file_parts <- c(file_parts, "allbins")

  message("[find_meta_tsvs] workdir: ", workdir)
  message("[find_meta_tsvs] expected parts from Excel: ", paste(expected_parts, collapse = ", "))
  message("[find_meta_tsvs] detected parts from files: ", paste(file_parts, collapse = ", "))

  if (!setequal(expected_parts, file_parts)) {
    stop(
      "Mismatch between Excel parts and metadata files in: ", workdir,
      " | expected parts = ", paste(expected_parts, collapse = ", "),
      " | detected parts = ", paste(file_parts, collapse = ", ")
    )
  }

  out <- character(0)
  if ("allbins" %in% expected_parts) out <- c(out, mainf)
  if ("stem1"  %in% expected_parts) out <- c(out, stem1)
  if ("stem2"  %in% expected_parts) out <- c(out, stem2)

  out[file.exists(out)]
}

get_part_name <- function(meta_tsv) {
  part_name <- basename(meta_tsv)
  part_name <- sub("^bin_metadata_with_cluster_raw_", "", part_name)
  part_name <- sub("^bin_metadata_with_cluster_raw", "allbins", part_name)
  part_name <- sub("\\.tsv$", "", part_name)
  part_name
}

get_label_for_part <- function(map_sub, part_name, sample_id = NA_character_) {
  sub_dt <- unique(
    map_sub[part == part_name, .(part, label)]
  )

  sub_dt <- sub_dt[!is.na(label) & label != ""]

  message("[get_label_for_part] sample_id: ", sample_id)
  message("[get_label_for_part] part_name: ", part_name)
  print(sub_dt)

  if (nrow(sub_dt) != 1) {
    stop(
      "Cannot uniquely determine label for sample_id = ", sample_id,
      ", part = ", part_name,
      ". Candidate labels: ", paste(unique(sub_dt$label), collapse = ", ")
    )
  }

  sub_dt$label[[1]]
}

#make_sample_ann_colors <- function(sample_levels) {
#  sample_levels <- unique(sample_levels)
#  n <- length(sample_levels)
#  pal <- if (n <= 8) hcl.colors(n, "Set 2") else hcl.colors(n, "Dynamic")
#  setNames(pal, sample_levels)
#}

make_sample_ann_colors <- function(sample_levels) {
  sample_levels <- unique(sample_levels)

  cols <- sapply(sample_levels, function(x) {
    if (grepl("BIO1", x, ignore.case = TRUE)) {
      "black"
    } else if (grepl("BIO2", x, ignore.case = TRUE)) {
      "#9B9B9B"
    } else {
      "black"  # fallback
    }
  })

  setNames(cols, sample_levels)
}

# ============================================================
# 1) Load mapping table
# Required columns in raw Excel:
# sample_id, um, k, cluster, domain, label
# Optional column:
# part
# ============================================================
map_dt <- as.data.table(read_excel(MAP_XLSX))

if ("Labels" %in% colnames(map_dt) && !("label" %in% colnames(map_dt))) {
  setnames(map_dt, "Labels", "label")
}
if ("Label" %in% colnames(map_dt) && !("label" %in% colnames(map_dt))) {
  setnames(map_dt, "Label", "label")
}
if ("Part" %in% colnames(map_dt) && !("part" %in% colnames(map_dt))) {
  setnames(map_dt, "Part", "part")
}
if ("Domain" %in% colnames(map_dt) && !("domain" %in% colnames(map_dt))) {
  setnames(map_dt, "Domain", "domain")
}
if ("Cluster" %in% colnames(map_dt) && !("cluster" %in% colnames(map_dt))) {
  setnames(map_dt, "Cluster", "cluster")
}
if ("Sample_ID" %in% colnames(map_dt) && !("sample_id" %in% colnames(map_dt))) {
  setnames(map_dt, "Sample_ID", "sample_id")
}

required_cols <- c("sample_id", "um", "k", "cluster", "domain", "label")
if (!all(required_cols %in% colnames(map_dt))) {
  stop(
    "Missing required columns in Excel. Required = ",
    paste(required_cols, collapse = ", "),
    " | Found = ",
    paste(colnames(map_dt), collapse = ", ")
  )
}

map_dt <- expand_mapping_table(map_dt)

map_dt <- map_dt[um == UM_USE & k == K_USE]

message("Samples in mapping table:")
print(unique(map_dt$sample_id))

message("Expanded mapping table:")
print(
  map_dt[, .(
    n_rows = .N,
    n_label_unique = uniqueN(label),
    labels = paste(unique(label), collapse = " | ")
  ), by = .(sample_id, part)][order(sample_id, part)]
)

message("Check sample/part/domain mapping:")
print(
  map_dt[, .(
    n_rows = .N,
    clusters = paste(sort(unique(cluster)), collapse = ","),
    domains  = paste(unique(domain), collapse = " | "),
    label    = paste(unique(label), collapse = " | ")
  ), by = .(sample_id, part)][order(sample_id, part)]
)

# ============================================================
# 2) Function to run one sample
# ============================================================
run_deg_one_sample <- function(sample_id, map_sub, base_dir, out_dir) {
  workdir <- file.path(base_dir, sample_id, UM_USE, "kmeans_k10_raw_out")
  counts_rds <- file.path(workdir, "../counts_bins_by_genes_sparse.rds")
  meta_tsvs  <- find_meta_tsvs(workdir, map_sub)

  message("======================================")
  message("Processing sample: ", sample_id)
  message("workdir: ", workdir)

  if (!file.exists(counts_rds)) stop("Missing counts_rds: ", counts_rds)

  message("[", sample_id, "] Metadata files found:")
  print(meta_tsvs)

  mat0 <- readRDS(counts_rds)
  mat0 <- as(mat0, "dgCMatrix")

  sample_out <- file.path(out_dir, sample_id)
  dir.create(sample_out, showWarnings = FALSE, recursive = TRUE)

  deg_all_parts <- list()

  for (meta_tsv in meta_tsvs) {
    part_name  <- get_part_name(meta_tsv)
    part_label <- get_label_for_part(map_sub, part_name, sample_id)

    message("--------------------------------------")
    message("Processing metadata: ", meta_tsv)
    message("Part name: ", part_name)
    message("Part label: ", part_label)

    meta <- fread(meta_tsv)
    stopifnot(all(c("bin_id", "cluster_raw") %in% colnames(meta)))
    meta[, bin_id := as.character(bin_id)]
    meta[, cluster_raw := as.integer(cluster_raw)]

    mat <- mat0

    n_match_cols <- sum(colnames(mat) %in% meta$bin_id)
    n_match_rows <- sum(rownames(mat) %in% meta$bin_id)

    message("[", part_label, "] Match to metadata: colnames=", n_match_cols,
            " rownames=", n_match_rows)

    if (n_match_cols == 0 && n_match_rows > 0) {
      message("[", part_label, "] Detected bins in rownames(mat). Transposing matrix ...")
      mat <- t(mat)
    }

    n_match_cols2 <- sum(colnames(mat) %in% meta$bin_id)
    if (n_match_cols2 == 0) {
      stop("[", part_label, "] Cannot match any bin_id to matrix colnames.")
    }

    sample_map <- copy(map_sub[part == part_name])
    setnames(sample_map, c("cluster", "domain"), c("cluster_raw", "cluster5"))

    meta <- merge(
      meta,
      unique(sample_map[, .(cluster_raw, cluster5)]),
      by = "cluster_raw",
      all.x = TRUE
    )

    if (any(is.na(meta$cluster5))) {
      bad <- sort(unique(meta[is.na(cluster5), cluster_raw]))
      stop("[", part_label, "] Unmapped cluster_raw found: ", paste(bad, collapse = ", "))
    }

    meta <- meta[cluster5 != "background"]
    meta[, cluster5 := factor(cluster5, levels = DOMAIN_ORDER)]

    message("[", part_label, "] Bin counts per domain:")
    print(table(meta$cluster5, useNA = "ifany"))

    bins_in_mat <- colnames(mat)
    idx <- match(bins_in_mat, meta$bin_id)

    keep <- !is.na(idx)
    mat <- mat[, keep, drop = FALSE]
    bins_in_mat <- colnames(mat)
    idx <- match(bins_in_mat, meta$bin_id)

    if (anyNA(idx)) {
      stop("[", part_label, "] Still found unmatched bins after filtering.")
    }

    meta_aligned <- meta[idx, ]
    stopifnot(all(meta_aligned$bin_id == bins_in_mat))

    seu <- CreateSeuratObject(
      counts = mat,
      assay = "Xenium",
      meta.data = as.data.frame(meta_aligned)
    )

    DefaultAssay(seu) <- "Xenium"
    Idents(seu) <- "cluster5"

    seu <- JoinLayers(seu)
    seu <- NormalizeData(seu)
    seu <- JoinLayers(seu)

    part_out <- file.path(sample_out, part_name)
    dir.create(part_out, showWarnings = FALSE, recursive = TRUE)

    deg_list <- list()
    for (cl in levels(Idents(seu))) {
      if (sum(Idents(seu) == cl) == 0) next

      message("[", part_label, "] Running DEG: ", cl, " vs rest")

      deg <- FindMarkers(
        seu,
        ident.1 = cl,
        ident.2 = NULL,
        only.pos = FALSE,
        assay = "Xenium",
        layer = "data"
      )

      deg <- as.data.table(deg, keep.rownames = "gene")
      deg[, sample := sample_id]
      deg[, sample_part := part_name]
      deg[, display_label := part_label]
      deg[, cluster := cl]

      setcolorder(
        deg,
        c("sample", "sample_part", "display_label", "cluster", "gene",
          setdiff(names(deg), c("sample", "sample_part", "display_label", "cluster", "gene")))
      )

      #fwrite(deg, file.path(part_out, paste0("DEG_", cl, "_vs_rest.tsv")), sep = "\t")
      deg_list[[cl]] <- deg
    }

    deg_all <- rbindlist(deg_list, use.names = TRUE, fill = TRUE)
    deg_all_sig <- deg_all[abs(avg_log2FC) >= 1 & p_val_adj < 0.05]

    #fwrite(deg_all, file.path(part_out, "DEG_all_5class_one_vs_rest.tsv"), sep = "\t")
    #fwrite(deg_all_sig, file.path(part_out, "DEG_all_5class_one_vs_rest_sig.tsv"), sep = "\t")

    deg_all_parts[[part_label]] <- deg_all
  }

  rbindlist(deg_all_parts, use.names = TRUE, fill = TRUE)
}

# ============================================================
# 3) Run all samples
# ============================================================
sample_ids <- unique(map_dt$sample_id)

deg_all_list <- lapply(sample_ids, function(sid) {
  map_sub <- map_dt[sample_id == sid]
  run_deg_one_sample(
    sample_id = sid,
    map_sub   = map_sub,
    base_dir  = BASE_DIR,
    out_dir   = OUT_DIR
  )
})

deg_merged <- rbindlist(deg_all_list, use.names = TRUE, fill = TRUE)
#fwrite(deg_merged, file.path(OUT_DIR, "DEG_merged_all_samples.tsv"), sep = "\t")

# ============================================================
# 4) Build merged heatmap matrix
# ============================================================
deg_merged[, cluster := as.character(cluster)]
deg_merged[, gene := as.character(gene)]
deg_merged[, sample := as.character(sample)]
deg_merged[, sample_part := as.character(sample_part)]
deg_merged[, display_label := as.character(display_label)]

deg_merged[, abs_lfc := abs(avg_log2FC)]
setorder(deg_merged, display_label, gene, cluster, -abs_lfc)
deg_unique <- deg_merged[, .SD[1], by = .(display_label, gene, cluster)]

deg_mean <- deg_unique[, .(
  avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
  n_display_label = .N
), by = .(gene, cluster)]

mat_dt <- dcast(deg_mean, gene ~ cluster, value.var = "avg_log2FC")

missing_domains <- setdiff(DOMAIN_ORDER, colnames(mat_dt))
if (length(missing_domains) > 0) {
  stop("Missing domains in merged DEG table: ", paste(missing_domains, collapse = ", "))
}

mat_dt <- mat_dt[, c("gene", DOMAIN_ORDER), with = FALSE]

mat <- as.matrix(mat_dt[, -1, with = FALSE])
rownames(mat) <- mat_dt$gene
mat[is.na(mat)] <- 0
mat[mat > 3] <- 3
mat[mat < -3] <- -3

gene_rank <- order(apply(abs(mat), 1, max, na.rm = TRUE), decreasing = TRUE)
mat <- mat[gene_rank, , drop = FALSE]

mat_plot <- if (!is.null(TOP_N) && TOP_N < nrow(mat)) {
  mat[seq_len(TOP_N), , drop = FALSE]
} else {
  mat
}

new_order <- c("epidermis", "phloem", "cambium", "parenchyma", "sclerenchyma")
mat_plot <- mat_plot[, new_order, drop = FALSE]

png(file.path(OUT_DIR, paste0("Merged_heatmap_top", nrow(mat_plot), ".png")),
    width = 1400, height = 2200, res = 200)
pheatmap(
  mat_plot,
  color = heat_colors,
  breaks = heat_breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,  
  #fontsize_row = 6,
  fontsize_col = 12,
  main = sprintf("", nrow(mat_plot))
)
dev.off()

pdf(file.path(OUT_DIR, paste0("Merged_heatmap_top", nrow(mat_plot), ".pdf")),
    width = 8.5, height = 12)
pheatmap(
  mat_plot,
  color = heat_colors,
  breaks = heat_breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  #fontsize_row = 6,
  fontsize_col = 12,
  main = sprintf("", nrow(mat_plot))
)
dev.off()

#fwrite(
#  data.table(gene = rownames(mat_plot), as.data.table(mat_plot)),
#  file.path(OUT_DIR, paste0("Merged_heatmap_matrix_top", nrow(mat_plot), ".tsv")),
#  sep = "\t"
#)

message("Done average heatmap.")
message("Merged DEG: ", file.path(OUT_DIR, "DEG_merged_all_samples.tsv"))

# ============================================================
# 5) Build sample-separated heatmap matrix
# ============================================================
deg_sep <- copy(deg_merged)

deg_sep[, sample := as.character(sample)]
deg_sep[, sample_part := as.character(sample_part)]
deg_sep[, cluster := as.character(cluster)]
deg_sep[, gene := as.character(gene)]
deg_sep[, display_label := as.character(display_label)]

deg_sep[, abs_lfc := abs(avg_log2FC)]
setorder(deg_sep, display_label, gene, cluster, -abs_lfc)
deg_sep_unique <- deg_sep[, .SD[1], by = .(display_label, gene, cluster)]

deg_sep_unique[, sample_domain := paste(display_label, cluster, sep = "__")]

mat_sep_dt <- dcast(deg_sep_unique, gene ~ sample_domain, value.var = "avg_log2FC")

mat_sep <- as.matrix(mat_sep_dt[, -1, with = FALSE])
rownames(mat_sep) <- mat_sep_dt$gene
mat_sep[is.na(mat_sep)] <- 0
mat_sep[mat_sep > 3] <- 3
mat_sep[mat_sep < -3] <- -3

sample_display_order <- sort(unique(deg_sep_unique$display_label))

col_order <- unlist(lapply(sample_display_order, function(sid) {
  paste(sid, DOMAIN_ORDER, sep = "__")
}))
col_order <- intersect(col_order, colnames(mat_sep))
mat_sep <- mat_sep[, col_order, drop = FALSE]

gene_rank_sep <- order(apply(abs(mat_sep), 1, max, na.rm = TRUE), decreasing = TRUE)
mat_sep <- mat_sep[gene_rank_sep, , drop = FALSE]

mat_sep_plot <- if (!is.null(TOP_N) && TOP_N < nrow(mat_sep)) {
  mat_sep[seq_len(TOP_N), , drop = FALSE]
} else {
  mat_sep
}

col_anno <- data.frame(
  #sample = sub("__.*$", "", colnames(mat_sep_plot)),
  domain = sub("^.*__", "", colnames(mat_sep_plot)),
  stringsAsFactors = FALSE
)
rownames(col_anno) <- colnames(mat_sep_plot)

#col_anno$sample <- factor(col_anno$sample, levels = sample_display_order)
col_anno$domain <- factor(col_anno$domain, levels = DOMAIN_ORDER)

#sample_levels <- levels(col_anno$sample)
#sample_levels <- sample_levels[!is.na(sample_levels)]

#ann_colors <- list(
  #domain = domain_colors,
  #sample = make_sample_ann_colors(sample_levels)
#)

ann_colors <- list(
  domain = domain_colors
)

#sample_counts <- table(col_anno$sample)
#gaps_col <- cumsum(as.integer(sample_counts))
#if (length(gaps_col) > 0) gaps_col <- gaps_col[-length(gaps_col)]

# ============================================================
# 5) Build sample-separated heatmap matrix
# Keep the previous mat_sep_plot preprocessing unchanged
# ============================================================

col_anno <- data.frame(
  # Uncomment this line to keep sample color annotations above the heatmap
  # sample = sub("__.*$", "", colnames(mat_sep_plot)),
  domain = sub("^.*__", "", colnames(mat_sep_plot)),
  stringsAsFactors = FALSE
)
rownames(col_anno) <- colnames(mat_sep_plot)

col_anno$domain <- factor(col_anno$domain, levels = DOMAIN_ORDER)

ann_colors <- list(
  domain = domain_colors
)


# ------------------------------------------------------------
# Correct the gaps_col calculation by directly comparing changes in column names
# ------------------------------------------------------------
# 1. Extract the sample name corresponding to each column in the plotted matrix
actual_samples <- sub("__.*$", "", colnames(mat_sep_plot))

# 2. Identify the index positions where adjacent columns belong to different samples
gaps_col <- which(actual_samples[-1] != actual_samples[-length(actual_samples)])
# ------------------------------------------------------------

png(file.path(OUT_DIR, paste0("Merged_heatmap_bySample_top", nrow(mat_sep_plot), ".png")),
    width = 2600, height = 2200, res = 220)
pheatmap(
  mat_sep_plot,
  color = heat_colors,
  breaks = heat_breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,  # pheatmap requires cluster_cols to be FALSE for gaps_col to take effect
  show_rownames = FALSE,
  annotation_col = col_anno,
  annotation_colors = ann_colors,
  gaps_col = gaps_col,   # Use the newly calculated sample boundaries
  #fontsize_row = 6,
  fontsize_col = 9,
  angle_col = 45,
  legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
  legend_labels = c("<= -3", "-2", "-1", "0", "1", "2", ">= 3"),  

  main = sprintf("", nrow(mat_sep_plot))
)
dev.off()

pdf(file.path(OUT_DIR, paste0("Merged_heatmap_bySample_top", nrow(mat_sep_plot), ".pdf")),
    width = 18, height = 12)
pheatmap(
  mat_sep_plot,
  color = heat_colors,
  breaks = heat_breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE, 
  show_rownames = FALSE,
  annotation_col = col_anno,
  annotation_colors = ann_colors,
  gaps_col = gaps_col,
  #fontsize_row = 6,
  fontsize_col = 9,
  angle_col = 45,
  legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
  legend_labels = c("<= -3", "-2", "-1", "0", "1", "2", ">= 3"),  
  main = sprintf("", nrow(mat_sep_plot))
)
dev.off()


#fwrite(
#  data.table(gene = rownames(mat_sep_plot), as.data.table(mat_sep_plot)),
#  file.path(OUT_DIR, paste0("Merged_heatmap_bySample_matrix_top", nrow(mat_sep_plot), ".tsv")),
#  sep = "\t"
#)

message("Sample-separated heatmap PNG: ",
        file.path(OUT_DIR, paste0("Merged_heatmap_bySample_top", nrow(mat_sep_plot), ".png")))
message("Sample-separated heatmap PDF: ",
        file.path(OUT_DIR, paste0("Merged_heatmap_bySample_top", nrow(mat_sep_plot), ".pdf")))


###### Cluster columns with sample and domain annotation

# ============================================================
# 6) Build sample-separated heatmap matrix
# ============================================================
deg_sep <- copy(deg_merged)

deg_sep[, sample := as.character(sample)]
deg_sep[, sample_part := as.character(sample_part)]
deg_sep[, cluster := as.character(cluster)]
deg_sep[, gene := as.character(gene)]
deg_sep[, display_label := as.character(display_label)]

deg_sep[, abs_lfc := abs(avg_log2FC)]
setorder(deg_sep, display_label, gene, cluster, -abs_lfc)
deg_sep_unique <- deg_sep[, .SD[1], by = .(display_label, gene, cluster)]

deg_sep_unique[, sample_domain := paste(display_label, cluster, sep = "__")]

mat_sep_dt <- dcast(deg_sep_unique, gene ~ sample_domain, value.var = "avg_log2FC")

mat_sep <- as.matrix(mat_sep_dt[, -1, with = FALSE])
rownames(mat_sep) <- mat_sep_dt$gene
mat_sep[is.na(mat_sep)] <- 0
mat_sep[mat_sep > 3] <- 3
mat_sep[mat_sep < -3] <- -3

gene_rank_sep <- order(apply(abs(mat_sep), 1, max, na.rm = TRUE), decreasing = TRUE)
mat_sep <- mat_sep[gene_rank_sep, , drop = FALSE]

mat_sep_plot <- if (!is.null(TOP_N) && TOP_N < nrow(mat_sep)) {
  mat_sep[seq_len(TOP_N), , drop = FALSE]
} else {
  mat_sep
}

col_anno <- data.frame(
  #sample = sub("__.*$", "", colnames(mat_sep_plot)),
  domain = sub("^.*__", "", colnames(mat_sep_plot)),
  stringsAsFactors = FALSE
)
rownames(col_anno) <- colnames(mat_sep_plot)

#col_anno$sample <- factor(col_anno$sample, levels = sort(unique(col_anno$sample)))
col_anno$domain <- factor(col_anno$domain, levels = DOMAIN_ORDER)

#sample_levels <- levels(col_anno$sample)
#sample_levels <- sample_levels[!is.na(sample_levels)]

#ann_colors <- list(
#  domain = domain_colors,
#  sample = make_sample_ann_colors(sample_levels)
#)

ann_colors <- list(
  domain = domain_colors
)
png(file.path(OUT_DIR, paste0("Merged_heatmap_bySample_top", nrow(mat_sep_plot), "_column.png")),
    width = 2600, height = 2200, res = 220)
pheatmap(
  mat_sep_plot,
  color = heat_colors,
  breaks = heat_breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = col_anno,
  annotation_colors = ann_colors,
  fontsize_col = 9,
  angle_col = 45,
  main = sprintf("", nrow(mat_sep_plot))
)
dev.off()

pdf(file.path(OUT_DIR, paste0("Merged_heatmap_bySample_top", nrow(mat_sep_plot), "_column.pdf")),
    width = 18, height = 12)
pheatmap(
  mat_sep_plot,
  color = heat_colors,
  breaks = heat_breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = col_anno,
  annotation_colors = ann_colors,
  fontsize_col = 9,
  angle_col = 45,
  main = sprintf("", nrow(mat_sep_plot))
)
dev.off()

#fwrite(
#  data.table(gene = rownames(mat_sep_plot), as.data.table(mat_sep_plot)),
#  file.path(OUT_DIR, paste0("Merged_heatmap_bySample_matrix_top", nrow(mat_sep_plot), "_column.tsv")),
#  sep = "\t"
#)

message("Sample-separated heatmap PNG: ",
        file.path(OUT_DIR, paste0("Merged_heatmap_bySample_top", nrow(mat_sep_plot), "_column.png")))
message("Sample-separated heatmap PDF: ",
        file.path(OUT_DIR, paste0("Merged_heatmap_bySample_top", nrow(mat_sep_plot), "_column.pdf")))


###### Cluster columns with domain annotation

col_anno <- data.frame(
  domain = sub("^.*__", "", colnames(mat_sep_plot)),
  stringsAsFactors = FALSE
)
rownames(col_anno) <- colnames(mat_sep_plot)

col_anno$domain <- factor(col_anno$domain, levels = DOMAIN_ORDER)

ann_colors <- list(
  domain = domain_colors
)

png(file.path(OUT_DIR, paste0("Merged_heatmap_bySample_top", nrow(mat_sep_plot), "_column_domain.png")),
    width = 2600, height = 2200, res = 220)
pheatmap(
  mat_sep_plot,
  color = heat_colors,
  breaks = heat_breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = col_anno,
  annotation_colors = ann_colors,
  fontsize_col = 9,
  angle_col = 45,
  main = sprintf("", nrow(mat_sep_plot))
)
dev.off()

pdf(file.path(OUT_DIR, paste0("Merged_heatmap_bySample_top", nrow(mat_sep_plot), "_column_domain.pdf")),
    width = 18, height = 12)
pheatmap(
  mat_sep_plot,
  color = heat_colors,
  breaks = heat_breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = col_anno,
  annotation_colors = ann_colors,
  fontsize_col = 9,
  angle_col = 45,
  main = sprintf("", nrow(mat_sep_plot))
)
dev.off()

#fwrite(
#  data.table(gene = rownames(mat_sep_plot), as.data.table(mat_sep_plot)),
#  file.path(OUT_DIR, paste0("Merged_heatmap_bySample_matrix_top", nrow(mat_sep_plot), "_column_domain.tsv")),
#  sep = "\t"
#)

message("Sample-separated heatmap PNG: ",
        file.path(OUT_DIR, paste0("Merged_heatmap_bySample_top", nrow(mat_sep_plot), "_column_domain.png")))
message("Sample-separated heatmap PDF: ",
        file.path(OUT_DIR, paste0("Merged_heatmap_bySample_top", nrow(mat_sep_plot), "_column_domain.pdf")))


