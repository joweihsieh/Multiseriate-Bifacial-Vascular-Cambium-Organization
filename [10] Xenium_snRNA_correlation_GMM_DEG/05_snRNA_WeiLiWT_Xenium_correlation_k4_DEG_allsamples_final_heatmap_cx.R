#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
})

setwd("/home/woodydrylab/FileShare/20260117_MultiLayer_Unifacial/two_cambium")

# ============================================================
# Input
# ============================================================
OUTDIR_BASE <- "GMM_on_mean_median_correlation/snRNA_DEG_from_mean_median_GMM_codex"

# mean_based 或 median_based
MODE_USE <- "mean_based"
# MODE_USE <- "median_based"

INDIR  <- file.path(OUTDIR_BASE, MODE_USE)
OUTDIR <- file.path(INDIR, "heatmap_log2FC")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

DEG_SIG_FILE <- file.path(INDIR, "DEG_Cphloem_Cxylem_one_vs_rest_sig.tsv")
DEG_OPP_FILE <- file.path(INDIR, "DEG_Cphloem_Cxylem_Opposite_Regulation.tsv")

FONTSIZE_ROW_ALL <- 8
FONTSIZE_ROW_OPP <- 7

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

pick_one_row_per_gene_comparison <- function(dt, fc_col = "avg_log2FC", padj_col = "p_val_adj") {
  dt <- copy(dt)

  if (!"gene_clean" %in% names(dt)) stop("Missing gene_clean")
  if (!"comparison" %in% names(dt)) stop("Missing comparison")
  if (!fc_col %in% names(dt)) stop("Missing FC column")
  if (!padj_col %in% names(dt)) stop("Missing padj column")

  dt <- dt[!is.na(gene_clean) & gene_clean != ""]
  dt[, gene_clean := clean_gene_id(gene_clean)]
  dt[, abs_fc := abs(get(fc_col))]
  setorderv(dt, c("comparison", "gene_clean", padj_col, "abs_fc"), c(1, 1, 1, -1), na.last = TRUE)
  dt <- dt[, .SD[1], by = .(comparison, gene_clean)]
  dt[, abs_fc := NULL]
  dt
}

make_gene_labels <- function(dt, gene_order) {
  labels <- as.character(gene_order)
  labels <- sub("^POTRI_", "Potri.", labels)
  make.unique(labels)
}

plot_logfc_heatmap <- function(deg_file,
                               prefix,
                               opposite_file = NULL,
                               show_row_names = FALSE,
                               fontsize_row = 8,
                               png_height = 1200,
                               png_width = 800) {
  if (!file.exists(deg_file)) {
    warning("File not found: ", deg_file)
    return(NULL)
  }

  dt <- fread(deg_file)

  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), names(dt))[1]
  if (length(fc_col) == 0 || is.na(fc_col)) {
    stop("Cannot find avg_log2FC / avg_logFC in: ", deg_file)
  }

  if (!"gene_clean" %in% names(dt)) {
    if ("gene" %in% names(dt)) {
      dt[, gene_clean := clean_gene_id(gene)]
    } else {
      stop("No gene_clean or gene column in: ", deg_file)
    }
  } else {
    dt[, gene_clean := clean_gene_id(gene_clean)]
  }

  if (!"comparison" %in% names(dt)) {
    stop("comparison column not found in: ", deg_file)
  }

  dt <- pick_one_row_per_gene_comparison(dt, fc_col = fc_col, padj_col = "p_val_adj")

  dt_wide <- dcast(
    dt,
    gene_clean ~ comparison,
    value.var = fc_col,
    fill = 0
  )

  gene_order <- dt_wide$gene_clean
  mat <- as.matrix(dt_wide[, setdiff(names(dt_wide), "gene_clean"), with = FALSE])
  row_labels <- make_gene_labels(dt, gene_order)
  rownames(mat) <- row_labels

  # 只有在真的提供 opposite_file 時才畫 annotation
  row_anno <- NULL
  anno_colors <- NULL
  regulation_vec <- rep(NA_character_, length(gene_order))

  if (!is.null(opposite_file) && file.exists(opposite_file)) {
    dt_opp <- fread(opposite_file)

    if (!"gene_clean" %in% names(dt_opp)) {
      if ("gene" %in% names(dt_opp)) {
        dt_opp[, gene_clean := clean_gene_id(gene)]
      } else {
        stop("No gene_clean or gene column in opposite file")
      }
    } else {
      dt_opp[, gene_clean := clean_gene_id(gene_clean)]
    }

    opp_genes <- unique(dt_opp$gene_clean)

    regulation_vec <- ifelse(gene_order %in% opp_genes, "Opposite", "Other")
    row_anno <- data.frame(Regulation = regulation_vec)
    rownames(row_anno) <- row_labels

    anno_colors <- list(
      Regulation = c(
        Opposite = "black",
        Other = "grey80"
      )
    )
  }

  fwrite(
    data.table(
      gene_clean = gene_order,
      label = row_labels,
      Regulation = regulation_vec
    ),
    file.path(OUTDIR, paste0(prefix, "_row_labels.tsv")),
    sep = "\t"
  )

  write.csv(
    mat,
    file = file.path(OUTDIR, paste0(prefix, "_matrix_log2FC.csv")),
    row.names = TRUE
  )

  png(
    filename = file.path(OUTDIR, paste0(prefix, "_heatmap.png")),
    width = png_width,
    height = png_height,
    res = 200
  )

  max_abs <- max(abs(mat), na.rm = TRUE)
  my_breaks <- seq(-max_abs, max_abs, length.out = 201)
  my_colors <- colorRampPalette(c("#005493", "white", "#941100"))(200)

  pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = show_row_names,
    fontsize_row = fontsize_row,
    border_color = NA,
    annotation_row = row_anno,
    annotation_colors = anno_colors,
    color = my_colors,
    breaks = my_breaks
  )

  dev.off()

  invisible(mat)
}

# ============================================================
# Run
# ============================================================
message("[1/2] Plotting one-vs-rest significant log2FC heatmap")
plot_logfc_heatmap(
  deg_file = DEG_SIG_FILE,
  prefix = paste0(gsub("_based$", "", MODE_USE), "_one_vs_rest_sig_ALL_log2FC"),
  opposite_file = DEG_OPP_FILE,
  show_row_names = FALSE,
  fontsize_row = FONTSIZE_ROW_ALL,
  png_height = 1200,
  png_width = 800
)

if (file.exists(DEG_OPP_FILE)) {
  message("[2/2] Plotting opposite-regulation log2FC heatmap")
  plot_logfc_heatmap(
    deg_file = DEG_OPP_FILE,
    prefix = paste0(gsub("_based$", "", MODE_USE), "_opposite_regulation_log2FC"),
    opposite_file = NULL,
    show_row_names = TRUE,
    fontsize_row = FONTSIZE_ROW_OPP,
    png_height = 900,
    png_width = 700
  )
} else {
  message("[2/2] Opposite regulation file not found, skip.")
}

message("Done: ", normalizePath(OUTDIR))