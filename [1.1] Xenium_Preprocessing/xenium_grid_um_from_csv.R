#!/usr/bin/env Rscript
# xenium_grid_um_from_csv.R

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(optparse)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="CSV with columns: x,y,gene", metavar="FILE"),
  make_option(c("-o","--outdir"), type="character", default="grid5um_out"),
  make_option(c("--bin_um"), type="double", default=5),
  make_option(c("--x0"), type="double", default=NA),
  make_option(c("--y0"), type="double", default=NA),
  make_option(c("--min_counts"), type="integer", default=0),
  make_option(c("--neighbor_k"), type="integer", default=0)
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || opt$input=="") stop("Please provide --input transcripts_xyz.csv")
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
dt <- fread(opt$input)

stopifnot(all(c("x","y","gene") %in% names(dt)))
dt <- dt[!is.na(x) & !is.na(y) & !is.na(gene) & gene!=""]

bin_um <- opt$bin_um
x0 <- opt$x0; y0 <- opt$y0
if (is.na(x0)) x0 <- min(dt$x)
if (is.na(y0)) y0 <- min(dt$y)

dt[, ix := as.integer(floor((x - x0)/bin_um))]
dt[, iy := as.integer(floor((y - y0)/bin_um))]

# counts per (bin,gene)
bg <- dt[, .(count=.N), by=.(ix,iy,gene)]

bins <- unique(bg[,.(ix,iy)])
setorder(bins, ix, iy); bins[, bin_row := .I]
setkey(bins, ix, iy); setkey(bg, ix, iy)
bg <- bins[bg]

genes <- sort(unique(bg$gene))
gmap <- data.table(gene=genes, gene_col=seq_along(genes))
setkey(gmap, gene); setkey(bg, gene)
bg <- gmap[bg]

mat <- sparseMatrix(
  i = bg$bin_row, j = bg$gene_col, x = bg$count,
  dims = c(nrow(bins), length(genes)),
  dimnames = list(paste0("bin_", bins$ix, "_", bins$iy), genes)
)

bin_meta <- data.table(
  bin_id = rownames(mat),
  ix = bins$ix, iy = bins$iy,
  x_center = x0 + (bins$ix + 0.5) * bin_um,
  y_center = y0 + (bins$iy + 0.5) * bin_um,
  total_counts = as.integer(Matrix::rowSums(mat)),
  n_genes = as.integer(Matrix::rowSums(mat>0))
)

# optional filtering
if (opt$min_counts > 0) {
  pass_counts <- bin_meta$total_counts >= opt$min_counts
  keep <- pass_counts

  if (opt$neighbor_k > 0) {
    if (opt$neighbor_k > 8) stop("--neighbor_k cannot exceed 8")
    idx_map <- bins[,.(ix,iy,idx=.I)]
    setkey(idx_map, ix, iy)
    offs <- data.table::CJ(dx=-1:1, dy=-1:1)[!(dx==0 & dy==0)]
    neigh_n <- integer(nrow(bins))
    for (k in seq_len(nrow(bins))) {
      neigh <- offs[,.(ix=bins$ix[k]+dx, iy=bins$iy[k]+dy)]
      nidx <- idx_map[neigh, on=.(ix,iy), idx]
      nidx <- nidx[!is.na(nidx)]
      neigh_n[k] <- sum(pass_counts[nidx])
    }
    keep <- pass_counts & (neigh_n >= opt$neighbor_k)
  }

  mat <- mat[keep,,drop=FALSE]
  bin_meta <- bin_meta[keep]
}

saveRDS(mat, file=file.path(opt$outdir, "counts_bins_by_genes_sparse.rds"))
fwrite(bin_meta, file=file.path(opt$outdir, "bin_metadata.tsv"), sep="\t")
fwrite(data.table(gene=colnames(mat)), file=file.path(opt$outdir, "genes.tsv"), sep="\t")

cat("Done. Outputs in:", opt$outdir, "\n")
