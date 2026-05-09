#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

setwd("/home/woodydrylab/FileShare/20260121_Xenium")

# =======================
# INPUTS (edit)
# =======================
LCM_CAMBIUM <- "/home/chingweilu/Multi_Bilayer_Unifacial/Adjusting_UMAP_neighbor/WeiL_LCM/quantification/Ptr_Cambium_gene_abundances.csv"
LCM_XYLEM   <- "/home/chingweilu/Multi_Bilayer_Unifacial/Adjusting_UMAP_neighbor/WeiL_LCM/quantification/Ptr_Xylem_gene_abundances.csv"
LCM_PHLOEM  <- "/home/chingweilu/Multi_Bilayer_Unifacial/Adjusting_UMAP_neighbor/WeiL_LCM/quantification/Ptr_Phloem_gene_abundances.csv"

# Xenium 300 gene list (one gene per line) OR read from your matrix colnames
XENIUM_300_TXT <- "/home/woodydrylab/FileShare/20260121_Xenium/output-XETG00360__0079031__BIO2_TISSUE_1__20260115__224443/grid05um_out/genes.tsv"

# output
OUTDIR <- "LCM_probe_representativeness_300_probes"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# =======================
# Helpers
# =======================
clean_gene_id <- function(x) {
  x <- sub("\\.v4\\.1$", "", x)
  x <- gsub("^Potri\\.", "Potri_", x)
  x
}

# entropy (0..log(3)); higher = more even across tissues
entropy3 <- function(p, c, x, eps=1e-12) {
  v <- pmax(c(p, c, x), eps)
  v <- v / sum(v)
  -sum(v * log(v))
}

# ternary -> 2D simplex coordinates
# vertices: Phloem=(0,0), Xylem=(1,0), Cambium=(0.5, sqrt(3)/2)
ternary_to_xy <- function(p, c, x) {
  s <- p + c + x
  p <- p / s; c <- c / s; x <- x / s
  X <- x + 0.5 * c
  Y <- (sqrt(3) / 2) * c
  list(X = X, Y = Y, p = p, c = c, x = x)
}

# =======================
# Load LCM and compute tissue means
# =======================
dt_c <- fread(LCM_CAMBIUM)
dt_x <- fread(LCM_XYLEM)
dt_p <- fread(LCM_PHLOEM)

# --- You need to adapt these two lines to your real column names
# Assumption: there is a gene id column named "Gene.ID" (or similar),
# and remaining columns are replicates / samples for that tissue.
gene_col <- "Gene.ID"   # <-- change if needed

dt_c[, gene_clean := clean_gene_id(get(gene_col))]
dt_x[, gene_clean := clean_gene_id(get(gene_col))]
dt_p[, gene_clean := clean_gene_id(get(gene_col))]

c_cols <- setdiff(names(dt_c), c(gene_col, "gene_clean"))
x_cols <- setdiff(names(dt_x), c(gene_col, "gene_clean"))
p_cols <- setdiff(names(dt_p), c(gene_col, "gene_clean"))

dt_c[, cambium_mean := rowMeans(.SD, na.rm=TRUE), .SDcols=c_cols]
dt_x[, xylem_mean   := rowMeans(.SD, na.rm=TRUE), .SDcols=x_cols]
dt_p[, phloem_mean  := rowMeans(.SD, na.rm=TRUE), .SDcols=p_cols]

# merge three tissues by gene
dt <- merge(dt_p[, .(gene_clean, phloem_mean)],
            dt_c[, .(gene_clean, cambium_mean)], by="gene_clean", all=TRUE)
dt <- merge(dt, dt_x[, .(gene_clean, xylem_mean)], by="gene_clean", all=TRUE)

# replace missing with 0
for (cc in c("phloem_mean","cambium_mean","xylem_mean")) dt[is.na(get(cc)), (cc) := 0]

# remove genes with all zeros
#dt <- dt[(phloem_mean + cambium_mean + xylem_mean) > 0]

# =======================
# Xenium 300 list
# =======================
xen300 <- read.csv(XENIUM_300_TXT)
xen300 <- xen300$gene
xen300 <- clean_gene_id(xen300)
xen300 <- xen300[nzchar(xen300)]
dt[, is_xenium300 := gene_clean %in% xen300]

# =======================
# Build ternary coords + QC metrics
# =======================
xy <- ternary_to_xy(dt$phloem_mean, dt$cambium_mean, dt$xylem_mean)
dt[, `:=`(p = xy$p, c = xy$c, x = xy$x, X = xy$X, Y = xy$Y)]

dt[, max_prop := pmax(p, c, x)]
dt[, ent := mapply(entropy3, p, c, x)]

# =======================
# Build rank-based ternary coords + QC metrics
# =======================

# ---- Make sure we have these columns:
stopifnot(all(c("gene_clean", "phloem_mean", "cambium_mean", "xylem_mean") %in% names(dt)))

# ---- Rank within each tissue (1 = highest expression)
# ties.method can be: "average", "min", "max", "first", "random", "dense"
TIES_METHOD <- "average"

dt[, rank_phloem  := frank(-phloem_mean,  ties.method = TIES_METHOD, na.last = "keep")]
dt[, rank_cambium := frank(-cambium_mean, ties.method = TIES_METHOD, na.last = "keep")]
dt[, rank_xylem   := frank(-xylem_mean,   ties.method = TIES_METHOD, na.last = "keep")]

# ---- Convert ranks to 0..1 scores (1 = top ranked, 0 = bottom ranked)
# Use N (number of genes) per tissue; here all share the same dt rows
N <- nrow(dt)
rank_to_score <- function(r, N) {
  # rank 1 -> 1; rank N -> 0
  (N - r) / (N - 1)
}

dt[, score_p := rank_to_score(rank_phloem,  N)]
dt[, score_c := rank_to_score(rank_cambium, N)]
dt[, score_x := rank_to_score(rank_xylem,   N)]

# ---- Optional: emphasize top ranks (nonlinear transform)
# GAMMA > 1 will emphasize top-ranked genes more; set to 1 to disable
GAMMA <- 1
dt[, `:=`(
  score_p = score_p^GAMMA,
  score_c = score_c^GAMMA,
  score_x = score_x^GAMMA
)]

# ---- Normalize to proportions for ternary
s <- dt$score_p + dt$score_c + dt$score_x
# avoid 0/0: if a gene is bottom in all tissues, s can be 0 (rare if N>1, but safe)
s[s == 0] <- 1e-12

dt[, p := score_p / s]
dt[, c := score_c / s]
dt[, x := score_x / s]

# ---- Ternary -> XY
xy <- ternary_to_xy(dt$p, dt$c, dt$x)
dt[, `:=`(X = xy$X, Y = xy$Y)]

# ---- QC metrics (still meaningful on proportions)
dt[, max_prop := pmax(p, c, x)]
dt[, ent := mapply(entropy3, p, c, x)]

# Also: rank-spread metric (optional): small = similar ranks across tissues
dt[, rank_sd := apply(cbind(rank_phloem, rank_cambium, rank_xylem), 1, sd)]

# =======================
# Plot 1: ternary (all genes vs 300 probes) using ranking
# =======================
png(file.path(OUTDIR, "LCM_ternary_rank_based_all_genes.png"), width=2200, height=1800, res=300)
#par(mar=c(1,1,2,1), bg="white")
par(mar=c(1,1,2,1), bg="transparent")

plot(NA, NA, xlim=c(-0.05, 1.05), ylim=c(-0.05, sqrt(3)/2 + 0.05),
     xlab="", ylab="", axes=FALSE, asp=1,
     main="Rank-based ternary (within-tissue expression ranking)")

VP <- c(0,0); VX <- c(1,0); VC <- c(0.5, sqrt(3)/2)
segments(VP[1],VP[2], VX[1],VX[2], col="grey30", lwd=2)
segments(VX[1],VX[2], VC[1],VC[2], col="grey30", lwd=2)
segments(VC[1],VC[2], VP[1],VP[2], col="grey30", lwd=2)

text(VP[1], VP[2]-0.03, "Phloem (rank)",  adj=c(0.5,1), cex=1.2)
text(VX[1], VX[2]-0.03, "Xylem (rank)",   adj=c(0.5,1), cex=1.2)
text(VC[1], VC[2]+0.03, "Cambium (rank)", adj=c(0.5,0), cex=1.2)

points(dt$X, dt$Y, pch=16, cex=0.35, col=rgb(0.2,0.2,0.2,0.15))

sel <- dt[is_xenium300 == TRUE]
points(sel$X, sel$Y, pch=16, cex=0.75, col=rgb(1,0,0,0.85))

legend("topright",
       legend=c("All other genes", "Xenium probes"),
       pch=16, pt.cex=c(0.9,1.1),
       col=c(rgb(0.2,0.2,0.2,0.35), rgb(1,0,0,0.85)),
       bty="n")
dev.off()

