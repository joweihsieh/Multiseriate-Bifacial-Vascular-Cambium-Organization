# 1. Use 34,699 genes from the LCM dataset.
# 2. Within each tissue (phloem / cambium / xylem), rank genes by mean expression.
# 3. Divide all genes into three expression tertiles: Low / Mid / High.
# 4. Map the Xenium 300-probe genes back to these tertile groups to examine their distribution.

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
# Total number of genes: 34,699
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


######################################################################


library(data.table)

# dt: gene_clean, phloem_mean, cambium_mean, xylem_mean, is_xenium300

# --- helper: assign Low/Mid/High by tertiles within each tissue
assign_tertile <- function(v) {
  qs <- quantile(v, probs=c(1/3, 2/3), na.rm=TRUE)
  cut(v, breaks=c(-Inf, qs[1], qs[2], Inf),
      labels=c("Low","Mid","High"), include.lowest=TRUE, right=TRUE)
}

dt[, phloem_grp  := assign_tertile(phloem_mean)]
dt[, cambium_grp := assign_tertile(cambium_mean)]
dt[, xylem_grp   := assign_tertile(xylem_mean)]

# --- count Xenium probes in each group for each tissue
tab_ph <- dt[is_xenium300==TRUE, .N, by=phloem_grp][order(phloem_grp)]
tab_ca <- dt[is_xenium300==TRUE, .N, by=cambium_grp][order(cambium_grp)]
tab_xy <- dt[is_xenium300==TRUE, .N, by=xylem_grp][order(xylem_grp)]

# make a long table for plotting stacked bars
long <- rbind(
  data.table(tissue="Phloem",  group=as.character(tab_ph$phloem_grp),  n=tab_ph$N),
  data.table(tissue="Cambium", group=as.character(tab_ca$cambium_grp), n=tab_ca$N),
  data.table(tissue="Xylem",   group=as.character(tab_xy$xylem_grp),   n=tab_xy$N)
)

# ensure Low/Mid/High order
long[, group := factor(group, levels=c("Low","Mid","High"))]
long[, tissue := factor(tissue, levels=c("Phloem","Cambium","Xylem"))]

# also compute proportion (optional)
long[, prop := n / sum(n), by=tissue]

write.table(
  long,
  file.path(OUTDIR, "Xenium300_LowMidHigh_by_Tissue_value.csv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

#a <- dt[(phloem_mean + cambium_mean + xylem_mean) == 0]
#b <- a[a$gene_clean %in%xen300,]
#dt[dt$gene_clean%in%b$gene_clean, ]
# --- plot stacked bar (base R)
png(file.path(OUTDIR, "Xenium300_LowMidHigh_byTissue.png"),
    width = 1800, height = 1400, res = 250)

par(mar=c(6,5,3,2))

# build matrix for barplot
mat_counts <- dcast(long, group ~ tissue, value.var="n", fill=0)
m <- as.matrix(mat_counts[, -1])
rownames(m) <- mat_counts$group

barplot(m, beside=FALSE, las=1,
        ylab="Number of Xenium probes (n)",
        main="Xenium 300 probes span low/mid/high expression within each LCM tissue")
legend("topright", legend=rownames(m), bty="n")
dev.off()



######################################################################
library(data.table)

# ============================================================
# Helper: assign Low / Mid / High by tissue-specific tertiles
# ============================================================
assign_LMH <- function(v) {
  qs <- quantile(v, probs = c(1/3, 2/3), na.rm = TRUE)
  cut(
    v,
    breaks = c(-Inf, qs[1], qs[2], Inf),
    labels = c("Low", "Mid", "High"),
    include.lowest = TRUE
  )
}

# ============================================================
# Assign groups within each tissue
# ============================================================
dt[, cambium_grp := assign_LMH(cambium_mean)]
dt[, xylem_grp   := assign_LMH(xylem_mean)]
dt[, phloem_grp  := assign_LMH(phloem_mean)]

# ============================================================
# Count Xenium probes in each group per tissue
# ============================================================
tab_cambium <- dt[is_xenium300 == TRUE, .N, by = cambium_grp]
tab_xylem   <- dt[is_xenium300 == TRUE, .N, by = xylem_grp]
tab_phloem  <- dt[is_xenium300 == TRUE, .N, by = phloem_grp]

# unify into long format
long <- rbind(
  data.table(tissue = "Cambium", group = tab_cambium$cambium_grp, n = tab_cambium$N),
  data.table(tissue = "Xylem",   group = tab_xylem$xylem_grp,     n = tab_xylem$N),
  data.table(tissue = "Phloem",  group = tab_phloem$phloem_grp,   n = tab_phloem$N)
)

# enforce order
long[, group  := factor(group, levels = c("Low", "Mid", "High"))]
long[, tissue := factor(tissue, levels = c("Phloem", "Cambium", "Xylem"))]

# also compute proportions (optional, but good for text)
long[, prop := n / sum(n), by = tissue]

print(long)

# ============================================================
# Stacked bar plot (base R)
# ============================================================

# colors: Low / Mid / High
grp_cols <- c(
  Low  = "#A6A6A6",  # grey
  Mid  = "#FDB863",  # orange
  High = "#D73027"   # red
)

# convert to matrix for barplot
mat <- dcast(long, group ~ tissue, value.var = "n", fill = 0)
m <- as.matrix(mat[, -1])
rownames(m) <- mat$group


png(file.path(OUTDIR, "Xenium300_LowMidHigh_by_Tissue.png"),
    width = 1800, height = 1400, res = 250)

par(mar = c(6, 5, 3, 2))

barplot(
  m,
  beside = FALSE,
  col = grp_cols[rownames(m)],
  border = NA,
  las = 1,
  ylab = "Number of Xenium probes",
  main = "Xenium probes in each LCM tissue"
)

legend(
  "bottomright",
  legend = rownames(m),
  fill = grp_cols[rownames(m)],
  bty = "n"
)
dev.off()
