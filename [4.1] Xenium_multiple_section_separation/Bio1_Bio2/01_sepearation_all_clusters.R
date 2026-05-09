setwd("/home/woodydrylab/FileShare/20260121_Xenium/output-XETG00360__0079031__BIO1_TISSUE_1_and_2__20260115__224443/grid05um_out/kmeans_k10_raw_out/")
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# Inputs (edit)
# -----------------------------
BIN_META <- "bin_metadata_with_cluster_raw.tsv"  # needs: bin_id, x_center, y_center
STEM1_BD <- "/home/woodydrylab/FileShare/20260121_Xenium/output-XETG00360__0079031__BIO1_TISSUE_1_and_2__20260115__224443/Selection_Upper_coordinates.csv"    # your polygon 1
STEM2_BD <- "/home/woodydrylab/FileShare/20260121_Xenium/output-XETG00360__0079031__BIO1_TISSUE_1_and_2__20260115__224443/Selection_Bottom_coordinates.csv"   # your polygon 2
OUT_TSV  <- "bin_assignment_by_stem.tsv"

# -----------------------------
# Helper: robustly extract x,y from a data.table
# -----------------------------
extract_xy <- function(dt) {
  # If there are explicit x/y-like columns, use them; otherwise use first two numeric columns.
  nms <- names(dt)
  nms_l <- tolower(nms)

  pick_by_name <- function(cands) {
    hit <- which(nms_l %in% cands)
    if (length(hit) > 0) return(hit[1])
    # try startswith x_/y_
    hit <- which(grepl(paste0("^", cands[1], "($|[_\\.])"), nms_l))
    if (length(hit) > 0) return(hit[1])
    return(NA_integer_)
  }

  ix <- pick_by_name(c("x","x_center","xcoord","x_coordinate","xcoordinate","x_um","xcenter"))
  iy <- pick_by_name(c("y","y_center","ycoord","y_coordinate","ycoordinate","y_um","ycenter"))

  if (!is.na(ix) && !is.na(iy)) {
    x <- suppressWarnings(as.double(dt[[ix]]))
    y <- suppressWarnings(as.double(dt[[iy]]))
  } else {
    # fallback: first two numeric-ish columns
    cols_num <- which(vapply(dt, function(v) {
      v2 <- suppressWarnings(as.double(v))
      sum(is.finite(v2)) > 0
    }, logical(1)))

    if (length(cols_num) < 2) stop("Cannot infer x/y: boundary file has <2 numeric columns.")
    x <- suppressWarnings(as.double(dt[[cols_num[1]]]))
    y <- suppressWarnings(as.double(dt[[cols_num[2]]]))
  }

  out <- data.table(x = x, y = y)
  out <- out[is.finite(x) & is.finite(y)]
  if (nrow(out) < 3) stop(sprintf("Boundary has <3 valid points after cleaning (n=%d).", nrow(out)))
  out
}

close_ring <- function(xy) {
  # ensure first=last
  if (!(isTRUE(all.equal(xy$x[1], xy$x[.N])) && isTRUE(all.equal(xy$y[1], xy$y[.N])))) {
    xy <- rbind(xy, xy[1])
  }
  xy
}

# -----------------------------
# Pure R point-in-polygon (ray casting)
# Returns TRUE if inside or on edge
# -----------------------------
point_in_poly <- function(px, py, vx, vy) {
  # vx,vy: polygon vertices, closed (first=last)
  n <- length(vx) - 1L
  inside <- rep(FALSE, length(px))

  # Vectorized over points, loop over edges (fast enough)
  for (i in 1:n) {
    x1 <- vx[i];   y1 <- vy[i]
    x2 <- vx[i+1]; y2 <- vy[i+1]

    # Check if edge crosses horizontal ray to the right of point
    cond1 <- ( (y1 > py) != (y2 > py) )
    # Avoid division by zero; y2-y1 could be 0, but cond1 will be FALSE then.
    xinters <- (x2 - x1) * (py - y1) / (y2 - y1) + x1

    cond2 <- px <= xinters
    inside <- xor(inside, cond1 & cond2)
  }

  # Treat points on boundary as inside:
  # quick boundary check using cross-product & bounding box with tolerance
  tol <- 1e-9
  on_edge <- rep(FALSE, length(px))
  for (i in 1:n) {
    x1 <- vx[i];   y1 <- vy[i]
    x2 <- vx[i+1]; y2 <- vy[i+1]

    # bounding box filter
    bb <- (px >= pmin(x1,x2)-tol) & (px <= pmax(x1,x2)+tol) &
          (py >= pmin(y1,y2)-tol) & (py <= pmax(y1,y2)+tol)
    if (!any(bb)) next

    # cross product == 0 for colinear
    cross <- (px[bb]-x1)*(y2-y1) - (py[bb]-y1)*(x2-x1)
    on_edge[bb] <- on_edge[bb] | (abs(cross) <= tol)
  }

  inside | on_edge
}

centroid_xy <- function(xy) {
  # simple centroid (mean of vertices) is fine for tie-break
  c(mean(xy$x), mean(xy$y))
}

# -----------------------------
# Load bins
# -----------------------------
bins <- fread(BIN_META)
stopifnot(all(c("bin_id","x_center","y_center") %in% names(bins)))

bins[, `:=`(
  x_center = as.double(x_center),
  y_center = as.double(y_center)
)]
bins <- bins[is.finite(x_center) & is.finite(y_center)]

# -----------------------------
# Load boundaries (robustly)
# -----------------------------
bd1_raw <- fread(STEM1_BD, header = TRUE)
bd2_raw <- fread(STEM2_BD, header = TRUE)

bd1 <- close_ring(extract_xy(bd1_raw))
bd2 <- close_ring(extract_xy(bd2_raw))

# vertices
v1x <- bd1$x; v1y <- bd1$y
v2x <- bd2$x; v2y <- bd2$y

# -----------------------------
# Assign stem by polygon
# -----------------------------
in1 <- point_in_poly(bins$x_center, bins$y_center, v1x, v1y)
in2 <- point_in_poly(bins$x_center, bins$y_center, v2x, v2y)

bins[, stem := NA_character_]
bins[in1 & !in2, stem := "stem1"]
bins[in2 & !in1, stem := "stem2"]

# overlap -> nearest centroid
ov <- which(in1 & in2)
if (length(ov) > 0) {
  c1 <- centroid_xy(bd1)
  c2 <- centroid_xy(bd2)
  d1 <- (bins$x_center[ov]-c1[1])^2 + (bins$y_center[ov]-c1[2])^2
  d2 <- (bins$x_center[ov]-c2[1])^2 + (bins$y_center[ov]-c2[2])^2
  bins$stem[ov] <- ifelse(d1 <= d2, "stem1", "stem2")
}

# -----------------------------
# Save + summary
# -----------------------------
fwrite(bins[, .(bin_id, x_center, y_center, stem)], OUT_TSV, sep = "\t")
print(table(bins$stem, useNA = "ifany"))
cat("Wrote:", OUT_TSV, "\n")



suppressPackageStartupMessages({
  library(ggplot2)
})

# -----------------------------
# Prepare polygon dfs for ggplot
# -----------------------------
poly1_df <- data.table(x = bd1$x, y = bd1$y, poly = "stem1", order = seq_len(nrow(bd1)))
poly2_df <- data.table(x = bd2$x, y = bd2$y, poly = "stem2", order = seq_len(nrow(bd2)))

# -----------------------------
# Downsample points for plotting (optional but recommended)
# -----------------------------
set.seed(1)
MAX_PLOT_POINTS <- 200000  # adjust (50k~300k is usually fine)

plot_bins <- bins
if (nrow(plot_bins) > MAX_PLOT_POINTS) {
  plot_bins <- plot_bins[sample.int(nrow(plot_bins), MAX_PLOT_POINTS)]
}

# Make NA explicit for coloring
plot_bins[, stem_plot := fifelse(is.na(stem), "NA", stem)]

# -----------------------------
# Plot
# -----------------------------
p <- ggplot() +
  geom_path(data = poly1_df, aes(x = x, y = y), linewidth = 0.6) +
  geom_path(data = poly2_df, aes(x = x, y = y), linewidth = 0.6) +
  geom_point(
    data = plot_bins,
    aes(x = x_center, y = y_center, color = stem_plot),
    size = 0.15, alpha = 0.5
  ) +
  coord_equal() +
  scale_y_reverse()+
  theme_classic() +
  labs(
    title = "Bin assignment sanity check",
    subtitle = sprintf("Plotted %s / %s bins", format(nrow(plot_bins), big.mark=","), format(nrow(bins), big.mark=",")),
    x = "x_center",
    y = "y_center",
    color = "stem"
  )

# Save
OUT_PNG <- sub("\\.tsv$", "_check.png", OUT_TSV)
ggsave(OUT_PNG, p, width = 7, height = 7, dpi = 300)
cat("Wrote:", OUT_PNG, "\n")


# -----------------------------
# Also write stem-specific bin_metadata tables for downstream radial ordering
# -----------------------------
# bins currently contains: bin_id, x_center, y_center, stem (and possibly more columns if BIN_META had them)

# If you want to keep ALL original columns from BIN_META, ensure bins is the full table.
# Here we assume bins is already the full BIN_META with an added 'stem' column.

fwrite(bins[stem == "stem1"], "bin_metadata_with_cluster_raw_stem1.tsv", sep="\t")
fwrite(bins[stem == "stem2"], "bin_metadata_with_cluster_raw_stem2.tsv", sep="\t")

cat("stem1 n=", nrow(bins[stem=="stem1"]), "\n", sep="")
cat("stem2 n=", nrow(bins[stem=="stem2"]), "\n", sep="")
