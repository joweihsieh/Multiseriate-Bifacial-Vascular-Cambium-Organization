#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

IN_TSV <- "/home/woodydrylab/FileShare/20260121_Xenium/Xenium_full_QC_summary_50um.tsv"
OUT_PNG <- "/home/woodydrylab/FileShare/20260121_Xenium/Xenium_50um_gmean_boxplot.png"
OLD_VALUE <- 25.91
USE_REGION <- NULL

dt <- fread(IN_TSV)
dt <- dt[is.finite(g_mean) & !is.na(g_mean)]

if (!is.null(USE_REGION)) {
  dt <- dt[region == USE_REGION]
}

if (nrow(dt) == 0) {
  stop("No data to plot.")
}

ymax <- max(dt$g_mean, na.rm = TRUE)
ymax <- max(ymax, OLD_VALUE + 10)

png(OUT_PNG, width = 800, height = 1200, res = 150)

par(
  mar = c(6, 5, 4, 2) + 0.1,
  font.lab = 2,
  cex.lab = 1.8,
  cex.axis = 1.5,
  cex.main = 1.8
)

boxplot(
  dt$g_mean,
  outline = FALSE,
  names = "",
  ylim = c(0, ymax),
  ylab = "#Genes per 50 µm bin",
  main = "",
  lwd = 2
)

points(
  jitter(rep(1, nrow(dt)), amount = 0.05),
  dt$g_mean,
  pch = 16,
  cex = 1.8
)

abline(
  h = OLD_VALUE,
  lty = 2,
  lwd = 3,
  col = "red"
)

arrows(
  x0 = 1,
  y0 = OLD_VALUE + 6,
  x1 = 1,
  y1 = OLD_VALUE,
  length = 0.15,
  lwd = 3,
  col = "red"
)

text(
  x = 1,
  y = OLD_VALUE + 7,
  labels = paste0("Previous Stereo-seq = ", OLD_VALUE),
  pos = 3,
  col = "red",
  cex = 1.6
)

dev.off()

cat("Saved to:", OUT_PNG, "\n")
