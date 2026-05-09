
######################### all in one with dots
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE_DIR <- "/home/woodydrylab/FileShare/20260121_Xenium"

files <- list.files(
  BASE_DIR,
  pattern = "transcriptome_distance_raw_anchor_vs_push.tsv",
  recursive = TRUE,
  full.names = TRUE
)

cat("Found", length(files), "files\n")
if (length(files) == 0) {
  stop("No transcriptome_distance_raw_anchor_vs_push.tsv files found.")
}

files <- files[c(1:3,5:11)]

all_dt <- rbindlist(
  lapply(files, function(f) {
    dt <- fread(f)
    dt[, sample := basename(dirname(f))]
    dt
  }),
  fill = TRUE
)

stopifnot(all(c("step_um", "transcriptome_cosine_dist_raw") %in% names(all_dt)))

all_dt <- all_dt[!is.na(step_um) & !is.na(transcriptome_cosine_dist_raw)]
all_dt[, step_um := as.character(step_um)]

step_levels <- sort(unique(as.numeric(all_dt$step_um)))
all_dt[, step_um := factor(step_um, levels = as.character(step_levels))]

print(table(all_dt$step_um))


# ============================================================
# Statistics
# ============================================================
kw <- kruskal.test(transcriptome_cosine_dist_raw ~ step_um, data = all_dt)

pw <- pairwise.wilcox.test(
  x = all_dt$transcriptome_cosine_dist_raw,
  g = all_dt$step_um,
  p.adjust.method = "BH",
  exact = FALSE
)

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return("< 1e-4")
  sprintf("%.3g", p)
}

p_to_stars <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

ann_list <- list()
pw_mat <- pw$p.value

if (!is.null(pw_mat)) {
  rn <- rownames(pw_mat)
  cn <- colnames(pw_mat)
  idx <- 1L

  for (i in seq_along(rn)) {
    for (j in seq_along(cn)) {
      pval <- pw_mat[i, j]
      if (!is.na(pval)) {
        ann_list[[idx]] <- data.table(
          group1 = rn[i],
          group2 = cn[j],
          p_adj = pval
        )
        idx <- idx + 1L
      }
    }
  }
}

if (length(ann_list) > 0) {
  ann_dt <- rbindlist(ann_list)
} else {
  ann_dt <- data.table(
    group1 = character(),
    group2 = character(),
    p_adj = numeric()
  )
}

y_max <- max(all_dt$transcriptome_cosine_dist_raw, na.rm = TRUE)
y_min <- min(all_dt$transcriptome_cosine_dist_raw, na.rm = TRUE)
y_range <- max(y_max - y_min, 1e-6)

if (nrow(ann_dt) > 0) {
  ann_dt[, x1 := match(group1, levels(all_dt$step_um))]
  ann_dt[, x2 := match(group2, levels(all_dt$step_um))]
  ann_dt <- ann_dt[!is.na(x1) & !is.na(x2)]
  ann_dt <- ann_dt[order(x1, x2)]
  ann_dt[, y := y_max + seq_len(.N) * (0.08 * y_range)]
  ann_dt[, label := paste0(
    vapply(p_adj, p_to_stars, character(1)),
    " (BH p=",
    vapply(p_adj, format_p, character(1)),
    ")"
  )]
}

# ============================================================
# Plot
# ============================================================
p <- ggplot(all_dt, aes(x = step_um, y = transcriptome_cosine_dist_raw)) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    height = 0,
    alpha = 0.25,
    size = 1
  ) +
  theme_classic() +
  labs(
    x = "Radial push distance (µm)",
    y = "Transcriptome cosine distance",
    title = "Anchor vs radial push transcriptome distance",
    subtitle = paste0("Kruskal-Wallis p = ", format_p(kw$p.value))
  )

if (nrow(ann_dt) > 0 && "label" %in% names(ann_dt)) {
  for (i in seq_len(nrow(ann_dt))) {
    seg_dt <- ann_dt[i]

    p <- p +
      geom_segment(
        data = seg_dt,
        aes(x = x1, xend = x2, y = y, yend = y),
        inherit.aes = FALSE
      ) +
      geom_segment(
        data = seg_dt,
        aes(x = x1, xend = x1, y = y - 0.015 * y_range, yend = y),
        inherit.aes = FALSE
      ) +
      geom_segment(
        data = seg_dt,
        aes(x = x2, xend = x2, y = y - 0.015 * y_range, yend = y),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = seg_dt,
        mapping = aes(
          x = (x1 + x2) / 2,
          y = y + 0.015 * y_range,
          label = label
        ),
        inherit.aes = FALSE,
        size = 3
      )
  }

  p <- p + coord_cartesian(
    ylim = c(y_min, max(ann_dt$y) + 0.12 * y_range)
  )
}

out_png <- file.path(BASE_DIR, "ALL_samples_transcriptome_distance_boxplot_with_points_stats.png")
ggsave(
  out_png,
  p,
  width = 7,
  height = 5.5,
  dpi = 600
)

# ============================================================
# Save stats
# ============================================================
out_txt <- file.path(BASE_DIR, "ALL_samples_transcriptome_distance_stats.txt")
zz <- file(out_txt, open = "wt")
writeLines("Kruskal-Wallis test", zz)
writeLines(capture.output(kw), zz)
writeLines("", zz)
writeLines("Pairwise Wilcoxon test (BH adjusted)", zz)
writeLines(capture.output(pw), zz)
close(zz)

cat("Saved plot:", out_png, "\n")
cat("Saved stats:", out_txt, "\n")

############### replot all boxplots again
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE_DIR <- "/home/woodydrylab/FileShare/20260121_Xenium"

files <- list.files(
  BASE_DIR,
  pattern = "transcriptome_distance_raw_anchor_vs_push.tsv",
  recursive = TRUE,
  full.names = TRUE
)

cat("Found files:\n")
print(files)

for (f in files) {

  dt <- fread(f)

  sample_name <- basename(dirname(f))

  p <- ggplot(dt, aes(x = factor(step_um), y = transcriptome_cosine_dist_raw, fill = factor(step_um))) +
    geom_boxplot(outlier.size = 0.4) +
    scale_fill_manual(
        values = c(
        "5"  = "#1f78b4",
        "10" = "#33a02c",
        "15" = "#e31a1c"
        )
      ) +
    theme_classic() +
    labs(
      x = "step (um)",
      y = "Cosine distance",
      title = sample_name
    )

  out_png <- file.path(dirname(f), "transcriptome_cosine_distance_raw_by_step.png")

  ggsave(
    out_png,
    p,
    width = 5.8,
    height = 4.2,
    dpi = 600
  )

  cat("Saved:", out_png, "\n")
}