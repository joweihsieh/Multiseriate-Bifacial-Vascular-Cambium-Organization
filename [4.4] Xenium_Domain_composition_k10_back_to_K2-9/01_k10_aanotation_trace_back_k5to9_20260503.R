suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

# ============================================================
# 0. Input / Output
# ============================================================
IN_FILE <- "/Users/jo-wei/Library/CloudStorage/Dropbox/YCL/Spatial_Transcriptomes/Xenium/20260121_Data/20260208_grid_based/5um_k5_to_k10_domain_composition/k10_backtrace_to_k5_k10_wide_prop_within_cluster.tsv"

OUT_DIR <- dirname(IN_FILE)

if (!file.exists(IN_FILE)) {
  stop("Input file does not exist: ", IN_FILE)
}

# ============================================================
# 1. Read wide proportion table
# ============================================================
dt <- fread(IN_FILE)

if ("V7" %in% names(dt)) {
  setnames(dt, "V7", "unassigned")
}

# ============================================================
# 2. Domain settings
# ============================================================
DOMAIN_LEVELS <- c(
  "unassigned",
  "epidermis",
  "phloem",
  "cambium",
  "parenchyma",
  "sclerenchyma",
  "background"
)

DOMAIN_COLORS <- c(
  unassigned   = "grey30",
  epidermis    = "#9590FF",
  phloem       = "#00B0F6",
  cambium      = "#F8766D",
  parenchyma   = "#39B600",
  sclerenchyma = "#D89000",
  background   = "grey30"
)

domain_cols <- intersect(DOMAIN_LEVELS, names(dt))

if (length(domain_cols) == 0) {
  stop("No domain proportion columns found in the input table.")
}

# ============================================================
# 3. Find majority domain for each cluster
# ============================================================
dt_majority <- melt(
  dt,
  id.vars = setdiff(names(dt), domain_cols),
  measure.vars = domain_cols,
  variable.name = "majority_domain",
  value.name = "majority_prop"
)

dt_majority <- dt_majority[
  dt_majority[, .I[which.max(majority_prop)],
              by = .(sample_id, sample_label, k_now, cluster_now)]$V1
]

dt_majority[is.na(majority_domain), majority_domain := "background"]
dt_majority[is.na(majority_prop), majority_prop := 0]

dt_majority[, majority_domain := factor(as.character(majority_domain), levels = DOMAIN_LEVELS)]

# ============================================================
# 4. Labels and polar geometry
# ============================================================
dt_majority[, k_label := factor(
  paste0("K = ", k_now),
  levels = paste0("K = ", sort(unique(k_now)))
)]

dt_majority[, cluster_label := paste0("C", sprintf("%02d", cluster_now))]

setorder(dt_majority, sample_label, k_now, majority_domain, cluster_now)

dt_majority[, N := .N, by = .(sample_label, k_now)]
dt_majority[, x_polar := (seq_len(.N) - 0.5) / N, by = .(sample_label, k_now)]
dt_majority[, bar_width := 0.85 / N]

dt_majority[, text_angle := 90 - (x_polar * 360)]
dt_majority[text_angle < -90, text_angle := text_angle + 180]

# ============================================================
# 5. Fill color setting
#    majority_prop < 0.65 -> dark gray
# ============================================================
dt_majority[, fill_group := as.character(majority_domain)]
dt_majority[majority_prop < 0.65, fill_group := "low_purity"]

FILL_COLORS <- c(
  DOMAIN_COLORS,
  low_purity = "grey"
)

FILL_BREAKS <- c(
  "unassigned",
  "epidermis",
  "phloem",
  "cambium",
  "parenchyma",
  "sclerenchyma",
  "background",
  "low_purity"
)

FILL_LABELS <- c(
  unassigned   = "unassigned",
  epidermis    = "epidermis",
  phloem       = "phloem",
  cambium      = "cambium",
  parenchyma   = "parenchyma",
  sclerenchyma = "sclerenchyma",
  background   = "background",
  low_purity   = "< 0.65"
)

# ============================================================
# 6. Base rose chart
# ============================================================
p_base <- ggplot(dt_majority, aes(x = x_polar, y = majority_prop)) +
  geom_col(
    aes(fill = fill_group, width = bar_width),
    color = "black",
    linewidth = 0.3,
    alpha = 0.85
  ) +
  coord_polar(theta = "x", start = 0) +
  facet_grid(k_label ~ sample_label) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-0.3, 1.3)) +
  scale_fill_manual(
    values = FILL_COLORS,
    breaks = FILL_BREAKS,
    labels = FILL_LABELS[FILL_BREAKS],
    name = "Majority domain",
    drop = FALSE
  ) +
  theme_void(base_size = 12) +
  labs(
    title = "Clustering Consistency Matrix (Nightingale Rose Chart)",
    subtitle = "Rows: Resolution (K) | Columns: Samples\nBar height represents the purity (%) of the majority domain. Clusters with majority proportion < 0.65 are shown in dark gray."
  ) +
  theme(
    strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.5),
    strip.text.x = element_text(face = "bold", size = 11, margin = margin(b = 5, t = 5)),
    strip.text.y = element_text(face = "bold", size = 11, angle = 0, margin = margin(l = 5, r = 5)),
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 5, t = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50", margin = margin(b = 15)),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.spacing = unit(1.2, "lines")
  ) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.85))
  )

# ============================================================
# 7. With cluster labels
# ============================================================
p_rose_rotated <- p_base +
  geom_text(
    aes(
      y = majority_prop + 0.15,
      label = cluster_label,
      angle = text_angle
    ),
    size = 3.2,
    fontface = "bold",
    color = "gray20"
  )

# ============================================================
# 8. Without cluster labels
# ============================================================
p_rose_rotated_no_text <- p_base

# ============================================================
# 9. Save only two PNG files
# ============================================================
n_samples <- uniqueN(dt_majority$sample_label)
n_k <- uniqueN(dt_majority$k_now)

plot_width <- max(10, 2.5 + 2.5 * n_samples)
plot_height <- max(8, 2.5 + 2.2 * n_k)

ggsave(
  filename = file.path(OUT_DIR, "consistency_rose_chart_rotated.png"),
  plot = p_rose_rotated,
  width = plot_width,
  height = plot_height,
  dpi = 300
)

ggsave(
  filename = file.path(OUT_DIR, "consistency_rose_chart_rotated.no.text.png"),
  plot = p_rose_rotated_no_text,
  width = plot_width,
  height = plot_height,
  dpi = 300
)

cat("Done.\n")
cat("Saved:\n")
cat(file.path(OUT_DIR, "consistency_rose_chart_rotated.png"), "\n")
cat(file.path(OUT_DIR, "consistency_rose_chart_rotated.no.text.png"), "\n")