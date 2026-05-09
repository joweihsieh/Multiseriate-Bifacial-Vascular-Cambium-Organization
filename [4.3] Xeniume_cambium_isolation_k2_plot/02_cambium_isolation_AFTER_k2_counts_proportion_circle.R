

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
  library(readxl)
})

setwd("/home/woodydrylab/FileShare/20260121_Xenium")

############################################################
# user settings
############################################################
XENIUM_BASE <- "/home/woodydrylab/FileShare/20260121_Xenium"
XLSX_PATH   <- "/home/woodydrylab/FileShare/20260121_Xenium/Cambium_k_single_furtherk2.xlsx"

OUTDIR <- "cambium_before_after_k2_windrose_auto"


dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

OUT_SUMMARY_TSV  <- file.path(OUTDIR, "cambium_k2_fusiform_ray_summary.tsv")
OUT_SUMMARY_XLSX <- file.path(OUTDIR, "cambium_k2_fusiform_ray_summary.xlsx")

OUT_COUNT_PNG <- file.path(OUTDIR, "circular_k2_fusiform_ray_counts.png")
OUT_COUNT_PDF <- file.path(OUTDIR, "circular_k2_fusiform_ray_counts.pdf")

OUT_PROP_PNG <- file.path(OUTDIR, "circular_k2_fusiform_ray_proportion.png")
OUT_PROP_PDF <- file.path(OUTDIR, "circular_k2_fusiform_ray_proportion.pdf")

OUT_COUNT_no_PNG <- file.path(OUTDIR, "circular_k2_fusiform_ray_counts_nolabel.png")
OUT_COUNT_no_PDF <- file.path(OUTDIR, "circular_k2_fusiform_ray_counts_nolabel.pdf")

OUT_PROP_no_PNG <- file.path(OUTDIR, "circular_k2_fusiform_ray_proportion_nolabel.png")
OUT_PROP_no_PDF <- file.path(OUTDIR, "circular_k2_fusiform_ray_proportion_nolabel.pdf")

FUSIFORM_COL <- "#E87D72"
RAY_COL      <- "#56BCC2"

############################################################
# helper
############################################################
pick_first_existing <- function(x, candidates) {
  hit <- candidates[candidates %in% x]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_bin_meta <- function(path) {
  dt <- fread(path)

  xcol <- pick_first_existing(
    names(dt),
    c("x", "X", "x_centroid", "x_center", "center_x", "pxl_col_in_fullres")
  )
  ycol <- pick_first_existing(
    names(dt),
    c("y", "Y", "y_centroid", "y_center", "center_y", "pxl_row_in_fullres")
  )
  ccol <- pick_first_existing(
    names(dt),
    c("cluster_raw", "cluster", "seurat_clusters", "kmeans_cluster")
  )

  if (is.na(xcol) || is.na(ycol) || is.na(ccol)) {
    stop("Cannot detect x/y/cluster columns in: ", path)
  }

  out <- dt[, .(
    x = as.numeric(get(xcol)),
    y = as.numeric(get(ycol)),
    cluster_k2 = as.character(get(ccol))
  )]

  out <- out[is.finite(x) & is.finite(y)]
  out
}

resolve_label_splitstem <- function(label_raw, subset_dir) {
  label_raw <- trimws(as.character(label_raw))

  if (grepl(",", label_raw, fixed = TRUE)) {
    sp <- trimws(strsplit(label_raw, ",", fixed = TRUE)[[1]])

    if (grepl("subset_stem1_", subset_dir)) {
      return(sp[1])
    }

    if (grepl("subset_stem2_", subset_dir)) {
      if (length(sp) >= 2) return(sp[2])
      return(paste0(sp[1], "_stem2"))
    }

    return(label_raw)
  }

  if (grepl("subset_stem1_", subset_dir)) {
    return(paste0(label_raw, "_stem1"))
  }

  if (grepl("subset_stem2_", subset_dir)) {
    return(paste0(label_raw, "_stem2"))
  }

  label_raw
}
make_target_table_one_row <- function(sample_id, label_raw, grid_dir, k, cambium_single) {
  k_fmt <- sprintf("%02d", as.integer(k))
  c_fmt <- sprintf("%02d", as.integer(cambium_single))

  base_dir <- file.path(XENIUM_BASE, sample_id, grid_dir)

  cand_single <- file.path(
    base_dir,
    paste0("subset_from_k", k_fmt, "_c", c_fmt),
    "kmeans_subset_k02_raw_out",
    "bin_metadata_with_cluster_raw.tsv"
  )

  cand_stem1 <- file.path(
    base_dir,
    paste0("subset_stem1_from_k", k_fmt, "_c", c_fmt),
    "kmeans_subset_k02_raw_out",
    "bin_metadata_with_cluster_raw.tsv"
  )

  cand_stem2 <- file.path(
    base_dir,
    paste0("subset_stem2_from_k", k_fmt, "_c", c_fmt),
    "kmeans_subset_k02_raw_out",
    "bin_metadata_with_cluster_raw.tsv"
  )

  out_list <- list()

  if (file.exists(cand_single)) {
    out_list[[length(out_list) + 1]] <- data.table(
      sample_id = sample_id,
      label_raw = label_raw,
      stem_part = "single",
      subset_dir = basename(dirname(dirname(cand_single))),
      bin_meta_path = cand_single
    )
  }

  if (file.exists(cand_stem1)) {
    out_list[[length(out_list) + 1]] <- data.table(
      sample_id = sample_id,
      label_raw = label_raw,
      stem_part = "stem1",
      subset_dir = basename(dirname(dirname(cand_stem1))),
      bin_meta_path = cand_stem1
    )
  }

  if (file.exists(cand_stem2)) {
    out_list[[length(out_list) + 1]] <- data.table(
      sample_id = sample_id,
      label_raw = label_raw,
      stem_part = "stem2",
      subset_dir = basename(dirname(dirname(cand_stem2))),
      bin_meta_path = cand_stem2
    )
  }

  if (length(out_list) == 0) return(NULL)
  rbindlist(out_list, fill = TRUE)
}

make_one_summary <- function(sample_id, stem_part, subset_dir, label, bin_meta_path) {
  bins <- read_bin_meta(bin_meta_path)

  tab <- bins[, .N, by = cluster_k2][order(-N, cluster_k2)]

  if (nrow(tab) != 2) {
    stop("Expected exactly 2 clusters in: ", bin_meta_path,
         " ; found ", nrow(tab))
  }

  fusiform_cluster <- as.character(tab$cluster_k2[1])
  ray_cluster      <- as.character(tab$cluster_k2[2])

  fusiform_bins <- tab$N[1]
  ray_bins      <- tab$N[2]
  total_bins    <- fusiform_bins + ray_bins

  data.table(
    sample_id = sample_id,
    stem_part = stem_part,
    subset_dir = subset_dir,
    label = label,
    fusiform_cluster = fusiform_cluster,
    ray_cluster = ray_cluster,
    fusiform_bins = fusiform_bins,
    ray_bins = ray_bins,
    total_bins = total_bins,
    fusiform_prop = fusiform_bins / total_bins,
    ray_prop = ray_bins / total_bins,
    bin_meta_path = bin_meta_path
  )
}

############################################################
# read Excel and build exact target list
############################################################
cfg <- as.data.table(read_excel(XLSX_PATH))

required_cols <- c("sample_id", "label", "um", "k", "cambium_single")
missing_cols <- setdiff(required_cols, names(cfg))
if (length(missing_cols) > 0) {
  stop("Excel 缺少必要欄位: ", paste(missing_cols, collapse = ", "))
}

cfg[, grid_dir := fifelse(
  grepl("^grid[0-9]+um_out$", as.character(um)),
  as.character(um),
  paste0("grid", sprintf("%02d", as.integer(gsub("[^0-9]", "", as.character(um)))), "um_out")
)]

target_list <- list()

for (i in seq_len(nrow(cfg))) {
  one <- make_target_table_one_row(
    sample_id = cfg$sample_id[i],
    label_raw = cfg$label[i],
    grid_dir = cfg$grid_dir[i],
    k = cfg$k[i],
    cambium_single = cfg$cambium_single[i]
  )

  if (is.null(one)) {
    warning("No matching subset k02 result found for: ",
            cfg$sample_id[i], " | ", cfg$grid_dir[i],
            " | k=", cfg$k[i], " | c=", cfg$cambium_single[i])
    next
  }

  target_list[[length(target_list) + 1]] <- one
}

if (length(target_list) == 0) {
  stop("No valid target k02 files found based on Excel.")
}

targets <- rbindlist(target_list, fill = TRUE)
targets <- unique(targets, by = c("sample_id", "stem_part", "subset_dir", "bin_meta_path"))

targets[, has_comma_label := grepl(",", as.character(label_raw), fixed = TRUE)]

targets <- targets[
  !(has_comma_label == TRUE & stem_part == "single")
]

targets[, label := mapply(resolve_label_splitstem, label_raw, subset_dir)]

targets[, has_comma_label := NULL]

message("Targets kept from Excel: ", nrow(targets))
print(targets[, .(sample_id, stem_part, label, subset_dir)])

############################################################
# summarize each target
############################################################
summary_list <- lapply(seq_len(nrow(targets)), function(i) {
  message("Processing: ", targets$sample_id[i], " | ",
          targets$stem_part[i], " | ", targets$subset_dir[i])

  make_one_summary(
    sample_id = targets$sample_id[i],
    stem_part = targets$stem_part[i],
    subset_dir = targets$subset_dir[i],
    label = targets$label[i],
    bin_meta_path = targets$bin_meta_path[i]
  )
})

final_dt <- rbindlist(summary_list, fill = TRUE)
setorder(final_dt, label, stem_part)

############################################################
# write summary
############################################################
fwrite(final_dt, OUT_SUMMARY_TSV, sep = "\t", quote = FALSE, na = "")

wb <- createWorkbook()
addWorksheet(wb, "k2_fusiform_ray_summary")
writeData(wb, "k2_fusiform_ray_summary", final_dt)
saveWorkbook(wb, OUT_SUMMARY_XLSX, overwrite = TRUE)

############################################################
# count wind-rose
############################################################
dt_count <- copy(final_dt)
dt_count[, sector_id := .I]

plot_dt_count <- rbindlist(list(
  dt_count[, .(label, sector_id, group = "ray", value = ray_bins)],
  dt_count[, .(label, sector_id, group = "fusiform", value = fusiform_bins)]
))
plot_dt_count[, group := factor(group, levels = c("ray", "fusiform"))]

shared_max <- max(dt_count$total_bins, na.rm = TRUE)
if (!is.finite(shared_max) || shared_max <= 0) shared_max <- 1

shared_breaks <- pretty(c(0, shared_max), n = 4)
shared_breaks <- shared_breaks[shared_breaks >= 0]
shared_breaks <- unique(shared_breaks)

label_radius_count <- max(shared_breaks) * 1.18

label_dt_count <- copy(dt_count)
label_dt_count[, x := sector_id]
label_dt_count[, y := label_radius_count]
label_dt_count[, angle := 90 - 360 * (sector_id - 0.5) / .N]
label_dt_count[, hjust := ifelse(angle < -90, 1, 0)]
label_dt_count[angle < -90, angle := angle + 180]

guide_dt_count <- data.table(y = shared_breaks[shared_breaks > 0])
guide_label_dt_count <- data.table(
  x = 0.55,
  y = shared_breaks[shared_breaks > 0],
  label = as.character(shared_breaks[shared_breaks > 0])
)

p_count <- ggplot(plot_dt_count, aes(x = sector_id, y = value, fill = group)) +
  geom_hline(
    data = guide_dt_count,
    aes(yintercept = y),
    inherit.aes = FALSE,
    color = "grey80",
    linewidth = 0.4
  ) +
  geom_col(width = 0.95, color = "white", linewidth = 0.6) +
  geom_text(
    data = label_dt_count,
    aes(x = x, y = y, label = label, angle = angle, hjust = hjust),
    inherit.aes = FALSE,
    size = 3.6
  ) +
  geom_text(
    data = guide_label_dt_count,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3,
    hjust = 0,
    color = "grey30"
  ) +
  coord_polar(theta = "x", clip = "off") +
  scale_x_continuous(
    breaks = dt_count$sector_id,
    labels = NULL,
    limits = c(0.5, nrow(dt_count) + 0.5)
  ) +
  scale_y_continuous(
    limits = c(0, label_radius_count * 1.05),
    breaks = shared_breaks,
    labels = NULL,
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = c(
    ray = RAY_COL,
    fusiform = FUSIFORM_COL
  )) +
  labs(fill = NULL) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "top",
    plot.margin = margin(20, 80, 20, 80)
  )

ggsave(OUT_COUNT_PNG, p_count, width = 10, height = 10, dpi = 300)
ggsave(OUT_COUNT_PDF, p_count, width = 10, height = 10)

############################################################
# proportion wind-rose
############################################################
dt_prop <- copy(final_dt)
dt_prop[, sector_id := .I]

plot_dt_prop <- rbindlist(list(
  dt_prop[, .(label, sector_id, group = "ray", value = ray_prop)],
  dt_prop[, .(label, sector_id, group = "fusiform", value = fusiform_prop)]
))
plot_dt_prop[, group := factor(group, levels = c("ray", "fusiform"))]

label_radius_prop <- 1.12

label_dt_prop <- copy(dt_prop)
label_dt_prop[, x := sector_id]
label_dt_prop[, y := label_radius_prop]
label_dt_prop[, angle := 90 - 360 * (sector_id - 0.5) / .N]
label_dt_prop[, hjust := ifelse(angle < -90, 1, 0)]
label_dt_prop[angle < -90, angle := angle + 180]

guide_dt_prop <- data.table(y = c(0.25, 0.5, 0.75, 1.0))
guide_label_dt_prop <- data.table(
  x = 0.55,
  y = c(0.25, 0.5, 0.75, 1.0),
  label = c("25%", "50%", "75%", "100%")
)

p_prop <- ggplot(plot_dt_prop, aes(x = sector_id, y = value, fill = group)) +
  geom_hline(
    data = guide_dt_prop,
    aes(yintercept = y),
    inherit.aes = FALSE,
    color = "grey80",
    linewidth = 0.4
  ) +
  geom_col(width = 0.95, color = "white", linewidth = 0.6) +
  geom_text(
    data = label_dt_prop,
    aes(x = x, y = y, label = label, angle = angle, hjust = hjust),
    inherit.aes = FALSE,
    size = 3.6
  ) +
  geom_text(
    data = guide_label_dt_prop,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3,
    hjust = 0,
    color = "grey30"
  ) +
  coord_polar(theta = "x", clip = "off") +
  scale_x_continuous(
    breaks = dt_prop$sector_id,
    labels = NULL,
    limits = c(0.5, nrow(dt_prop) + 0.5)
  ) +
  scale_y_continuous(
    limits = c(0, 1.18),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = NULL,
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = c(
    ray = RAY_COL,
    fusiform = FUSIFORM_COL
  )) +
  labs(fill = NULL) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "top",
    plot.margin = margin(20, 80, 20, 80)
  )

ggsave(OUT_PROP_PNG, p_prop, width = 10, height = 10, dpi = 300)
ggsave(OUT_PROP_PDF, p_prop, width = 10, height = 10)

message("Done.")
message("Summary TSV : ", OUT_SUMMARY_TSV)
message("Summary XLSX: ", OUT_SUMMARY_XLSX)
message("Count PNG   : ", OUT_COUNT_PNG)
message("Count PDF   : ", OUT_COUNT_PDF)
message("Prop PNG    : ", OUT_PROP_PNG)
message("Prop PDF    : ", OUT_PROP_PDF)



p_no_count <- ggplot(plot_dt_count, aes(x = sector_id, y = value, fill = group)) +
  geom_hline(
    data = guide_dt_count,
    aes(yintercept = y),
    inherit.aes = FALSE,
    color = "grey80",
    linewidth = 0.4
  ) +
  geom_col(width = 0.95, color = "white", linewidth = 0.6) +
  coord_polar(theta = "x", clip = "off") +
  scale_x_continuous(
    breaks = dt_count$sector_id,
    labels = NULL,
    limits = c(0.5, nrow(dt_count) + 0.5)
  ) +
  scale_y_continuous(
    limits = c(0, label_radius_count * 1.05),
    breaks = shared_breaks,
    labels = NULL,
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = c(
    ray = RAY_COL,
    fusiform = FUSIFORM_COL
  )) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "top",
    plot.margin = margin(20, 20, 20, 20)
  )


p_no_prop <- ggplot(plot_dt_prop, aes(x = sector_id, y = value, fill = group)) +
  geom_hline(
    data = guide_dt_prop,
    aes(yintercept = y),
    inherit.aes = FALSE,
    color = "grey80",
    linewidth = 0.4
  ) +
  geom_col(width = 0.95, color = "white", linewidth = 0.6) +
  coord_polar(theta = "x", clip = "off") +
  scale_x_continuous(
    breaks = dt_prop$sector_id,
    labels = NULL,
    limits = c(0.5, nrow(dt_prop) + 0.5)
  ) +
  scale_y_continuous(
    limits = c(0, 1.18),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = NULL,
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = c(
    ray = RAY_COL,
    fusiform = FUSIFORM_COL
  )) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "top",
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave(OUT_COUNT_no_PNG, p_no_count, width = 10, height = 10, dpi = 300)
ggsave(OUT_COUNT_no_PDF, p_no_count, width = 10, height = 10)

ggsave(OUT_PROP_no_PNG, p_no_prop, width = 10, height = 10, dpi = 300)
ggsave(OUT_PROP_no_PDF, p_no_prop, width = 10, height = 10)
