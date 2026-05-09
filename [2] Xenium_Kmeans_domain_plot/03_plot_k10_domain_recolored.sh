#!/usr/bin/env bash
# Re-plot one sample's existing 5 um / k = 10 K-means result
# using manual domain annotations and fixed colors.
# Make sure you are inside the Xenium output folder first
# conda activate r_spatial
# cd /home/woodydrylab/FileShare/20260121_Xenium/output-XXXX...

for um in 5
do
  grid_dir=$(printf "grid%02dum_out" "$um")

  if [ ! -d "${grid_dir}" ]; then
    echo "[skip] ${grid_dir} not found"
    continue
  fi

  echo "=============================="
  echo "Grid = ${um} um  (${grid_dir})"
  echo "=============================="

  for k in 10
  do
    echo "  Running k = $k"
    outdir="${grid_dir}/kmeans_k${k}_raw_out"

    Rscript - <<EOF
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(readxl)
})

outdir <- "${outdir}"
k_run <- ${k}
bin_um <- ${um}

PAIR_XLSX   <- "/home/woodydrylab/FileShare/20260121_Xenium/Cambium_Ray_pairs.xlsx"
DOMAIN_XLSX <- "/home/woodydrylab/FileShare/20260121_Xenium/k10_5domain.xlsx"

sample_id_run <- basename(getwd())
grid_name <- sprintf("grid%02dum_out", bin_um)

cat("sample_id_run =", sample_id_run, "\\n")
cat("grid_name     =", grid_name, "\\n")
cat("k_run         =", k_run, "\\n")
cat("PAIR_XLSX     =", PAIR_XLSX, "\\n")
cat("DOMAIN_XLSX   =", DOMAIN_XLSX, "\\n")
cat("outdir        =", outdir, "\\n")

if (!dir.exists(outdir)) {
  stop("Output directory not found: ", outdir)
}
if (!file.exists(PAIR_XLSX)) {
  stop("PAIR_XLSX not found: ", PAIR_XLSX)
}
if (!file.exists(DOMAIN_XLSX)) {
  stop("DOMAIN_XLSX not found: ", DOMAIN_XLSX)
}

# ============================================================
# fixed domain colors
# ============================================================
domain_colors <- c(
  "sclerenchyma" = "#D89000",
  "epidermis"    = "#9590FF",
  "phloem"       = "#00B0F6",
  "parenchyma"   = "#39B600",
  "cambium"      = "#F8766D",
  "background"   = "#333333"
)

background_color <- "black"

# ============================================================
# helper
# ============================================================
std_colnames <- function(dt) {
  nms <- names(dt)
  nms2 <- tolower(nms)
  nms2 <- gsub("[[:space:]]+", "_", nms2)
  nms2 <- gsub("[^a-z0-9_]", "_", nms2)
  nms2 <- gsub("_+", "_", nms2)
  nms2 <- gsub("^_|_$", "", nms2)
  setnames(dt, old = nms, new = nms2)
  dt
}

rename_first_match <- function(dt, target, candidates) {
  hit <- intersect(candidates, names(dt))
  if (!(target %in% names(dt)) && length(hit) >= 1) {
    setnames(dt, hit[1], target)
  }
  dt
}

parse_pair <- function(x) {
  if (length(x) == 0) return(NULL)

  x <- as.character(x[[1]])
  x <- trimws(x)

  if (is.na(x) || x == "" || x %in% c("NA", "NaN")) {
    return(NULL)
  }

  vals <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  vals <- vals[vals != ""]
  if (length(vals) == 0) return(NULL)

  vals
}

make_bg_theme <- function() {
  theme_void() +
    theme(
      plot.background   = element_rect(fill = "black", color = NA),
      panel.background  = element_rect(fill = "black", color = NA),
      legend.background = element_rect(fill = "black", color = NA),
      legend.key        = element_rect(fill = "black", color = NA),
      legend.text       = element_text(color = "white"),
      legend.title      = element_text(color = "white"),
      plot.title        = element_text(color = "white", hjust = 0.5)
    )
}

detect_stem_modes <- function(outdir, sample_id_run) {
  meta_default <- file.path(outdir, "bin_metadata_with_cluster_raw.tsv")
  meta_stem1   <- file.path(outdir, "bin_metadata_with_cluster_raw_stem1.tsv")
  meta_stem2   <- file.path(outdir, "bin_metadata_with_cluster_raw_stem2.tsv")

  if (file.exists(meta_stem1) && file.exists(meta_stem2)) {
    return(data.table(
      stem_tag = c("stem1", "stem2"),
      meta_file = c(meta_stem1, meta_stem2),
      out_suffix = c("_stem1", "_stem2"),
      label_single = c("stem1", "stem2")
    ))
  }

  if (file.exists(meta_default)) {
    return(data.table(
      stem_tag = "default",
      meta_file = meta_default,
      out_suffix = "",
      label_single = sample_id_run
    ))
  }

  stop("No metadata file found in outdir: ", outdir)
}

find_best_pair_row <- function(pair_sub, label_single, stem_tag) {
  if (nrow(pair_sub) == 0) return(pair_sub)

  if ("label" %in% names(pair_sub)) {
    lab <- trimws(as.character(pair_sub\$label))

    hit <- pair_sub[tolower(lab) == tolower(label_single)]
    if (nrow(hit) > 0) return(hit[1])

    if (stem_tag %in% c("stem1", "stem2")) {
      hit <- pair_sub[grepl(paste0(stem_tag, "\$"), lab, ignore.case = TRUE)]
      if (nrow(hit) > 0) return(hit[1])
    }
  }

  pair_sub[1]
}

find_best_domain_rows <- function(domain_sub_all, label_single, stem_tag) {
  if (nrow(domain_sub_all) == 0) return(domain_sub_all)

  if ("label" %in% names(domain_sub_all)) {
    lab <- trimws(as.character(domain_sub_all\$label))

    hit <- domain_sub_all[tolower(lab) == tolower(label_single)]
    if (nrow(hit) > 0) return(hit)

    if (stem_tag %in% c("stem1", "stem2")) {
      hit <- domain_sub_all[grepl(paste0(stem_tag, "\$"), lab, ignore.case = TRUE)]
      if (nrow(hit) > 0) return(hit)
    }
  }

  domain_sub_all
}

# ============================================================
# Read pair Excel
# ============================================================
pair_dt <- as.data.table(read_excel(PAIR_XLSX))
pair_dt <- std_colnames(pair_dt)

pair_dt <- rename_first_match(pair_dt, "sample_id",    c("sample_id", "sample", "sampleid"))
pair_dt <- rename_first_match(pair_dt, "um",           c("um", "grid", "grid_name"))
pair_dt <- rename_first_match(pair_dt, "k",            c("k"))
pair_dt <- rename_first_match(pair_dt, "label",        c("label"))
pair_dt <- rename_first_match(pair_dt, "pair_cambium", c("pair_cambium"))
pair_dt <- rename_first_match(pair_dt, "pair_ray",     c("pair_ray"))

required_pair_cols <- c("sample_id", "um", "k", "pair_cambium", "pair_ray")
if (!all(required_pair_cols %in% names(pair_dt))) {
  stop(
    "Missing required columns in PAIR_XLSX. Required = ",
    paste(required_pair_cols, collapse = ", "),
    " | Found = ",
    paste(names(pair_dt), collapse = ", ")
  )
}

pair_dt[, sample_id := trimws(as.character(sample_id))]
pair_dt[, um := trimws(as.character(um))]
pair_dt[, k := as.integer(k)]
if ("label" %in% names(pair_dt)) pair_dt[, label := trimws(as.character(label))]
pair_dt[, pair_cambium := trimws(as.character(pair_cambium))]
pair_dt[, pair_ray := trimws(as.character(pair_ray))]

pair_sub <- pair_dt[
  sample_id == sample_id_run & um == grid_name & k == k_run
]

if (nrow(pair_sub) == 0) {
  stop(
    "No matching row in PAIR_XLSX for sample_id = ", sample_id_run,
    ", um = ", grid_name,
    ", k = ", k_run
  )
}

cat("Matched pair_sub:\\n")
print(pair_sub)

# ============================================================
# Read domain Excel
# ============================================================
domain_dt <- as.data.table(read_excel(DOMAIN_XLSX))
domain_dt <- std_colnames(domain_dt)

domain_dt <- rename_first_match(domain_dt, "sample_id", c("sample_id", "sample", "sampleid"))
domain_dt <- rename_first_match(domain_dt, "um",        c("um", "grid", "grid_name"))
domain_dt <- rename_first_match(domain_dt, "k",         c("k"))
domain_dt <- rename_first_match(domain_dt, "label",     c("label"))
domain_dt <- rename_first_match(domain_dt, "cluster",   c("cluster", "cluster_raw", "cluster_id"))
domain_dt <- rename_first_match(domain_dt, "domain",    c("domain", "tissue", "annotation"))

required_domain_cols <- c("sample_id", "um", "k", "cluster", "domain")
if (!all(required_domain_cols %in% names(domain_dt))) {
  stop(
    "Missing required columns in DOMAIN_XLSX. Required = ",
    paste(required_domain_cols, collapse = ", "),
    " | Found = ",
    paste(names(domain_dt), collapse = ", ")
  )
}

domain_dt[, sample_id := trimws(as.character(sample_id))]
domain_dt[, um        := trimws(as.character(um))]
domain_dt[, k         := as.integer(k)]
if ("label" %in% names(domain_dt)) domain_dt[, label := trimws(as.character(label))]
domain_dt[, cluster   := trimws(as.character(cluster))]
domain_dt[, domain    := trimws(tolower(as.character(domain)))]
domain_dt[domain %in% c("bg", "backgroud", "backgrond"), domain := "background"]

domain_sub_all <- domain_dt[
  sample_id == sample_id_run & um == grid_name & k == k_run
]

if (nrow(domain_sub_all) == 0) {
  stop(
    "No matching rows in DOMAIN_XLSX for sample_id = ", sample_id_run,
    ", um = ", grid_name,
    ", k = ", k_run
  )
}

bad_domains <- setdiff(unique(domain_sub_all\$domain), names(domain_colors))
if (length(bad_domains) > 0) {
  stop(
    "Unknown domain(s) in DOMAIN_XLSX: ",
    paste(bad_domains, collapse = ", "),
    ". Allowed = ",
    paste(names(domain_colors), collapse = ", ")
  )
}

cat("Matched domain_sub_all:\\n")
print(domain_sub_all)

# ============================================================
# determine stem modes by existing metadata files
# ============================================================
stem_modes <- detect_stem_modes(outdir, sample_id_run)

cat("stem_modes:\\n")
print(stem_modes)

if ("label" %in% names(pair_sub)) {
  labs_pair <- unique(trimws(as.character(pair_sub\$label)))
  labs_pair <- labs_pair[!is.na(labs_pair) & labs_pair != ""]

  if (length(labs_pair) >= 2) {
    hit1 <- labs_pair[grepl("stem1\$", labs_pair, ignore.case = TRUE)]
    hit2 <- labs_pair[grepl("stem2\$", labs_pair, ignore.case = TRUE)]

    if (nrow(stem_modes) == 2) {
      if (length(hit1) >= 1) stem_modes[stem_tag == "stem1", label_single := hit1[1]]
      if (length(hit2) >= 1) stem_modes[stem_tag == "stem2", label_single := hit2[1]]
    }
  }
}

cat("stem_modes after label matching:\\n")
print(stem_modes)

# ============================================================
# main loop per stem
# ============================================================
for (ii in seq_len(nrow(stem_modes))) {

  stem_tag      <- stem_modes\$stem_tag[ii]
  bin_meta_file <- stem_modes\$meta_file[ii]
  out_suffix    <- stem_modes\$out_suffix[ii]
  label_single  <- stem_modes\$label_single[ii]

  cat("\\n====================================\\n")
  cat("Processing stem_tag   =", stem_tag, "\\n")
  cat("Processing label      =", label_single, "\\n")
  cat("Metadata file         =", bin_meta_file, "\\n")
  cat("====================================\\n")

  if (!file.exists(bin_meta_file)) {
    stop("Missing metadata file: ", bin_meta_file)
  }

  # ----------------------------
  # subset pair/domain with safe fallback
  # ----------------------------
  pair_one <- find_best_pair_row(pair_sub, label_single, stem_tag)
  if (nrow(pair_one) == 0) {
    stop("No usable pair row found for stem_tag = ", stem_tag)
  }

  domain_one <- find_best_domain_rows(domain_sub_all, label_single, stem_tag)
  if (nrow(domain_one) == 0) {
    stop("No usable domain rows found for stem_tag = ", stem_tag)
  }

  cat("pair_one:\\n")
  print(pair_one)

  cat("domain_one:\\n")
  print(domain_one)

  if (anyDuplicated(domain_one\$cluster)) {
    dup_clusters <- unique(domain_one\$cluster[duplicated(domain_one\$cluster)])
    stop(
      "Duplicate cluster annotation in DOMAIN_XLSX for stem_tag = ",
      stem_tag, "; clusters = ", paste(dup_clusters, collapse = ", ")
    )
  }

  cluster2domain <- setNames(domain_one\$domain, domain_one\$cluster)
  cluster2color  <- setNames(domain_colors[domain_one\$domain], domain_one\$cluster)

  # ----------------------------
  # parse groups
  # ----------------------------
  cluster_groups <- list()

  grp_cambium <- parse_pair(pair_one\$pair_cambium)
  grp_ray     <- parse_pair(pair_one\$pair_ray)

  if (!is.null(grp_cambium)) cluster_groups[["pair_cambium"]] <- list(grp_cambium)
  if (!is.null(grp_ray))     cluster_groups[["pair_ray"]]     <- list(grp_ray)

  cat("Parsed cluster groups:\\n")
  print(cluster_groups)

  # ----------------------------
  # read metadata
  # ----------------------------
  bin_meta <- fread(bin_meta_file)

  required_meta_cols <- c("x_center", "y_center", "cluster_raw")
  if (!all(required_meta_cols %in% names(bin_meta))) {
    stop(
      "Missing required columns in metadata. Required = ",
      paste(required_meta_cols, collapse = ", "),
      " | Found = ",
      paste(names(bin_meta), collapse = ", ")
    )
  }

  bin_meta[, cluster_raw := trimws(as.character(cluster_raw))]
  #clusters <- sort(unique(bin_meta\$cluster_raw))
  clusters_unique <- unique(bin_meta\$cluster_raw)
  clusters <- clusters_unique[order(as.numeric(clusters_unique))]

  missing_anno <- setdiff(clusters, names(cluster2domain))
  if (length(missing_anno) > 0) {
    warning(
      "These clusters exist in metadata but not in DOMAIN_XLSX for stem_tag = ",
      stem_tag, ": ", paste(missing_anno, collapse = ", "),
      ". They will use fallback colors."
    )
    fallback_pal <- hue_pal()(length(missing_anno))
    names(fallback_pal) <- missing_anno
    cluster2color  <- c(cluster2color, fallback_pal)
    cluster2domain <- c(cluster2domain, setNames(rep("unannotated", length(missing_anno)), missing_anno))
  }

  bin_meta[, cluster_raw := factor(cluster_raw, levels = clusters)]

  # ----------------------------
  # all clusters plot
  # ----------------------------
  p_all_domain <- ggplot(bin_meta, aes(x = x_center, y = y_center, fill = cluster_raw)) +
    geom_raster() +
    coord_equal(expand = FALSE) +
    scale_y_reverse() +
    scale_fill_manual(values = cluster2color[clusters], drop = FALSE) +
    labs(
      fill = "Cluster",
      title = paste0(
        "Xenium: ", bin_um, " µm; k = ", k_run,
        out_suffix, " (", label_single, ")"
      )
    ) +
    make_bg_theme()

  ggsave(
    file.path(outdir, paste0("kmeans_grid", bin_um, "um_k", k_run, out_suffix, "_all_domain_colored.png")),
    plot = p_all_domain, width = 6, height = 6, dpi = 300, bg = "black"
  )

  # ----------------------------
  # combined plot of all non-background clusters
  # ----------------------------
  non_bg_clusters <- clusters[cluster2domain[clusters] != "background"]

  if (length(non_bg_clusters) > 0) {

    dt_non_bg <- bin_meta[cluster_raw %in% non_bg_clusters]

    p_non_bg <- ggplot() +
      geom_raster(
        data = bin_meta,
        aes(x = x_center, y = y_center),
        fill = background_color
      ) +
      geom_raster(
        data = dt_non_bg,
        aes(x = x_center, y = y_center, fill = cluster_raw)
      ) +
      coord_equal(expand = FALSE) +
      scale_y_reverse() +
      scale_fill_manual(
        values = cluster2color[non_bg_clusters],
        breaks = non_bg_clusters,
        drop = FALSE
      ) +
      labs(
        fill = "Cluster",
        title = paste0(
          "Xenium: ", bin_um, " µm; k = ", k_run,
          out_suffix, "; all non-background clusters",
          " (", label_single, ")"
        )
      ) +
      make_bg_theme()

    out_png_non_bg <- file.path(
      outdir,
      paste0(
        "kmeans_grid", bin_um, "um_k", k_run,
        out_suffix, "_all_non_background_domain_colored.png"
      )
    )

    ggsave(
      out_png_non_bg,
      plot = p_non_bg, width = 6, height = 6, dpi = 300, bg = "black"
    )

    message("Saved: ", out_png_non_bg)
  } else {
    message("[skip] no non-background clusters found.")
  }

  # ----------------------------
  # pair plots
  # ----------------------------
  if (length(cluster_groups) > 0) {
    for (grp_label in names(cluster_groups)) {

      grp <- cluster_groups[[grp_label]][[1]]
      grp <- intersect(grp, clusters)

      if (length(grp) == 0) {
        message("[skip] ", grp_label, " has no valid clusters in this run.")
        next
      }

      dt_grp <- bin_meta[cluster_raw %in% grp]
      grp_domains <- unique(cluster2domain[grp])
      grp_name <- paste(grp, collapse = "+")

      p_grp <- ggplot() +
        geom_raster(
          data = bin_meta,
          aes(x = x_center, y = y_center),
          fill = background_color
        ) +
        geom_raster(
          data = dt_grp,
          aes(x = x_center, y = y_center, fill = cluster_raw)
        ) +
        coord_equal(expand = FALSE) +
        scale_y_reverse() +
        scale_fill_manual(values = cluster2color[grp], breaks = grp, drop = FALSE) +
        labs(
          fill = "Cluster",
          title = paste0(
            "Xenium: ", bin_um, " µm; k = ", k_run,
            out_suffix, "; ", grp_label, " = ", grp_name,
            "; domain = ", paste(grp_domains, collapse = ", "),
            " (", label_single, ")"
          )
        ) +
        make_bg_theme()

      out_png <- file.path(
        outdir,
        paste0(
          "kmeans_grid", bin_um, "um_k", k_run,
          out_suffix, "_",
          grp_label, "_", grp_name, "_domain_colors.png"
        )
      )

      ggsave(
        out_png,
        plot = p_grp, width = 6, height = 6, dpi = 300, bg = "black"
      )

      message("Saved: ", out_png)
    }
  }

  # ----------------------------
  # single cluster plots
  # ----------------------------
  single_dir <- file.path(outdir, paste0("single_clusters_domain_colored", out_suffix))
  dir.create(single_dir, showWarnings = FALSE, recursive = TRUE)

  for (cl in clusters) {

    dt_cl <- bin_meta[cluster_raw == cl]
    cl_domain <- unname(cluster2domain[cl])
    cl_color  <- unname(cluster2color[cl])

    if (is.na(cl_domain) || length(cl_domain) == 0) cl_domain <- "unannotated"
    if (is.na(cl_color)  || length(cl_color)  == 0) cl_color  <- "black"

    p_single <- ggplot() +
      geom_raster(
        data = bin_meta,
        aes(x = x_center, y = y_center),
        fill = background_color
      ) +
      geom_raster(
        data = dt_cl,
        aes(x = x_center, y = y_center),
        fill = cl_color
      ) +
      coord_equal(expand = FALSE) +
      scale_y_reverse() +
      labs(
        title = paste0(
          "Xenium: ", bin_um, " µm; k = ", k_run,
          out_suffix, "; cluster ", cl,
          "; domain = ", cl_domain,
          " (", label_single, ")"
        )
      ) +
      make_bg_theme()

    out_png <- file.path(
      single_dir,
      paste0(
        "kmeans_grid", bin_um, "um_k", k_run,
        out_suffix, "_cluster", cl, "_", cl_domain, ".png"
      )
    )

    ggsave(
      out_png,
      plot = p_single, width = 6, height = 6, dpi = 300, bg = "black"
    )

    message("Saved single cluster: ", out_png)
  }

  # ----------------------------
  # summary table
  # ----------------------------
  summary_dt <- data.table(
    sample_id = sample_id_run,
    um = grid_name,
    k = k_run,
    stem = stem_tag,
    label = label_single,
    cluster = clusters,
    domain = unname(cluster2domain[clusters]),
    color = unname(cluster2color[clusters])
  )

  fwrite(
    summary_dt,
    file.path(outdir, paste0("cluster_domain_color_map_k", k_run, out_suffix, ".tsv")),
    sep = "\\t"
  )
}

cat("Done.\\n")
EOF

  done
done
