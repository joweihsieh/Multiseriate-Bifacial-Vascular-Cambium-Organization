#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="/home/woodydrylab/FileShare/20260121_Xenium"

find "$BASE_DIR" \
  \( -path "*/ray_match3/*_real_vs_permB_topK_metrics.tsv" -o -path "*/ray_match_final/ray_match3/*_real_vs_permB_topK_metrics.tsv" \) \
  -type f | while read -r INPUT_TSV
do
  echo "Processing: $INPUT_TSV"

  Rscript - "$INPUT_TSV" <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
INPUT_TSV <- args[1]

if (!file.exists(INPUT_TSV)) {
  stop("File not found: ", INPUT_TSV)
}

topk_all_dt <- fread(INPUT_TSV)

required_cols <- c("score", "perm_group")
miss <- setdiff(required_cols, names(topk_all_dt))
if (length(miss) > 0) {
  stop("Missing required columns in ", INPUT_TSV, ": ", paste(miss, collapse = ", "))
}

top_k <- suppressWarnings(max(topk_all_dt$topk_rank, na.rm = TRUE))
if (!is.finite(top_k)) top_k <- NA_integer_

out_prefix <- sub("_real_vs_permB_topK_metrics\\.tsv$", "", INPUT_TSV)


p_topk_score <- ggplot(topk_all_dt, aes(x = score, color = perm_group, fill = perm_group)) +
  geom_density(alpha = 0.25, size = 1) +

  scale_color_manual(values = c(
    "real" = "limegreen",     
    "permB" = "grey"     
  )) +
  scale_fill_manual(values = c(
    "real" = "limegreen",
    "permB" = "grey"
  )) +

  theme_bw(base_size = 12) +
  labs(
    title = if (is.finite(top_k)) {
      paste0("Real vs permutation-B top ", top_k, ": best_score distribution")
    } else {
      "Real vs permutation-B: best_score distribution"
    },
    x = "best_score",
    y = "Density"
  ) +
  theme(
    plot.background   = element_rect(fill = "transparent", color = NA),
    panel.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key        = element_rect(fill = "transparent", color = NA)
  )

pdf_file <- paste0(out_prefix, "_score_density_transparent.pdf")
png_file <- paste0(out_prefix, "_score_density_transparent.png")

ggsave(pdf_file, p_topk_score, width = 7, height = 5, bg = "transparent")
ggsave(png_file, p_topk_score, width = 7, height = 5, dpi = 300, bg = "transparent")

message("Saved: ", pdf_file)
message("Saved: ", png_file)
EOF

done

echo "Done."