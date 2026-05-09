# Make sure you are inside the Xenium output folder first
# conda activate r_spatial
# cd /home/woodydrylab/FileShare/20260121_Xenium/output-XXXX...
# Running a single sample takes more than 5 hours.

#for um in 1 2 3 4 5 6 7 8 9 10
for um in 5
do
  grid_dir="grid${um}um_out"

  # Skip if the grid_dir does not exist to avoid errors
  if [ ! -d "${grid_dir}" ]; then
    echo "[skip] ${grid_dir} not found"
    continue
  fi

  echo "=============================="
  echo "Grid = ${um} um  (${grid_dir})"
  echo "=============================="

  for k in 2 3 4 5 6 7 8 9 10
  do
    echo "  Running k = $k"
    outdir="${grid_dir}/kmeans_k${k}_raw_out"

    Rscript /home/woodydrylab/FileShare/20260121_Xenium/xenium_um_grid_kmeans_raw.R \
      --matrix ${grid_dir}/counts_bins_by_genes_sparse.rds \
      --binmeta ${grid_dir}/bin_metadata.tsv \
      --outdir ${outdir} \
      --k ${k} \
      --n_pcs 30

    Rscript - <<EOF
library(data.table)
library(ggplot2)
library(scales)

outdir <- "${outdir}"
k <- ${k}
bin_um <- ${um}

bin_meta <- fread(file.path(outdir, "bin_metadata_with_cluster_raw.tsv"))
bin_meta[, cluster_raw := factor(cluster_raw)]

clusters <- levels(bin_meta[["cluster_raw"]])
pal <- hue_pal()(length(clusters))
names(pal) <- clusters

# all clusters
p_all <- ggplot(bin_meta, aes(x = x_center, y = y_center, fill = cluster_raw)) +
  geom_tile(width = bin_um, height = bin_um) +
  coord_equal() +
  scale_y_reverse() +
  scale_fill_manual(values = pal, drop = FALSE) +
  labs(fill = "Cluster", title = paste0("Xenium: ", bin_um, " µm; k = ", k)) +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.key = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    plot.title = element_text(color = "white", hjust = 0.5)
  )

ggsave(file.path(outdir, paste0("kmeans_grid", bin_um, "um_k", k, "_all.png")),
       plot = p_all, width = 6, height = 6, dpi = 300, bg = "black")

# each cluster separately (same color as p_all)
for (cl in clusters) {
  dt_cl <- bin_meta[cluster_raw == cl]
  cl_color <- pal[[cl]]

  p_one <- ggplot() +
    geom_tile(data = bin_meta, aes(x = x_center, y = y_center),
              width = bin_um, height = bin_um, fill = "black") +
    geom_tile(data = dt_cl, aes(x = x_center, y = y_center),
              width = bin_um, height = bin_um, fill = cl_color) +
    coord_equal() +
    scale_y_reverse() +
    labs(title = paste0("Xenium: ", bin_um, " µm; k = ", k, "; cluster = ", cl)) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      plot.title = element_text(color = "white", hjust = 0.5)
    )

  ggsave(file.path(outdir, paste0("kmeans_grid", bin_um, "um_k", k, "_cluster_", cl, ".png")),
         plot = p_one, width = 6, height = 6, dpi = 300, bg = "black")
}
EOF

  done
done
