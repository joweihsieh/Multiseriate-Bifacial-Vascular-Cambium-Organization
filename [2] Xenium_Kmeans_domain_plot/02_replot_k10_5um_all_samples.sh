#!/usr/bin/env bash
# Re-plot existing 5 um / k = 10 outputs for all Xenium samples.
# Make sure you are inside the Xenium output folder first
# conda activate r_spatial
# cd /home/woodydrylab/FileShare/20260121_Xenium/output-XXXX...


#################
# Run all output folders
for d in /home/woodydrylab/FileShare/20260121_Xenium/output-*
do
  cd "$d"
  bash /home/woodydrylab/FileShare/20260121_Xenium/03_plot_k10_domain_recolored.sh
done

################# download
cd /home/woodydrylab/FileShare/20260121_Xenium

tar -czf xenium_selected_plots.tar.gz \
  $(find output-* \
    -path "*/kmeans_k10_raw_out/single_clusters_domain_colored*/*.png" \
    -o -path "*/kmeans_k10_raw_out/*_domain_colored.png")



tar -xzf xenium_selected_plots.tar.gz
