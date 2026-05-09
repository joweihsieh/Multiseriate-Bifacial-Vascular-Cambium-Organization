# LCM Correlation Scripts

This folder contains two related scripts:

## 1. `01_correlation_allsamples.R`

Purpose:
- Compute LCM-vs-Xenium correlations for all available `output-*` samples and `grid*um_out` folders.
- Write correlation columns into `bin_metadata_LCM.tsv`.
- Generate one PNG per correlation column under `corr_maps_spearman_pearson/`.

Main inputs:
- Xenium bin matrices: `counts_bins_by_genes_sparse.rds`
- Xenium bin metadata: `bin_metadata.tsv`
- LCM reference tables: the three `Ptr_*_gene_abundances.csv` files

Main outputs:
- `bin_metadata_LCM.tsv`
- `corr_maps_spearman_pearson/*.png`

## 2. `02_correlation_plot_only_with_a_line_smooth_cyan_clean.R`

Purpose:
- Plot precomputed `corr_LCM_*` columns from `bin_metadata_LCM.tsv`.
- Show the heatmap plus a smoothed middle x-slice line.
- For Cambium plots, optionally draw dashed vertical guide lines from the lower curve to the slice line.

Main inputs:
- `bin_metadata_LCM.tsv`
- `Selection_Bottom_coordinates.csv`
- `Selection_Upper_coordinates.csv`

Main outputs:
- PNG figures under `corr_maps_byLCMthreshold_with_mid_xslice_split/`

## Notes

- I kept the original filenames unchanged.
- Only the plot-only script needed Chinese comment cleanup; the all-samples script was already in English.
- Both scripts parse successfully after the comment updates.

