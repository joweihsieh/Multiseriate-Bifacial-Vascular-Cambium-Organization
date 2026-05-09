# Xenium K-means

This folder contains two different stages of the Xenium K-means workflow:

1. Run K-means on grid-based expression matrices across multiple bin sizes and `k` values
2. Re-plot the existing 5 um / `k = 10` results with manual domain annotations and fixed colors

## Recommended entry points

- `01_run_kmeans_grid_sweep.sh`
  Run K-means for grid sizes 1-10 um and `k = 2-10`.
- `02_replot_k10_5um_all_samples.sh`
  Re-plot existing `5 um / k = 10` outputs for all samples.
- `03_plot_k10_domain_recolored.sh`
  Re-plot one sample's existing `5 um / k = 10` output using domain annotations from Excel.

## Previous filenames

These were the original filenames before renaming:

- `2. Kmeans_full_1umto10um_k2to10.sh`
- `2. Kmeans_full_5um_k10_twomore_clusters_assignedcolors_fivedomain_all_20260405.sh`
- `plot_few_clusters_recolor_20260405.sh`

## Workflow summary

### 1. Run the full K-means sweep

```bash
bash 01_run_kmeans_grid_sweep.sh
```

This expects to be run inside a Xenium sample output folder such as:

```bash
cd /home/woodydrylab/FileShare/20260121_Xenium/output-XXXX...
```

### 2. Re-plot annotated 5 um / k = 10 results for all samples

```bash
bash 02_replot_k10_5um_all_samples.sh
```

This loops over all `output-*` folders and calls the single-sample re-plot script.

### 3. Re-plot one sample only

```bash
cd /home/woodydrylab/FileShare/20260121_Xenium/output-XXXX...
bash /home/woodydrylab/FileShare/20260121_Xenium/03_plot_k10_domain_recolored.sh
```
