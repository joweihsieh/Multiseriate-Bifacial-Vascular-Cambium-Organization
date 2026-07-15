## Instructions to run on data

### 1. Obtain and organize the input data

This repository contains analysis scripts but does not define a single bundled input-data directory. The large Xenium, LCM RNA-seq, and single-nucleus/single-cell RNA-seq datasets should be obtained separately and placed on local or shared storage.

A practical directory layout is:

```text
project_data/
├── xenium/
│   ├── output-sample-01/
│   │   ├── transcripts.parquet
│   │   ├── Selection_Bottom_coordinates.csv
│   │   └── Selection_Upper_coordinates.csv
│   ├── output-sample-02/
│   └── ...
├── lcm/
│   ├── Ptr_*_gene_abundances.csv
│   └── ...
├── snRNA_or_scRNA/
│   ├── plotting_TenX_Ptr.csv
│   ├── geneUMI_TenX_Ptr.csv
│   └── Cell Ranger or Seurat outputs
├── annotations/
│   ├── k10_5domain.xlsx
│   └── gene_annotation.tsv
└── results/
```

This layout is only a recommendation. The current scripts do not automatically discover this structure, so the paths in the scripts must be updated before execution.


### 2. Xenium transcript preprocessing

The first reusable command-line workflow is in `[1.1] Xenium_Preprocessing`.

Export transcript coordinates from a Xenium `transcripts.parquet` file:

```bash
cd "[1.1] Xenium_Preprocessing"

python 01_extract_transcripts_from_parquet.py \
  --input "/path/to/xenium/output-sample/transcripts.parquet" \
  --output "/path/to/xenium/output-sample/transcripts_xyz.csv" \
  --gene-prefix Potri_
```

Generate grid-based count matrices for the configured bin sizes:

```bash
bash 02_make_multi_bin_grids.sh \
  "/path/to/xenium/output-sample/transcripts_xyz.csv"
```

To analyze only one bin size, run the R script directly:

```bash
Rscript xenium_grid_um_from_csv.R \
  --input "/path/to/xenium/output-sample/transcripts_xyz.csv" \
  --outdir "/path/to/xenium/output-sample/grid5um_out" \
  --bin_um 5 \
  --neighbor_k 2
```

Repeat this step for each Xenium sample.

### 3. Bin-level quality control

After the grid matrices have been generated, edit `XENIUM_BASE` in the scripts under `[1.2] Xenium_bin_Transcript_Counts`, and run:

```bash
Rscript "[1.2] Xenium_bin_Transcript_Counts/01_plot_5um_bin_histograms.R"
Rscript "[1.2] Xenium_bin_Transcript_Counts/02_make_5um_qc_summary.R"
Rscript "[1.2] Xenium_bin_Transcript_Counts/03_make_50um_qc_summary.R"
Rscript "[1.2] Xenium_bin_Transcript_Counts/04_plot_50um_gmean_boxplot.R"
```

These scripts expect the relevant grid directory and, where multiple tissue sections occur in the same Xenium output, the upper and lower coordinate-selection CSV files.

### 4. K-means clustering and domain plots

The K-means workflow is under `[2] Xenium_Kmeans_domain_plot`.

Run the clustering script from inside one Xenium sample output directory:

```bash
cd "/path/to/xenium/output-sample"

bash "/path/to/Multiseriate-Bifacial-Vascular-Cambium-Organization/[2] Xenium_Kmeans_domain_plot/01_run_kmeans_grid_sweep.sh"
```

In the currently committed shell script, the active loop is set to the 5-µm grid. To run all grid sizes, change:

```bash
for um in 5
```

to:

```bash
for um in 1 2 3 4 5 6 7 8 9 10
```

The script runs K-means for `k = 2` through `k = 10`. The reported 5-µm multi-cluster analysis generally required less than 5 hours, although runtime depends on sample size and hardware.

For recoloring and domain plotting, update the absolute paths and annotation workbook in the shell scripts, then run the relevant scripts in order:

```bash
bash "[2] Xenium_Kmeans_domain_plot/02_replot_k10_5um_all_samples.sh"
bash "[2] Xenium_Kmeans_domain_plot/03_plot_k10_domain_recolored.sh"
bash "[2] Xenium_Kmeans_domain_plot/04_plot_few_clusters_recolor.sh"
```

### 5. Downstream workflows

The remaining numbered folders are downstream or complementary analyses. After editing their input and output paths, run the scripts in numerical order within the selected folder.

| Folder | Main purpose | General execution pattern |
|---|---|---|
| `[0] Xenium_Probe` | Probe annotation and LCM expression/rank summaries | `Rscript` for scripts `00`–`02` |
| `[3] Xenium_LCM_RNAseq_correlation` | Xenium–LCM expression correlation and spatial maps | `Rscript` for scripts `01`–`02` |
| `[4.1] Xenium_multiple_section_separation` | Separate multiple tissue sections within shared Xenium outputs | Run the scripts inside the relevant sample-pair subfolder |
| `[4.2] Xenium_cambium_single` | Cambium-only clustering and assessment of ray/fusiform separation | Run shell script `01`, then R script `02` |
| `[4.3] Xeniume_cambium_isolation_k2_plot` | Further K=2 clustering of isolated cambium and cluster-composition summaries | Run shell script `01`, then R script `02` |
| `[4.4] Xenium_Domain_composition_k10_back_to_K2-9` | Trace K=10 domain annotations back to K=2–9 solutions | Run R script `01` |
| `[5.1] Ray_lineage` | Ray-lineage scoring, top-ranked features, and permutation plots | Run R scripts `01`–`02`, then shell script `03` |
| `[5.2] Xenium_scRNA_ray_fusiform_lineage` | Xenium–scRNA organizer distances and statistical comparisons | Run R script `01`, then R script `02` |
| `[6.1] Simulation` | Simulated uniseriate, multiseriate, and lineage-segregated models | Run R scripts `01`–`04` |
| `[6.2] Fusiform_cambium_k2to4` | Fusiform-cambium subset clustering and inner/outer distributions | Run scripts `01`–`04` in order |
| `[7] Xenium_Cosine_Similarity` | Spatial cosine-similarity analyses for each ROI | Run R scripts `01`–`03` |
| `[8] Xenium_DEG` | Domain-level differential expression, heatmaps, annotation, and GO analysis | Run R scripts `01`–`02` |
| `[9] snRNA_overlapping_rate` | Overlap-rate analysis between single-nucleus datasets | Run R script `01` |
| `[10] Xenium_snRNA_correlation_GMM_DEG` | snRNA annotation, Xenium–snRNA correlation, GMM assignment, DEG, and heatmaps | Run R scripts `01`–`05` in order |

These modules do not all need to be run for every reproduction task. Select the modules corresponding to the figure or analysis being reproduced.

### 6. Xenium–scRNA distance analysis example

Edit the following variables in `[5.2] Xenium_scRNA_ray_fusiform_lineage/01_Xenium_scRNA_distance.R`:

```text
XENIUM_BASE
DOMAIN_XLSX
PLOTTING_FILE
UMI_FILE
OUTDIR
```

Then run:

```bash
Rscript "[5.2] Xenium_scRNA_ray_fusiform_lineage/01_Xenium_scRNA_distance.R"
```

Next, set `INFILE` in `02_Xenium_scRNA_distance_test_only.R` to the `observed_metrics.tsv` generated by the first script, set a writable `OUTDIR`, and run:

```bash
Rscript "[5.2] Xenium_scRNA_ray_fusiform_lineage/02_Xenium_scRNA_distance_test_only.R"
```

The checked-in statistical script first removes values outside `Q1 - 1.5 × IQR` and `Q3 + 1.5 × IQR` independently within each comparison group, and then performs unpaired pairwise Welch t-tests with `p.adjust.method = "none"`. The generated t-test table therefore contains unadjusted two-sided P values. If one-tailed P values are reported in the manuscript, the direction check and conversion should be performed in code and written to a separate output table rather than applied manually.

### 7. Minimal installation test without external data

The main simulation script can be used to confirm that R and the essential plotting packages work without requiring the primary datasets:

```bash
Rscript "[6.1] Simulation/01_simulation_3_models_theory_final.R"
```

A successful run creates:

```text
theoretical_patterns_final_cambium_xylem/
├── 00_base_patterns.png
└── 01_k2_to_k6_compare.png
```

The terminal should finish with a message indicating that the PNG outputs were written to the simulation output directory.

---

## Expected output

### General output behavior

Output locations are controlled by variables such as `OUTDIR`, by command-line `--outdir` arguments, or by the current working directory. There is not yet one standardized repository-wide `results/` directory. A successful run usually prints `Done.` or reports the output directory in the terminal.

The workflows generate four main output types:

```text
.tsv / .csv   tabular results and metadata
.rds          R matrices or serialized analysis objects
.png / .pdf   figures
.log / .txt   terminal logs, summaries, or environment information
```

### Expected outputs by workflow

| Folder | Principal expected outputs |
|---|---|
| `[0] Xenium_Probe` | Annotated probe/gene tables, LCM expression-tertile summaries, rank-based summaries, and ternary plots. |
| `[1.1] Xenium_Preprocessing` | `transcripts_xyz.csv`; for each grid size, a `grid*um_out/` directory containing `counts_bins_by_genes_sparse.rds`, `bin_metadata.tsv`, and `genes.tsv`. |
| `[1.2] Xenium_bin_Transcript_Counts` | Per-sample transcript-count histograms; `Xenium_full_QC_summary.tsv`; `Xenium_full_QC_summary_50um.tsv`; 50-µm histograms and summary boxplots. |
| `[2] Xenium_Kmeans_domain_plot` | `grid*um_out/kmeans_k*_raw_out/` directories; `bin_metadata_with_cluster_raw.tsv`; all-cluster and single-cluster spatial PNGs; recolored domain plots; optionally `xenium_selected_plots.tar.gz`. |
| `[3] Xenium_LCM_RNAseq_correlation` | `bin_metadata_LCM.tsv`; correlation maps in `corr_maps_spearman_pearson/`; thresholded/smoothed correlation maps in `corr_maps_byLCMthreshold_with_mid_xslice_split/`. |
| `[4.1] Xenium_multiple_section_separation` | Section-specific metadata/count files and spatial plots for the separated tissue sections. |
| `[4.2] Xenium_cambium_single` | Cambium-only clustering outputs, spatial cluster plots, and summaries describing the K value at which ray and fusiform groups separate. |
| `[4.3] Xeniume_cambium_isolation_k2_plot` | K=2 cambium-subset cluster assignments, spatial plots, and cluster count/proportion summaries. |
| `[4.4] Xenium_Domain_composition_k10_back_to_K2-9` | Tables linking K=10 domain annotations to K=2–9 clusters, plus composition summaries and plots. |
| `[5.1] Ray_lineage` | Ray/fusiform lineage-score tables, top-ranked features or genes, permutation/background results, density plots, and overlay figures. |
| `[5.2] Xenium_scRNA_ray_fusiform_lineage` | Xenium initial-cell assignments, group bin counts, organizer and Xenium non-zero-mean expression tables, filtering summaries, observed/permuted distance tables, boxplots, an RDS result object, outlier summaries, and pairwise test tables. |
| `[6.1] Simulation` | Simulated spatial-pattern and K-means comparison PNGs, including `00_base_patterns.png` and `01_k2_to_k6_compare.png` from the main three-model script. |
| `[6.2] Fusiform_cambium_k2to4` | Fusiform-cambium subset K-means directories, recolored spatial plots, and inner/outer histogram or density plots for selected K values. |
| `[7] Xenium_Cosine_Similarity` | Per-ROI cosine-similarity tables and spatial plots, plus combined summary boxplots. |
| `[8] Xenium_DEG` | Differential-expression tables, five-domain heatmaps, gene-annotation tables, and GO-enrichment or GO-annotation outputs. |
| `[9] snRNA_overlapping_rate` | Pairwise overlap counts/rates between single-nucleus datasets and associated summary tables or plots. |
| `[10] Xenium_snRNA_correlation_GMM_DEG` | Annotated snRNA objects/tables, Xenium–snRNA correlation matrices and plots, GMM assignments, DEG tables, and final heatmaps. |

