# Xenium Preprocessing

This folder is organized as a small two-step workflow:

1. Extract transcript coordinates from `transcripts.parquet`
2. Build bin-by-gene matrices for one or more bin sizes

## Files

- `01_extract_transcripts_from_parquet.py`
  Extract `x`, `y`, and `gene` columns from Xenium `transcripts.parquet`, filter genes by prefix, and write `transcripts_xyz.csv`.
- `xenium_grid_um_from_csv.R`
  Convert `transcripts_xyz.csv` into:
  - `counts_bins_by_genes_sparse.rds`
  - `bin_metadata.tsv`
  - `genes.tsv`
- `02_make_multi_bin_grids.sh`
  Batch-run the R script across multiple bin sizes.
- `1. Transcript_5um.R`
  Legacy note pointing to the new workflow.

## Recommended workflow

### 1. Export transcripts from parquet

```bash
python 01_extract_transcripts_from_parquet.py \
  --input transcripts.parquet \
  --output transcripts_xyz.csv \
  --gene-prefix Potri_
```

### 2. Build multiple grid sizes

```bash
bash 02_make_multi_bin_grids.sh transcripts_xyz.csv
```

This generates:

- `grid1um_out`
- `grid2um_out`
- `grid3um_out`
- `grid4um_out`
- `grid5um_out`
- `grid6um_out`
- `grid7um_out`
- `grid8um_out`
- `grid9um_out`
- `grid10um_out`

If you want a different neighbor threshold:

```bash
bash 02_make_multi_bin_grids.sh transcripts_xyz.csv 2
```

## Single bin size example

If you only want one grid size, run the R script directly:

```bash
Rscript xenium_grid_um_from_csv.R \
  --input transcripts_xyz.csv \
  --outdir grid5um_out \
  --bin_um 5 \
  --neighbor_k 2
```
