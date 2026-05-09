#!/usr/bin/env Rscript

# Legacy note:
# This file previously mixed Python, shell commands, and an embedded R backup.
# The workflow has been split into task-specific scripts in the same folder:
#
# 1. 01_extract_transcripts_from_parquet.py
#    Read Xenium transcripts.parquet, keep x/y/gene columns, filter genes,
#    and write transcripts_xyz.csv.
#
# 2. xenium_grid_um_from_csv.R
#    Convert transcripts_xyz.csv into sparse bin-by-gene matrices and metadata
#    for any requested bin size.
#
# 3. 02_make_multi_bin_grids.sh
#    Batch-run xenium_grid_um_from_csv.R across multiple bin sizes.
#
# Suggested workflow:
# python 01_extract_transcripts_from_parquet.py --input transcripts.parquet
# bash 02_make_multi_bin_grids.sh transcripts_xyz.csv

message("This is a legacy note file.")
message("Use 01_extract_transcripts_from_parquet.py, xenium_grid_um_from_csv.R, and 02_make_multi_bin_grids.sh instead.")
