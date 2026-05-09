#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def pick_column(columns, candidates, label):
    for name in candidates:
        if name in columns:
            return name
    raise ValueError(
        f"Could not find a {label} column. Tried: {', '.join(candidates)}"
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract x/y/gene columns from Xenium transcripts.parquet and write "
            "a simplified transcripts_xyz.csv file."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to transcripts.parquet",
    )
    parser.add_argument(
        "--output",
        default="transcripts_xyz.csv",
        help="Output CSV path (default: transcripts_xyz.csv)",
    )
    parser.add_argument(
        "--gene-prefix",
        default="Potri_",
        help="Keep only genes starting with this prefix (default: Potri_)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    df = pd.read_parquet(input_path)

    x_col = pick_column(df.columns, ["x_location", "x"], "x")
    y_col = pick_column(df.columns, ["y_location", "y"], "y")
    g_col = pick_column(df.columns, ["feature_name", "gene"], "gene")

    df = df[[x_col, y_col, g_col]].rename(
        columns={x_col: "x", y_col: "y", g_col: "gene"}
    )

    df = df[df["gene"].astype(str).str.startswith(args.gene_prefix, na=False)]
    df.to_csv(output_path, index=False)

    print(
        f"Wrote {output_path} with columns {list(df.columns)} and {len(df)} rows."
    )


if __name__ == "__main__":
    main()
