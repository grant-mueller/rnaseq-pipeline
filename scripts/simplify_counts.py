#!/usr/bin/env python3
"""
simplify_counts.py
==================
Post-process the raw featureCounts output:
  1. Strip the 5 annotation meta-columns (Chr, Start, End, Strand, Length).
  2. Rename BAM-path column headers to clean sample names.
  3. Write a clean TSV: GeneID  sample1  sample2 ...

Usage:
    python simplify_counts.py --input results/counts/counts_matrix.txt \
                               --output results/counts/counts_clean.txt
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


# featureCounts always emits these five annotation columns after GeneID
META_COLUMNS = ["Chr", "Start", "End", "Strand", "Length"]


def simplify_counts(input_path: str, output_path: str) -> None:
    """Strip meta-columns and rename BAM paths to sample names."""

    # ── Read featureCounts output (first row is a comment, skip it) ──────────
    df = pd.read_csv(input_path, sep="\t", comment="#")

    # Verify expected structure
    expected_first_cols = ["Geneid"] + META_COLUMNS
    actual_first_cols = list(df.columns[: len(expected_first_cols)])
    if actual_first_cols != expected_first_cols:
        print(
            f"WARNING: Unexpected column order.\n"
            f"  Expected: {expected_first_cols}\n"
            f"  Got:      {actual_first_cols}",
            file=sys.stderr,
        )

    # ── Drop annotation meta-columns ─────────────────────────────────────────
    cols_to_drop = [c for c in META_COLUMNS if c in df.columns]
    df = df.drop(columns=cols_to_drop)

    # ── Rename BAM path columns → sample names ───────────────────────────────
    # featureCounts column headers look like: results/star/ctrl_rep1/ctrl_rep1.Aligned...bam
    # We extract the sample name as the directory basename.
    rename_map: dict = {}
    for col in df.columns[1:]:  # skip Geneid
        # Try to extract sample name from path like .../star/<sample>/<sample>.*.bam
        match = re.search(r"/star/([^/]+)/", col)
        if match:
            rename_map[col] = match.group(1)
        else:
            # Fall back: use the stem of the filename without extension
            rename_map[col] = Path(col).stem.split(".")[0]

    df = df.rename(columns=rename_map)

    # ── Write output ─────────────────────────────────────────────────────────
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)

    n_genes = len(df)
    n_samples = len(df.columns) - 1  # exclude Geneid
    print(f"✓ Cleaned counts: {n_genes} genes × {n_samples} samples → {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",  required=True, help="Raw featureCounts output")
    parser.add_argument("--output", required=True, help="Cleaned count matrix TSV")
    args = parser.parse_args()
    simplify_counts(args.input, args.output)


if __name__ == "__main__":
    main()
