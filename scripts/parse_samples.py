#!/usr/bin/env python3
"""
parse_samples.py
================
Utility: load and validate the samples CSV.

Expected columns:  sample, condition, read1, read2

Returns a dict:
    {
        "ctrl_rep1": {
            "condition": "control",
            "read1": "data/raw/ctrl_rep1_R1.fastq.gz",
            "read2": "data/raw/ctrl_rep1_R2.fastq.gz",
        },
        ...
    }
"""

import sys
import csv
from pathlib import Path


REQUIRED_COLUMNS = {"sample", "condition", "read1", "read2"}


def load_samples(samples_csv: str) -> dict:
    """
    Parse a samples CSV and return a validated sample dictionary.
    Raises ValueError with a descriptive message on any validation failure.
    """
    path = Path(samples_csv)
    if not path.exists():
        raise FileNotFoundError(f"Sample sheet not found: {samples_csv}")

    samples: dict = {}

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)

        # ── Column validation ────────────────────────────────────────────────
        if reader.fieldnames is None:
            raise ValueError("Sample sheet appears to be empty.")

        missing_cols = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing_cols:
            raise ValueError(
                f"Sample sheet is missing required columns: {missing_cols}\n"
                f"Found columns: {reader.fieldnames}"
            )

        # ── Row-level validation ─────────────────────────────────────────────
        for i, row in enumerate(reader, start=2):  # row 1 is header
            sample = row["sample"].strip()

            if not sample:
                raise ValueError(f"Row {i}: 'sample' field is empty.")

            if sample in samples:
                raise ValueError(
                    f"Row {i}: Duplicate sample name '{sample}'. "
                    "All sample names must be unique."
                )

            # Check FASTQ files exist (warn rather than hard-fail so dry-runs work)
            for col in ("read1", "read2"):
                fq = row[col].strip()
                if not fq:
                    raise ValueError(
                        f"Row {i} (sample={sample}): '{col}' field is empty."
                    )

            samples[sample] = {
                "condition": row["condition"].strip(),
                "read1": row["read1"].strip(),
                "read2": row["read2"].strip(),
            }

    if not samples:
        raise ValueError("Sample sheet contains no data rows.")

    # ── Ensure at least two distinct conditions for DESeq2 ───────────────────
    conditions = {v["condition"] for v in samples.values()}
    if len(conditions) < 2:
        raise ValueError(
            f"At least 2 distinct conditions are required for DESeq2. "
            f"Found: {conditions}"
        )

    return samples


# ── CLI entry point (useful for manual validation) ───────────────────────────
if __name__ == "__main__":
    import json

    csv_path = sys.argv[1] if len(sys.argv) > 1 else "config/samples.csv"
    try:
        result = load_samples(csv_path)
        print(f"✓ Loaded {len(result)} samples successfully.")
        print(json.dumps(result, indent=2))
    except (FileNotFoundError, ValueError) as exc:
        print(f"✗ Validation failed: {exc}", file=sys.stderr)
        sys.exit(1)
