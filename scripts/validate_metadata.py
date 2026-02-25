#!/usr/bin/env python3
"""
validate_metadata.py
====================
Standalone script to validate:
  1. samples.csv – all required columns, no duplicates, FASTQ files accessible.
  2. config.yaml – required keys are present, reference files exist.

Run before kicking off the pipeline:
    python scripts/validate_metadata.py --config config/config.yaml

Exit code 0 = all checks passed.
Exit code 1 = one or more checks failed.
"""

import argparse
import sys
from pathlib import Path

import yaml

from parse_samples import load_samples


REQUIRED_CONFIG_KEYS = [
    "samples",
    "genome",
    "trimming",
    "star",
    "featurecounts",
    "deseq2",
    "results",
]

REQUIRED_GENOME_KEYS = ["fasta", "gtf", "star_index", "sjdb_overhang"]


def check_file(path: str, label: str, errors: list) -> None:
    """Append an error if path does not exist."""
    if not Path(path).exists():
        errors.append(f"  ✗ {label} not found: {path}")
    else:
        print(f"  ✓ {label}: {path}")


def validate_config(config_path: str) -> list:
    """Return a list of error strings (empty = all good)."""
    errors: list = []

    # ── Load YAML ─────────────────────────────────────────────────────────────
    try:
        with open(config_path) as fh:
            cfg = yaml.safe_load(fh)
    except Exception as exc:
        return [f"Failed to parse config YAML: {exc}"]

    # ── Required top-level keys ───────────────────────────────────────────────
    print("\n[1] Checking config.yaml keys …")
    for key in REQUIRED_CONFIG_KEYS:
        if key not in cfg:
            errors.append(f"  ✗ Missing key: '{key}'")
        else:
            print(f"  ✓ {key}")

    if "genome" in cfg:
        for gkey in REQUIRED_GENOME_KEYS:
            if gkey not in cfg["genome"]:
                errors.append(f"  ✗ Missing genome key: '{gkey}'")

    # ── Check reference files exist ───────────────────────────────────────────
    print("\n[2] Checking reference files …")
    if "genome" in cfg:
        check_file(cfg["genome"]["fasta"], "Genome FASTA", errors)
        check_file(cfg["genome"]["gtf"],   "GTF annotation", errors)
        # STAR index directory is allowed to not exist (it will be built)
        idx = cfg["genome"]["star_index"]
        if not Path(idx).exists():
            print(f"  ⚠ STAR index not found (will be built): {idx}")

    # ── Validate sample sheet ─────────────────────────────────────────────────
    print("\n[3] Validating sample sheet …")
    if "samples" in cfg:
        try:
            samples = load_samples(cfg["samples"])
            print(f"  ✓ {len(samples)} samples loaded successfully.")

            # Check FASTQ files
            missing_fq: list = []
            for sname, meta in samples.items():
                for rkey in ("read1", "read2"):
                    if not Path(meta[rkey]).exists():
                        missing_fq.append(f"    {sname} {rkey}: {meta[rkey]}")
            if missing_fq:
                errors.append(
                    "  ✗ Missing FASTQ files:\n" + "\n".join(missing_fq)
                )
            else:
                print("  ✓ All FASTQ files found.")
        except (FileNotFoundError, ValueError) as exc:
            errors.append(f"  ✗ Sample sheet error: {exc}")

    return errors


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        default="config/config.yaml",
        help="Path to config.yaml (default: config/config.yaml)",
    )
    args = parser.parse_args()

    print(f"Validating pipeline metadata: {args.config}")
    errors = validate_config(args.config)

    if errors:
        print("\n═══ VALIDATION FAILED ═══")
        for err in errors:
            print(err, file=sys.stderr)
        sys.exit(1)
    else:
        print("\n═══ All validation checks passed ✓ ═══")
        sys.exit(0)


if __name__ == "__main__":
    main()
