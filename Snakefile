"""
RNA-seq Differential Expression Pipeline
=========================================
Workflow:
  1.  FastQC   – raw read QC
  2.  fastp    – adapter trimming & quality filtering
  3.  FastQC   – trimmed read QC
  4.  STAR     – genome index (if absent) + alignment
  5.  featureCounts – gene-level quantification
  6.  DESeq2   – differential expression analysis
  7.  MultiQC  – aggregated QC report
"""

import os
from pathlib import Path

# ── Load configuration ──────────────────────────────────────────────────────
configfile: "config/config.yaml"

# ── Parse sample sheet via utility script ───────────────────────────────────
# scripts/parse_samples.py returns a validated dict keyed by sample name
import sys
sys.path.insert(0, "scripts")
from parse_samples import load_samples

SAMPLES = load_samples(config["samples"])
SAMPLE_NAMES = list(SAMPLES.keys())
RESULTS = config["results"]

# ── Helper functions ─────────────────────────────────────────────────────────
def get_raw_read(wildcards, read):
    """Return path to raw R1 or R2 for a given sample wildcard."""
    return SAMPLES[wildcards.sample][read]

def all_bams():
    return expand(
        "{results}/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        results=RESULTS,
        sample=SAMPLE_NAMES,
    )

# ── Include modular rule files ───────────────────────────────────────────────
include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/quantify.smk"
include: "rules/deseq2.smk"
include: "rules/multiqc.smk"

# ── Default target ───────────────────────────────────────────────────────────
rule all:
    input:
        # Raw FastQC
        expand(
            "{results}/fastqc/raw/{sample}_{read}_fastqc.html",
            results=RESULTS,
            sample=SAMPLE_NAMES,
            read=["R1", "R2"],
        ),
        # Trimmed FastQC
        expand(
            "{results}/fastqc/trimmed/{sample}_{read}_fastqc.html",
            results=RESULTS,
            sample=SAMPLE_NAMES,
            read=["R1", "R2"],
        ),
        # Alignments
        all_bams(),
        # featureCounts matrix
        f"{RESULTS}/counts/counts_matrix.txt",
        # DESeq2 outputs
        f"{RESULTS}/deseq2/DE_results.csv",
        f"{RESULTS}/deseq2/normalized_counts.csv",
        f"{RESULTS}/deseq2/MA_plot.pdf",
        f"{RESULTS}/deseq2/volcano_plot.pdf",
        f"{RESULTS}/deseq2/heatmap.pdf",
        # MultiQC
        f"{RESULTS}/multiqc/multiqc_report.html",
