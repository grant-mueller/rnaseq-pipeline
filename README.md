# RNA-seq Differential Expression Pipeline

[![CI](https://github.com/grant-mueller/rnaseq-pipeline/actions/workflows/ci.yaml/badge.svg)](https://github.com/grant-mueller/rnaseq-pipeline/actions/workflows/ci.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A production-ready, fully containerised RNA-seq differential expression pipeline built with **Snakemake**, **STAR**, **featureCounts**, **DESeq2**, and **MultiQC**.

---

## Overview
```
Raw FASTQ  →  FastQC  →  fastp (trim)  →  FastQC (post-trim)
           →  STAR align  →  featureCounts  →  DESeq2  →  MultiQC
```

| Step | Tool | Version |
|------|------|---------|
| QC | FastQC | 0.12.1 |
| Trimming | fastp | 0.23.4 |
| Alignment | STAR | 2.7.11a |
| Quantification | featureCounts (Subread) | 2.0.6 |
| Differential expression | DESeq2 | 1.42.0 |
| Aggregate QC report | MultiQC | 1.21 |
| Workflow orchestration | Snakemake | 8.5.3 |

---

## Example Dataset

This pipeline was developed and validated using publicly available *C. elegans* RNA-seq data from GEO accession [GSE123054](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123054).

**Study:** Transcriptome analysis of extruded germline from *C. elegans* at different temperatures.

**Comparison:** Empty vector (EV) RNAi-treated worms at 15°C vs 20°C (n=3 replicates per condition).

| Sample | Condition | SRA Accession |
|--------|-----------|---------------|
| EV_15C_rep1 | EV_15C | SRR8257107 |
| EV_15C_rep2 | EV_15C | SRR8257108 |
| EV_15C_rep3 | EV_15C | SRR8257109 |
| EV_20C_rep1 | EV_20C | SRR8257110 |
| EV_20C_rep2 | EV_20C | SRR8257111 |
| EV_20C_rep3 | EV_20C | SRR8257112 |

**Reference genome:** *C. elegans* WBcel235, Ensembl release 111.

**Results:** 164 differentially expressed genes at padj < 0.05, |LFC| > 1 (120 upregulated, 44 downregulated at 15°C relative to 20°C).

---

## Repository Layout
```
rnaseq-pipeline/
├── Snakefile                   # Top-level workflow entry point
├── Dockerfile                  # Full containerised environment
├── config/
│   ├── config.yaml             # Pipeline parameters + paths
│   ├── samples.csv             # Sample sheet (edit this!)
│   └── multiqc_config.yaml     # MultiQC report settings
├── rules/
│   ├── qc.smk                  # FastQC rules
│   ├── trim.smk                # fastp trimming
│   ├── align.smk               # STAR index + alignment
│   ├── quantify.smk            # featureCounts
│   ├── deseq2.smk              # DESeq2 rule
│   └── multiqc.smk             # MultiQC rule
├── scripts/
│   ├── parse_samples.py        # Sample sheet loader + validator
│   ├── simplify_counts.py      # Clean featureCounts output
│   ├── validate_metadata.py    # Pre-flight validation script
│   └── deseq2_analysis.R       # Complete DESeq2 R script
├── envs/
│   ├── qc.yaml                 # FastQC + MultiQC conda env
│   ├── trim.yaml               # fastp conda env
│   ├── align.yaml              # STAR + samtools conda env
│   ├── quantify.yaml           # featureCounts conda env
│   ├── r.yaml                  # R + DESeq2 conda env
│   └── python.yaml             # Python utilities conda env
├── data/
│   ├── raw/                    # Place paired-end FASTQs here
│   └── reference/              # Genome FASTA + GTF go here
├── results/                    # All outputs written here (git-ignored)
└── .github/
    └── workflows/
        └── ci.yaml             # GitHub Actions CI
```

---

## Installation

### Option A – Docker (recommended, zero dependency setup)
```bash
# Build image (~15 min first time)
docker build -t rnaseq-pipeline .

# Run the pipeline
docker run --rm \
  -v $(pwd)/data:/pipeline/data \
  -v $(pwd)/results:/pipeline/results \
  -v $(pwd)/config:/pipeline/config \
  rnaseq-pipeline \
  --configfile config/config.yaml \
  --cores 8
```

### Option B – Conda environments
```bash
# Install Snakemake
conda create -n snakemake -c conda-forge -c bioconda snakemake=8.5.3
conda activate snakemake

# Snakemake manages all other environments automatically
snakemake --use-conda --cores 8
```

---

## Running the Pipeline

### 1. Download public data (optional)

To reproduce the example dataset, install SRA tools and download the reads:
```bash
conda create -n sra-tools -c bioconda sra-tools entrez-direct -y
conda activate sra-tools

# Download all 6 samples
cat accessions.txt | while read acc; do
  fasterq-dump $acc --outdir data/raw --split-files --threads 4 --progress
  gzip data/raw/${acc}_1.fastq data/raw/${acc}_2.fastq
done
```

### 2. Prepare reference genome
```bash
mkdir -p data/reference
cd data/reference

# C. elegans WBcel235 genome and annotation (Ensembl release 111)
wget https://ftp.ensembl.org/pub/release-111/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.111.gtf.gz
gunzip *.gz
```

### 3. Edit the sample sheet

Edit `config/samples.csv` to match your data:
```csv
sample,condition,read1,read2
EV_15C_rep1,EV_15C,data/raw/SRR8257107_1.fastq.gz,data/raw/SRR8257107_2.fastq.gz
EV_15C_rep2,EV_15C,data/raw/SRR8257108_1.fastq.gz,data/raw/SRR8257108_2.fastq.gz
EV_15C_rep3,EV_15C,data/raw/SRR8257109_1.fastq.gz,data/raw/SRR8257109_2.fastq.gz
EV_20C_rep1,EV_20C,data/raw/SRR8257110_1.fastq.gz,data/raw/SRR8257110_2.fastq.gz
EV_20C_rep2,EV_20C,data/raw/SRR8257111_1.fastq.gz,data/raw/SRR8257111_2.fastq.gz
EV_20C_rep3,EV_20C,data/raw/SRR8257112_1.fastq.gz,data/raw/SRR8257112_2.fastq.gz
```

### 4. Edit config.yaml

Review `config/config.yaml` and confirm these key settings:
```yaml
genome:
  fasta: "data/reference/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
  gtf:   "data/reference/Caenorhabditis_elegans.WBcel235.111.gtf"
  sjdb_overhang: 75   # read_length - 1

deseq2:
  reference_level: "EV_20C"   # control condition
```

### 5. Validate metadata
```bash
python scripts/validate_metadata.py --config config/config.yaml
```

### 6. Dry-run to check the DAG
```bash
snakemake --use-conda --cores 8 --dry-run
```

### 7. Execute the pipeline
```bash
# Local execution
snakemake --use-conda --cores 8

# Cluster (SLURM) execution
snakemake --use-conda --cores 200 \
  --executor slurm \
  --default-resources slurm_account=myaccount mem_mb=8000 runtime=120

# Docker execution
docker run --rm \
  -v $(pwd)/data:/pipeline/data \
  -v $(pwd)/results:/pipeline/results \
  rnaseq-pipeline --configfile config/config.yaml --cores 8
```

---

## Expected Outputs
```
results/
├── fastqc/
│   ├── raw/          # Per-sample FastQC HTML + ZIP (raw reads)
│   └── trimmed/      # Per-sample FastQC HTML + ZIP (trimmed reads)
├── trimmed/
│   ├── {sample}_R1.fastq.gz
│   ├── {sample}_R2.fastq.gz
│   └── {sample}_fastp.json / .html
├── star/
│   └── {sample}/
│       ├── {sample}.Aligned.sortedByCoord.out.bam  (.bai index)
│       └── {sample}.Log.final.out
├── counts/
│   ├── counts_matrix.txt          # Raw featureCounts output
│   └── counts_clean.txt           # Cleaned matrix (GeneID x samples)
├── deseq2/
│   ├── DE_results.csv             # Full results table (all genes)
│   ├── normalized_counts.csv      # VST-normalised expression matrix
│   ├── MA_plot.pdf
│   ├── volcano_plot.pdf
│   └── heatmap.pdf                # Top 50 DE genes
└── multiqc/
    └── multiqc_report.html        # Interactive QC summary
```

### Key DESeq2 output columns (`DE_results.csv`)

| Column | Description |
|--------|-------------|
| `gene` | Gene identifier |
| `baseMean` | Mean normalised count across all samples |
| `log2FoldChange` | Shrunk LFC (apeglm) |
| `lfcSE` | Standard error of LFC |
| `pvalue` | Wald test p-value |
| `padj` | Benjamini-Hochberg adjusted p-value |

---

## Adjusting Strandedness

Check strandedness with RSeQC or Salmon before running. Set `featurecounts.strandedness` in `config.yaml`:

| Value | Protocol |
|-------|----------|
| `0` | Unstranded |
| `1` | Stranded (same sense as read) |
| `2` | Reverse-stranded (most dUTP / Illumina TruSeq) |

---

## Reproducibility

All tool versions are pinned in `envs/*.yaml`. The Docker image bakes all dependencies at exact versions. The GitHub Actions CI runs on every push to validate the DAG and rebuild the image.

To record an exact software manifest after a run:
```bash
conda env export -n base > envs/exact_versions.yaml
```

---

## License

MIT – see [LICENSE](LICENSE).
