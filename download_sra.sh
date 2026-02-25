#!/bin/bash
# =============================================================================
# download_sra.sh
# Downloads and compresses paired-end FASTQ files from SRA accessions.
# Usage: bash download_sra.sh
# Reads accession numbers one per line from accessions.txt
# =============================================================================

set -euo pipefail  # exit on error, undefined variable, or pipe failure

ACCESSIONS_FILE="accessions.txt"
OUTDIR="data/raw"

# Check accessions.txt exists before starting
if [ ! -f "$ACCESSIONS_FILE" ]; then
    echo "Error: $ACCESSIONS_FILE not found. Please create it first."
    exit 1
fi

mkdir -p "$OUTDIR"

while read accession; do
    echo "Downloading $accession ..."
    fasterq-dump "$accession" \
        --outdir "$OUTDIR" \
        --split-files \
        --threads 8 \
        --progress

    echo "Compressing $accession ..."
    gzip "${OUTDIR}/${accession}_1.fastq"
    gzip "${OUTDIR}/${accession}_2.fastq"

    echo "Done: $accession"
    echo "---"
done < "$ACCESSIONS_FILE"

echo "All downloads complete."
