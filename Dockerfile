# =============================================================================
# Dockerfile – RNA-seq Differential Expression Pipeline
# =============================================================================
# Base: Ubuntu 22.04 LTS
# Installs: Snakemake, STAR, featureCounts (Subread), fastp, FastQC,
#           MultiQC, samtools, R + DESeq2 + tidyverse + plotting packages
# =============================================================================

FROM ubuntu:22.04

LABEL maintainer="bioinformatics@example.com" \
      description="RNA-seq DE pipeline – STAR | featureCounts | DESeq2"

# ── Non-interactive frontend to avoid apt prompts ────────────────────────────
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# ── System dependencies ───────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    bzip2 \
    ca-certificates \
    git \
    unzip \
    default-jre \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgfortran5 \
    locales \
    && rm -rf /var/lib/apt/lists/*

RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8

# ── Miniforge (conda-forge by default, no ToS issues) ────────────────────────
ENV CONDA_DIR=/opt/conda
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
        -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniforge.sh

ENV PATH=${CONDA_DIR}/bin:$PATH

# ── Conda channels ────────────────────────────────────────────────────────────
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# ── Install bioinformatics tools via conda ────────────────────────────────────
# Split into logical groups to leverage Docker layer caching
RUN conda install -y -n base \
    "python>=3.11,<3.13" \
    "snakemake>=8.0,<9.0" \
    && conda clean -afy

RUN conda install -y -n base \
    "star>=2.7" \
    "samtools>=1.19" \
    && conda clean -afy

RUN conda install -y -n base \
    "subread>=2.0" \
    "fastp>=0.23" \
    "fastqc>=0.12" \
    "multiqc>=1.21" \
    && conda clean -afy

RUN conda install -y -n base \
    "pandas>=2.0" \
    "pyyaml>=6.0" \
    && conda clean -afy

# ── R + Bioconductor packages ─────────────────────────────────────────────────
RUN conda install -y -n base \
    "r-base>=4.3" \
    "bioconductor-deseq2>=1.42" \
    "bioconductor-apeglm>=1.24" \
    "r-ggplot2>=3.4" \
    "r-ggrepel>=0.9" \
    "r-pheatmap>=1.0" \
    "r-rcolorbrewer>=1.1" \
    "r-dplyr>=1.1" \
    "r-tidyr>=1.3" \
    "r-optparse>=1.7" \
    "r-tibble>=3.2" \
    && conda clean -afy

# ── Set working directory ─────────────────────────────────────────────────────
WORKDIR /pipeline

# ── Copy pipeline files ───────────────────────────────────────────────────────
COPY . /pipeline/

# ── Verify key tool versions at build time ────────────────────────────────────
RUN echo "=== Tool versions ===" && \
    snakemake --version && \
    STAR --version && \
    featureCounts -v 2>&1 | head -2 && \
    fastp --version 2>&1 | head -1 && \
    fastqc --version && \
    multiqc --version && \
    samtools --version | head -1 && \
    Rscript -e 'packageVersion("DESeq2")' && \
    echo "=== All tools available ==="

ENTRYPOINT ["snakemake"]
CMD ["--help"]
