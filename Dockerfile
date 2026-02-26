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

# ── Miniconda ─────────────────────────────────────────────────────────────────
ENV CONDA_DIR=/opt/conda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniconda.sh

ENV PATH=${CONDA_DIR}/bin:$PATH

# ── Conda channels ────────────────────────────────────────────────────────────
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# ── Install bioinformatics tools via conda ────────────────────────────────────
# Split into logical groups to leverage Docker layer caching
RUN conda install -y -n base -c conda-forge -c bioconda \
    "snakemake>=8.0,<9.0" \
    && conda clean -afy

RUN conda install -y -n base \
    star=2.7.11a \
    samtools=1.19 \
    && conda clean -afy

RUN conda install -y -n base \
    subread=2.0.6 \
    fastp=0.23.4 \
    fastqc=0.12.1 \
    multiqc=1.21 \
    && conda clean -afy

RUN conda install -y -n base \
    pandas=2.1.4 \
    pyyaml=6.0.1 \
    && conda clean -afy

# ── R + Bioconductor packages ─────────────────────────────────────────────────
RUN conda install -y -n base \
    r-base=4.3.2 \
    bioconductor-deseq2=1.42.0 \
    bioconductor-apeglm=1.24.0 \
    r-ggplot2=3.4.4 \
    r-ggrepel=0.9.4 \
    r-pheatmap=1.0.12 \
    r-rcolorbrewer=1.1_3 \
    r-dplyr=1.1.4 \
    r-tidyr=1.3.0 \
    r-optparse=1.7.3 \
    r-tibble=3.2.1 \
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
