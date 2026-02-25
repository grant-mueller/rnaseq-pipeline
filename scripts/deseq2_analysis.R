#!/usr/bin/env Rscript
# =============================================================================
# deseq2_analysis.R
# =============================================================================
# Complete DESeq2 differential expression analysis script.
#
# Inputs:
#   --counts        Path to clean featureCounts TSV (GeneID × samples)
#   --samples       Path to samples.csv (sample, condition, read1, read2)
#   --outdir        Output directory for all results
#   --condition_col Column in samples.csv used as the grouping factor
#   --reference     Reference / control condition level
#   --alpha         Adjusted p-value significance cutoff (default 0.05)
#   --lfc_threshold Log2 fold-change threshold for volcano/MA (default 1)
#
# Outputs:
#   DE_results.csv        – full DESeq2 results table (all genes)
#   normalized_counts.csv – VST-normalised count matrix
#   MA_plot.pdf           – MA plot (log fold-change vs mean expression)
#   volcano_plot.pdf      – Volcano plot (−log10 p-adj vs log2FC)
#   heatmap.pdf           – Heatmap of top 50 DE genes (VST counts)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
})

# ── Argument parsing ──────────────────────────────────────────────────────────
option_list <- list(
  make_option("--counts",        type = "character", help = "Clean counts TSV"),
  make_option("--samples",       type = "character", help = "Samples CSV"),
  make_option("--outdir",        type = "character", help = "Output directory"),
  make_option("--condition_col", type = "character", default = "condition"),
  make_option("--reference",     type = "character", default = "control"),
  make_option("--alpha",         type = "double",    default = 0.05),
  make_option("--lfc_threshold", type = "double",    default = 1.0)
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required args
for (arg in c("counts", "samples", "outdir")) {
  if (is.null(opt[[arg]])) {
    stop(paste("Missing required argument: --", arg, sep = ""))
  }
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

message("═══ DESeq2 Analysis ═══")
message(sprintf("Counts:        %s", opt$counts))
message(sprintf("Samples:       %s", opt$samples))
message(sprintf("Output dir:    %s", opt$outdir))
message(sprintf("Condition col: %s", opt$condition_col))
message(sprintf("Reference:     %s", opt$reference))
message(sprintf("Alpha:         %s", opt$alpha))
message(sprintf("LFC threshold: %s", opt$lfc_threshold))

# ── Load data ─────────────────────────────────────────────────────────────────

# Count matrix: rows = genes, cols = samples
message("\n[1] Loading count matrix …")
counts_df <- read.table(opt$counts, header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE)
message(sprintf("    %d genes × %d samples", nrow(counts_df), ncol(counts_df)))

# Sample metadata
message("[2] Loading sample metadata …")
col_data <- read.csv(opt$samples, row.names = 1, check.names = FALSE)
col_data <- col_data[colnames(counts_df), , drop = FALSE]  # align row order

# Coerce condition column to factor with specified reference level first
col_data[[opt$condition_col]] <- factor(
  col_data[[opt$condition_col]],
  levels = c(opt$reference,
             setdiff(unique(col_data[[opt$condition_col]]), opt$reference))
)
message(sprintf("    Conditions: %s",
                paste(levels(col_data[[opt$condition_col]]), collapse = " vs ")))

# ── Build DESeqDataSet ────────────────────────────────────────────────────────
message("\n[3] Building DESeqDataSet …")
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_df),   # ensure integer counts
  colData   = col_data,
  design    = as.formula(paste("~", opt$condition_col))
)

# Pre-filter: remove genes with < 10 reads in total across all samples
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
message(sprintf("    Genes after low-count filter: %d", nrow(dds)))

# ── Run DESeq2 ────────────────────────────────────────────────────────────────
message("\n[4] Running DESeq2 …")
dds <- DESeq(dds)

# ── Extract results ───────────────────────────────────────────────────────────
message("\n[5] Extracting results …")
res <- results(
  dds,
  alpha       = opt$alpha,
  lfcThreshold = 0,              # initial extraction without threshold
  independentFiltering = TRUE
)

# Apply lfcShrink for more accurate fold-change estimates (apeglm method)
res_shrunk <- lfcShrink(
  dds,
  coef   = resultsNames(dds)[2],  # second coefficient = treatment vs reference
  type   = "apeglm",
  res    = res
)

# Convert to data.frame and add gene column
res_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj)

# Write full results
results_path <- file.path(opt$outdir, "DE_results.csv")
write.csv(res_df, results_path, row.names = FALSE)
message(sprintf("    DE results saved: %s", results_path))

sig_up   <- sum(res_df$padj < opt$alpha & res_df$log2FoldChange >  opt$lfc_threshold, na.rm = TRUE)
sig_down <- sum(res_df$padj < opt$alpha & res_df$log2FoldChange < -opt$lfc_threshold, na.rm = TRUE)
message(sprintf("    Significant (padj < %.2f, |LFC| > %.1f): %d up, %d down",
                opt$alpha, opt$lfc_threshold, sig_up, sig_down))

# ── Variance-stabilising transformation (VST) for visualisation ───────────────
message("\n[6] Computing VST normalised counts …")
vst_data <- vst(dds, blind = FALSE)
vst_mat  <- assay(vst_data)

norm_path <- file.path(opt$outdir, "normalized_counts.csv")
write.csv(as.data.frame(vst_mat) %>% tibble::rownames_to_column("gene"),
          norm_path, row.names = FALSE)
message(sprintf("    Normalised counts saved: %s", norm_path))

# ── Plot helpers ──────────────────────────────────────────────────────────────
THEME_BASE <- theme_bw(base_size = 12) +
  theme(
    plot.title   = element_text(size = 14, face = "bold"),
    legend.title = element_text(face = "bold")
  )

# ── MA Plot ───────────────────────────────────────────────────────────────────
message("\n[7] Generating MA plot …")
ma_df <- res_df %>%
  mutate(
    sig = case_when(
      padj < opt$alpha & log2FoldChange >  opt$lfc_threshold ~ "Up",
      padj < opt$alpha & log2FoldChange < -opt$lfc_threshold ~ "Down",
      TRUE ~ "NS"
    ),
    sig = factor(sig, levels = c("Up", "Down", "NS"))
  )

# Label the top 10 most significant genes
top_genes <- ma_df %>%
  filter(sig != "NS") %>%
  slice_min(padj, n = 10, with_ties = FALSE)

ma_plot <- ggplot(ma_df, aes(x = baseMean, y = log2FoldChange, colour = sig)) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_hline(yintercept = c(-opt$lfc_threshold, opt$lfc_threshold),
             linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, colour = "black") +
  scale_x_log10() +
  scale_colour_manual(values = c("Up" = "#d62728", "Down" = "#1f77b4", "NS" = "grey60")) +
  geom_label_repel(data = top_genes, aes(label = gene),
                   size = 2.5, max.overlaps = 20, show.legend = FALSE) +
  labs(
    title  = "MA Plot",
    x      = "Mean normalised counts (log10)",
    y      = "Log2 fold change",
    colour = "DE status"
  ) +
  THEME_BASE

ggsave(file.path(opt$outdir, "MA_plot.pdf"),
       ma_plot, width = 8, height = 6)
message("    MA plot saved.")

# ── Volcano Plot ──────────────────────────────────────────────────────────────
message("\n[8] Generating volcano plot …")
volcano_df <- ma_df %>%
  mutate(neg_log10_padj = -log10(pmax(padj, 1e-300)))

# Re-derive top_genes from volcano_df so neg_log10_padj column is available
top_genes_volcano <- volcano_df %>%
  filter(sig != "NS") %>%
  slice_min(padj, n = 10, with_ties = FALSE)

volcano_plot <- ggplot(volcano_df,
                       aes(x = log2FoldChange, y = neg_log10_padj, colour = sig)) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_vline(xintercept = c(-opt$lfc_threshold, opt$lfc_threshold),
             linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(opt$alpha),
             linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("Up" = "#d62728", "Down" = "#1f77b4", "NS" = "grey60")) +
  geom_label_repel(data = top_genes_volcano, aes(label = gene),
                   size = 2.5, max.overlaps = 20, show.legend = FALSE) +
  labs(
    title  = "Volcano Plot",
    x      = "Log2 fold change",
    y      = expression(-log[10](adjusted~p-value)),
    colour = "DE status"
  ) +
  THEME_BASE

ggsave(file.path(opt$outdir, "volcano_plot.pdf"),
       volcano_plot, width = 8, height = 6)
message("    Volcano plot saved.")

# ── Heatmap ───────────────────────────────────────────────────────────────────
message("\n[9] Generating heatmap (top 50 DE genes) …")

# Select top 50 significant genes by adjusted p-value
top50_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  slice_min(padj, n = 50, with_ties = FALSE) %>%
  pull(gene)

if (length(top50_genes) < 2) {
  message("    WARNING: Fewer than 2 significant genes; skipping heatmap.")
  file.create(file.path(opt$outdir, "heatmap.pdf"))
} else {
  # Subset and scale VST matrix (z-score per gene)
  heat_mat <- vst_mat[top50_genes, ]
  heat_mat_scaled <- t(scale(t(heat_mat)))  # scale rows

  # Annotation sidebar
  anno_col <- data.frame(
    Condition = col_data[[opt$condition_col]],
    row.names = rownames(col_data)
  )
  anno_colours <- list(
    Condition = setNames(
      colorRampPalette(brewer.pal(8, "Set2"))(nlevels(col_data[[opt$condition_col]])),
      levels(col_data[[opt$condition_col]])
    )
  )

  pdf(file.path(opt$outdir, "heatmap.pdf"), width = 10, height = 12)
  pheatmap(
    heat_mat_scaled,
    annotation_col  = anno_col,
    annotation_colors = anno_colours,
    color           = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    cluster_rows    = TRUE,
    cluster_cols    = TRUE,
    show_rownames   = TRUE,
    show_colnames   = TRUE,
    fontsize_row    = 7,
    fontsize_col    = 9,
    main            = "Top 50 Differentially Expressed Genes (VST, z-score)"
  )
  dev.off()
  message("    Heatmap saved.")
}

message("\n═══ DESeq2 analysis complete ✓ ═══")
