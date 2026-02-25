# rules/deseq2.smk  –  DESeq2 differential expression

rule deseq2:
    """
    Run the DESeq2 R script.
    Inputs: clean counts matrix + samples CSV.
    Outputs: DE results table, normalized counts, and three plots.
    """
    input:
        counts  = "{results}/counts/counts_clean.txt",
        samples = config["samples"],
    output:
        results    = "{results}/deseq2/DE_results.csv",
        norm_counts = "{results}/deseq2/normalized_counts.csv",
        ma_plot    = "{results}/deseq2/MA_plot.pdf",
        volcano    = "{results}/deseq2/volcano_plot.pdf",
        heatmap    = "{results}/deseq2/heatmap.pdf",
    params:
        outdir          = "{results}/deseq2",
        condition_col   = config["deseq2"]["condition_column"],
        reference_level = config["deseq2"]["reference_level"],
        alpha           = config["deseq2"]["alpha"],
        lfc_threshold   = config["deseq2"]["lfc_threshold"],
    conda: "../envs/r.yaml"
    log: "{results}/logs/deseq2/deseq2.log"
    shell:
        """
        Rscript scripts/deseq2_analysis.R \
            --counts         {input.counts} \
            --samples        {input.samples} \
            --outdir         {params.outdir} \
            --condition_col  {params.condition_col} \
            --reference      {params.reference_level} \
            --alpha          {params.alpha} \
            --lfc_threshold  {params.lfc_threshold} \
            > {log} 2>&1
        """
