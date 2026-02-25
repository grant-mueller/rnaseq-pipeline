# rules/quantify.smk  –  featureCounts gene-level quantification

rule featurecounts:
    """
    Count reads per gene across all samples in a single featureCounts call.
    Output is a tab-separated matrix with gene IDs as rows and samples as columns.
    """
    input:
        bams = all_bams(),
        gtf  = config["genome"]["gtf"],
    output:
        counts  = "{results}/counts/counts_matrix.txt",
        summary = "{results}/counts/counts_matrix.txt.summary",
    params:
        strandedness = config["featurecounts"]["strandedness"],
        extra        = config["featurecounts"]["extra_flags"],
        threads      = config["featurecounts"]["threads"],
    conda: "../envs/quantify.yaml"
    log: "{results}/logs/featurecounts/featurecounts.log"
    shell:
        """
        featureCounts \
            -T {params.threads} \
            -a {input.gtf} \
            -o {output.counts} \
            -s {params.strandedness} \
            {params.extra} \
            {input.bams} \
            > {log} 2>&1
        """


rule simplify_counts:
    """
    Strip featureCounts metadata columns (Chr, Start, End, Strand, Length)
    leaving only GeneID + one count column per sample.
    The Python script also renames columns to clean sample names.
    """
    input:
        "{results}/counts/counts_matrix.txt",
    output:
        "{results}/counts/counts_clean.txt",
    conda: "../envs/python.yaml"
    log: "{results}/logs/featurecounts/simplify_counts.log"
    shell:
        """
        python scripts/simplify_counts.py \
            --input  {input} \
            --output {output} \
            > {log} 2>&1
        """
