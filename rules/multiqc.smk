# rules/multiqc.smk  –  aggregate QC with MultiQC

rule multiqc:
    """
    Collect all QC outputs from FastQC, fastp, STAR and featureCounts
    into a single interactive HTML report.
    """
    input:
        # Raw FastQC zips
        expand(
            "{results}/fastqc/raw/{sample}_{read}_fastqc.zip",
            results=RESULTS, sample=SAMPLE_NAMES, read=["R1", "R2"],
        ),
        # Trimmed FastQC zips
        expand(
            "{results}/fastqc/trimmed/{sample}_{read}_fastqc.zip",
            results=RESULTS, sample=SAMPLE_NAMES, read=["R1", "R2"],
        ),
        # fastp JSON reports
        expand(
            "{results}/trimmed/{sample}_fastp.json",
            results=RESULTS, sample=SAMPLE_NAMES,
        ),
        # STAR final logs
        expand(
            "{results}/star/{sample}/{sample}.Log.final.out",
            results=RESULTS, sample=SAMPLE_NAMES,
        ),
        # featureCounts summary
        f"{RESULTS}/counts/counts_matrix.txt.summary",
    output:
        html = "{results}/multiqc/multiqc_report.html",
    params:
        outdir     = "{results}/multiqc",
        search_dir = "{results}",
        cfg        = config["multiqc"]["config"],
    conda: "../envs/qc.yaml"
    log: "{results}/logs/multiqc/multiqc.log"
    shell:
        """
        multiqc \
            --config {params.cfg} \
            --outdir {params.outdir} \
            --filename multiqc_report \
            --force \
            {params.search_dir} \
            > {log} 2>&1
        """
