# rules/trim.smk  –  Read trimming with fastp

rule fastp:
    """
    Trim adapters and low-quality bases with fastp.
    Produces trimmed FASTQ files and a JSON/HTML QC report.
    """
    input:
        r1 = lambda wc: get_raw_read(wc, "read1"),
        r2 = lambda wc: get_raw_read(wc, "read2"),
    output:
        r1   = "{results}/trimmed/{sample}_R1.fastq.gz",
        r2   = "{results}/trimmed/{sample}_R2.fastq.gz",
        json = "{results}/trimmed/{sample}_fastp.json",
        html = "{results}/trimmed/{sample}_fastp.html",
    params:
        min_len = config["trimming"]["min_length"],
        quality = config["trimming"]["quality"],
    threads: config["trimming"]["threads"]
    conda: "../envs/trim.yaml"
    log: "{results}/logs/fastp/{sample}.log"
    shell:
        """
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --json {output.json} \
            --html {output.html} \
            --qualified_quality_phred {params.quality} \
            --length_required {params.min_len} \
            --detect_adapter_for_pe \
            --correction \
            --thread {threads} \
            2> {log}
        """
