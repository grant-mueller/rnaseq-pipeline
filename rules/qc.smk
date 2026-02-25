rule fastqc_raw:
    input:
        r1 = lambda wc: get_raw_read(wc, "read1"),
        r2 = lambda wc: get_raw_read(wc, "read2"),
    output:
        html_r1 = "{results}/fastqc/raw/{sample}_R1_fastqc.html",
        zip_r1  = "{results}/fastqc/raw/{sample}_R1_fastqc.zip",
        html_r2 = "{results}/fastqc/raw/{sample}_R2_fastqc.html",
        zip_r2  = "{results}/fastqc/raw/{sample}_R2_fastqc.zip",
    params:
        outdir = "{results}/fastqc/raw",
    threads: 2
    conda: "../envs/qc.yaml"
    log: "{results}/logs/fastqc_raw/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc --threads {threads} --outdir {params.outdir} \
            {input.r1} {input.r2} > {log} 2>&1

        r1_base=$(basename {input.r1} .fastq.gz)
        r2_base=$(basename {input.r2} .fastq.gz)

        mv {params.outdir}/${{r1_base}}_fastqc.html {output.html_r1}
        mv {params.outdir}/${{r1_base}}_fastqc.zip  {output.zip_r1}
        mv {params.outdir}/${{r2_base}}_fastqc.html {output.html_r2}
        mv {params.outdir}/${{r2_base}}_fastqc.zip  {output.zip_r2}
        """


rule fastqc_trimmed:
    input:
        r1 = "{results}/trimmed/{sample}_R1.fastq.gz",
        r2 = "{results}/trimmed/{sample}_R2.fastq.gz",
    output:
        html_r1 = "{results}/fastqc/trimmed/{sample}_R1_fastqc.html",
        zip_r1  = "{results}/fastqc/trimmed/{sample}_R1_fastqc.zip",
        html_r2 = "{results}/fastqc/trimmed/{sample}_R2_fastqc.html",
        zip_r2  = "{results}/fastqc/trimmed/{sample}_R2_fastqc.zip",
    params:
        outdir = "{results}/fastqc/trimmed",
    threads: 2
    conda: "../envs/qc.yaml"
    log: "{results}/logs/fastqc_trimmed/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc --threads {threads} --outdir {params.outdir} \
            {input.r1} {input.r2} > {log} 2>&1
        """
