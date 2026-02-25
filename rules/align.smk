# rules/align.smk  –  STAR genome indexing + alignment

rule star_index:
    """
    Build STAR genome index.
    Only runs when the index directory does not yet exist.
    """
    input:
        fasta = config["genome"]["fasta"],
        gtf   = config["genome"]["gtf"],
    output:
        # Touch a sentinel file so Snakemake tracks completion
        directory(config["genome"]["star_index"]),
    params:
        sjdb_overhang = config["genome"]["sjdb_overhang"],
    threads: config["star"]["threads"]
    conda: "../envs/align.yaml"
    log: f"{RESULTS}/logs/star_index/star_index.log"
    shell:
        """
        mkdir -p {output}
        STAR \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.sjdb_overhang} \
            --runThreadN {threads} \
            > {log} 2>&1
        """


rule star_align:
    """
    Align paired-end trimmed reads to the genome with STAR.
    Outputs a coordinate-sorted BAM and STAR log files (picked up by MultiQC).
    """
    input:
        r1    = "{results}/trimmed/{sample}_R1.fastq.gz",
        r2    = "{results}/trimmed/{sample}_R2.fastq.gz",
        index = config["genome"]["star_index"],
    output:
        bam   = "{results}/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        log   = "{results}/star/{sample}/{sample}.Log.final.out",
    params:
        prefix      = "{results}/star/{sample}/{sample}.",
        extra_flags = config["star"]["extra_flags"],
    threads: config["star"]["threads"]
    conda: "../envs/align.yaml"
    log: "{results}/logs/star/{sample}.log"
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix} \
            {params.extra_flags} \
            > {log} 2>&1

        # Index the BAM for downstream tools
        samtools index {output.bam}
        """
