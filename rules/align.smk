# rules/align.smk  –  STAR genome indexing + alignment + duplicate marking
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


rule markdup:
    """
    Mark PCR duplicates with samtools markdup.
    Requires name-sorted → fixmate → coord-sorted → markdup pipeline.
    Outputs a duplicate-marked BAM and a metrics file parsed by MultiQC.
    The original STAR BAM is preserved; featureCounts uses the markdup BAM.
    """
    input:
        bam = "{results}/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        bam     = "{results}/star/{sample}/{sample}.markdup.bam",
        metrics = "{results}/star/{sample}/{sample}.markdup.metrics.txt",
    threads: 4
    conda: "../envs/align.yaml"
    log: "{results}/logs/markdup/{sample}.log"
    shell:
        """
        # Step 1 – sort by name (required for fixmate)
        samtools sort -n \
            -@ {threads} \
            -o {wildcards.results}/star/{wildcards.sample}/{wildcards.sample}.namesorted.bam \
            {input.bam} \
            >> {log} 2>&1

        # Step 2 – fill in mate-score tags (required for markdup)
        samtools fixmate \
            -m \
            -@ {threads} \
            {wildcards.results}/star/{wildcards.sample}/{wildcards.sample}.namesorted.bam \
            {wildcards.results}/star/{wildcards.sample}/{wildcards.sample}.fixmate.bam \
            >> {log} 2>&1

        # Step 3 – re-sort by coordinate
        samtools sort \
            -@ {threads} \
            -o {wildcards.results}/star/{wildcards.sample}/{wildcards.sample}.coordsorted.bam \
            {wildcards.results}/star/{wildcards.sample}/{wildcards.sample}.fixmate.bam \
            >> {log} 2>&1

        # Step 4 – mark duplicates and write metrics
        samtools markdup \
            -@ {threads} \
            -s \
            -f {output.metrics} \
            {wildcards.results}/star/{wildcards.sample}/{wildcards.sample}.coordsorted.bam \
            {output.bam} \
            >> {log} 2>&1

        # Step 5 – index the final markdup BAM
        samtools index {output.bam} >> {log} 2>&1

        # Step 6 – clean up intermediates
        rm {wildcards.results}/star/{wildcards.sample}/{wildcards.sample}.namesorted.bam \
           {wildcards.results}/star/{wildcards.sample}/{wildcards.sample}.fixmate.bam \
           {wildcards.results}/star/{wildcards.sample}/{wildcards.sample}.coordsorted.bam
        """
