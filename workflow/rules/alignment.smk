rule get_reference_genome:
    output:
        fa="resources/reference.fasta.gz"
    params:
        url=config["reference"]["fasta"]
    conda:
        "../envs/alignment.yaml"
    threads:
        4
    shell:
        "wget -O - {params.url} | bgzip -c -@ {threads} > {output.fa}"

rule get_reference_indeces:
    input:
        "resources/reference.fasta.gz"
    output:
        "resources/reference.fasta.gz.fai",
        "resources/reference.fasta.gz.gzi"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools faidx {output}"

rule get_reference_dict:
    input:
        "resources/reference.fasta.gz"
    output:
        "resources/reference.dict"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools dict {input} -o {output}"

rule get_bwa_indeces:
    input:
        fasta="resources/reference.fasta.gz",
        fai="resources/reference.fasta.gz.fai",
        gzi="resources/reference.fasta.gz.gzi",
        dict="resources/reference.dict"
    output:
        indeces=[
            "resources/reference.fasta.gz.amb",
            "resources/reference.fasta.gz.ann",
            "resources/reference.fasta.gz.bwt",
            "resources/reference.fasta.gz.pac",
            "resources/reference.fasta.gz.sa"
        ]
    conda:
        "../envs/alignment.yaml"
    shell:
        "bwa index resources/reference.fasta.gz {input.fasta}"

rule align_reads:
    input:
        r1="output/{sample}/trimming/{sample}_R1_trimmed.fastq.gz",
        r2="output/{sample}/trimming/{sample}_R2_trimmed.fastq.gz",
        fasta="resources/reference.fasta.gz",
        fai="resources/reference.fasta.gz.fai",
        gzi="resources/reference.fasta.gz.gzi",
        dict="resources/reference.dict",
        indeces=[
            "resources/reference.fasta.gz.amb",
            "resources/reference.fasta.gz.ann",
            "resources/reference.fasta.gz.bwt",
            "resources/reference.fasta.gz.pac",
            "resources/reference.fasta.gz.sa"
        ]
    output:
        bam="output/{sample}/alignment/{sample}.bam",
        bai="output/{sample}/alignment/{sample}.bam.bai"
    params:
        rg=r'@RG\tID:{sample}\tLB:{sample}\tPL:illumina\tPM:nextseq\tSM:{sample}' 
    conda:
        "../envs/alignment.yaml"
    threads:
        4
    shell:
        """
        bwa mem -R '{params.rg}' -t {threads} {input.fasta} {input.r1} {input.r2} | samtools view -hb -@ {threads} - | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """