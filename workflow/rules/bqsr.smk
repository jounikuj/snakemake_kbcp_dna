rule get_target_regions:
    input:
        ref="resources/reference.dict",
        targets=config["sequencing"]["gene_panel"]
    output:
        "resources/targets.interval_list"
    conda:
        "../envs/bqsr.yaml"
    shell:
        "picard BedToIntervalList I={input.targets} O={output} SD={input.ref}"

rule get_dbsnp_variants:
    output:
        vcf="resources/dbsnp_variants.vcf.gz",
        tbi="resources/dbsnp_variants.vcf.gz.tbi"
    params:
        url=config["variant_calling"]["known_sites"]["dbsnp"]
    conda:
        "../envs/bqsr.yaml"
    threads:
        4
    shell:
        """
        wget -O - {params.url} | bgzip -c -@ {threads} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule get_mills_variants:
    output:
        vcf="resources/mills_variants.vcf.gz",
        tbi="resources/mills_variants.vcf.gz.tbi"
    params:
        url=config["variant_calling"]["known_sites"]["mills"]
    conda:
        "../envs/bqsr.yaml"
    threads:
        4
    shell:
        """
        wget -O - {params.url} | zcat - | bgzip -c -@ {threads} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

if config["deduplication"]:
    rule run_gatk_baserecalibrator:
        input:
            bam="output/{sample}/alignment/{sample}_dedup.bam",
            bai="output/{sample}/alignment/{sample}_dedup.bam.bai",
            ref=[
                "resources/reference.fasta.gz",
                "resources/reference.fasta.gz.fai",
                "resources/reference.dict"
            ],
            dbsnp=[
                "resources/dbsnp_variants.vcf.gz",
                "resources/dbsnp_variants.vcf.gz.tbi"
            ],
            mills=[
                "resources/mills_variants.vcf.gz",
                "resources/mills_variants.vcf.gz.tbi"
            ],
            targets="resources/targets.interval_list"
        output:
            recal="output/{sample}/alignment/{sample}_recal.table"
        conda:
            "../envs/bqsr.yaml"
        shell:
            "gatk BaseRecalibrator -I {input.bam} -R {input.ref[0]} --known-sites {input.dbsnp[0]} --known-sites {input.mills[0]} -L {input.targets} -O {output.recal}"
else:
    rule run_gatk_baserecalibrator:
        input:
            bam="output/{sample}/alignment/{sample}.bam",
            bai="output/{sample}/alignment/{sample}.bam.bai",
            ref=[
                "resources/reference.fasta.gz",
                "resources/reference.fasta.gz.fai",
                "resources/reference.dict"
            ],
            dbsnp=[
                "resources/dbsnp_variants.vcf.gz",
                "resources/dbsnp_variants.vcf.gz.tbi"
            ],
            mills=[
                "resources/mills_variants.vcf.gz",
                "resources/mills_variants.vcf.gz.tbi"
            ],
            targets="resources/targets.interval_list"
        output:
            recal="output/{sample}/alignment/{sample}_recal.table"
        conda:
            "../envs/bqsr.yaml"
        shell:
            "gatk BaseRecalibrator -I {input.bam} -R {input.ref[0]} --known-sites {input.dbsnp[0]} --known-sites {input.mills[0]} -L {input.targets} -O {output.recal}"

if config["deduplication"]:
    rule run_gatk_applybqsr:
        input:
            bam="output/{sample}/alignment/{sample}_dedup.bam",
            bai="output/{sample}/alignment/{sample}_dedup.bam.bai",
            bqsr="output/{sample}/alignment/{sample}_recal.table",
            ref=[
                "resources/reference.fasta.gz",
                "resources/reference.fasta.gz.fai",
                "resources/reference.fasta.gz.gzi"
            ],
            targets="resources/targets.interval_list"
        output:
            bam="output/{sample}/alignment/{sample}_bqsr.bam",
            bai="output/{sample}/alignment/{sample}_bqsr.bam.bai"
        conda:
            "../envs/bqsr.yaml"
        shell:
            """
            gatk ApplyBQSR -R {input.ref[0]} -I {input.bam} --bqsr-recal-file {input.bqsr} -L {input.targets} --create-output-bam-index false -O {output.bam}
            samtools index {output.bam} 
            """
else:
    rule run_gatk_applybqsr:
        input:
            bam="output/{sample}/alignment/{sample}.bam",
            bai="output/{sample}/alignment/{sample}.bam.bai",
            bqsr="output/{sample}/alignment/{sample}_recal.table",
            ref=[
                "resources/reference.fasta.gz",
                "resources/reference.fasta.gz.fai",
                "resources/reference.fasta.gz.gzi"
            ],
            targets="resources/targets.interval_list"
        output:
            bam="output/{sample}/alignment/{sample}_bqsr.bam",
            bai="output/{sample}/alignment/{sample}_bqsr.bam.bai"
        conda:
            "../envs/bqsr.yaml"
        shell:
            """
            gatk ApplyBQSR -R {input.ref[0]} -I {input.bam} --bqsr-recal-file {input.bqsr} -L {input.targets} --create-output-bam-index false -O {output.bam}
            samtools index {output.bam} 
            """