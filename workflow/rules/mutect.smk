rule get_panel_of_normals:
    output:
        "resources/panel_of_normals.vcf.gz"
    params:
        url=config["variant_calling"]["panel_of_normals"]
    conda:
        "../envs/mutect.yaml"
    threads:
        4
    shell:
        """
        wget -O - {params.url} | zcat - | awk '{{gsub(/^chr/,""); print}}' - | bgzip -c -@ {threads} > {output}
        tabix -p vcf {output}
        """

rule get_gnomad_variants:
    output:
        vcf="resources/gnomad_variants.vcf.gz",
        tbi="resources/gnomad_variants.vcf.gz.tbi"
    params:
        url=config["variant_calling"]["known_sites"]["gnomad"]
    conda:
        "../envs/mutect.yaml"
    threads:
        4
    shell:
        """
        wget -O - {params.url} | zcat - | bgzip -c -@ {threads} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule run_gatk_mutect:
    input:
        bam="output/{sample}/alignment/{sample}_bqsr.bam",
        bai="output/{sample}/alignment/{sample}_bqsr.bam.bai",
        ref=[
            "resources/reference.fasta.gz",
            "resources/reference.fasta.gz.fai",
            "resources/reference.fasta.gz.gzi",
            "resources/reference.dict",
            "resources/reference.fasta.gz",
        ],      
        targets="resources/targets.interval_list",
        germline="resources/gnomad_variants.vcf.gz",
        pon="resources/panel_of_normals.vcf.gz"
    output:
        vcf="output/{sample}/mutect/{sample}_nonfiltered.vcf"
    conda:
        "../envs/mutect.yaml"
    shell:
        "gatk Mutect2 -R {input.ref[0]} -I {input.bam} -L {input.targets} --germline-resource {input.germline} --panel-of-normals {input.pon} -O {output.vcf}"

rule run_gatk_filtermutectcalls:
    input:
        vcf="output/{sample}/mutect/{sample}_nonfiltered.vcf",
        ref=[
            "resources/reference.fasta.gz",
            "resources/reference.fasta.gz.fai",
            "resources/reference.fasta.gz.gzi",
            "resources/reference.dict",
            "resources/reference.fasta.gz"
        ]
    output:
        vcf="output/{sample}/mutect/{sample}_filtered.vcf"
    conda:
        "../envs/mutect.yaml"
    shell:
        "gatk FilterMutectCalls -V {input.vcf} -R {input.ref[0]} -O {output.vcf}"