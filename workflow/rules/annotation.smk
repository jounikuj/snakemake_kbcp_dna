rule get_snpeff_database:
    """
    Download up-to-date database to used virtual environment.
    """
    output:
        temp("resources/snpeff/lock.txt")
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        snpEff download -v hg38
        touch {output}
        """

rule run_snpeff:
    """
    Run SnpEff and annotate detected variants with basic information, such
    as affected gene, amino acid change etc.
    """
    input:
        vcf="output/{sample}/mutect/{sample}_filtered.vcf",
        lock="resources/snpeff/lock.txt"
    output:
        vcf="output/{sample}/snpeff/{sample}_filtered_snpeff.vcf"
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        snpEff ann -noStats -canon hg38 {input.vcf} > {output.vcf}
        """

rule download_clinvar_records:
    """
    Download ClinVar records from NCBI FTP server. Downloaded file is compressed
    and indexed after downloading.
    """
    output:
        vcf="resources/clinvar.vcf.gz"
    params:
        url=config["annotation"]["clinvar"]
    conda:
        "../envs/annotation.yaml"
    threads:
        4
    shell:
        """
        wget -O - {params.url} | zcat - | bgzip -c -@ {threads} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule add_clinvar_annotation:
    """
    Compare ClinVar records to detected mutations and add ClinVar annotation to
    VCF file.
    """
    input:
        vcf="output/{sample}/snpeff/{sample}_filtered_snpeff.vcf",
        clinvar="resources/clinvar.vcf.gz"
    output:
        vcf="output/{sample}/snpeff/{sample}_filtered_snpeff_clinvar.vcf"
    conda:
        "../envs/annotation.yaml"
    shell:
        "SnpSift annotate -info CLNSIG {input.clinvar} {input.vcf} > {output}"

rule generate_oncokb_input:
    """
    Generates MAF file that can be used to run OncoKB MafAnnotator.
    """
    input:
        "output/{sample}/snpeff/{sample}_filtered_snpeff_clinvar.vcf"
    output:
        "output/{sample}/oncokb/{sample}_oncokb.input"
    conda:
        "../envs/annotation.yaml"
    script:
        "../scripts/generate_oncokb_input.py"

rule run_oncokb:
    """
    Run OncoKB MafAnnotator to retrieve OncoKB database matches.
    """
    input:
        "output/{sample}/oncokb/{sample}_oncokb.input"
    output:
        "output/{sample}/oncokb/{sample}_oncokb.output"
    params:
        token=config["annotation"]["oncokb"]["token"],
        cancer=config["annotation"]["oncokb"]["cancer"]
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        MafAnnotator.py -i {input} -o {output} -b {params.token} -t {params.cancer}
        """

rule parse_oncokb_output:
    """
    Generates MAF file that can be used to run OncoKB MafAnnotator.
    """
    input:
        "output/{sample}/oncokb/{sample}_oncokb.output"
    output:
        "output/{sample}/oncokb/{sample}_oncokb.bed"
    conda:
        "../envs/annotation.yaml"
    script:
        "../scripts/parse_oncokb_output.py"

rule add_oncokb_annotation:
    """
    Add OncoKB annotation to VCF file.
    """
    input:
        vcf="output/{sample}/snpeff/{sample}_filtered_snpeff_clinvar.vcf",
        bed="output/{sample}/oncokb/{sample}_oncokb.bed",
        header="resources/oncokb/header.txt"
    output:
        bed=[
            temp("output/{sample}/oncokb/{sample}_oncokb.bed.gz"),
            temp("output/{sample}/oncokb/{sample}_oncokb.bed.gz.tbi")
        ],
        vcf=[
            temp("output/{sample}/oncokb/{sample}_filtered_snpeff_clinvar.vcf.gz"),
            temp("output/{sample}/oncokb/{sample}_filtered_snpeff_clinvar.vcf.gz.tbi")
        ],
        oncokb="output/{sample}/oncokb/{sample}_filtered_snpeff_clinvar_oncokb.vcf"
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        bgzip -c {input.bed} > {output.bed[0]}
        tabix -p bed {output.bed[0]}
        bgzip -c {input.vcf} > {output.vcf[0]}
        tabix -p vcf {output.vcf[0]}
        bcftools annotate -a {output.bed[0]} -c CHROM,POS,-,REF,ALT,OncoKB_variant,OncoKB_function,OncoKB_oncogenic -h {input.header} {output.vcf[0]} > {output.oncokb}
        """