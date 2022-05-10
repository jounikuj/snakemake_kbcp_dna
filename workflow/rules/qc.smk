rule run_fastqc:
    input:
        "output/{sample}/alignment/{sample}_bqsr.bam"
    output:
        "output/{sample}/fastqc/{sample}_bqsr_fastqc.html"
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc --outdir output/{wildcards.sample}/fastqc {input}"