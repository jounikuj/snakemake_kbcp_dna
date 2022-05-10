rule generate_report:
    input: 
        "output/{sample}/snpeff/{sample}_filtered_snpeff_clinvar.vcf"
    output: 
        "output/{sample}/{sample}_results.tsv"
    conda:
        "../envs/results.yaml"
    script:
        "../scripts/generate_report.py"