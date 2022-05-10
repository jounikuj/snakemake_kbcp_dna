rule run_picard_markduplicates:
    input:
        bam="output/{sample}/alignment/{sample}.bam",
        bai="output/{sample}/alignment/{sample}.bam.bai"
    output:
        bam="output/{sample}/alignment/{sample}_dedup.bam",
        bai="output/{sample}/alignment/{sample}_dedup.bam.bai",
    log:
        "output/{sample}/alignment/{sample}_markduplicates_metrics.txt"
    conda:
        "../envs/deduplication.yaml"
    shell:
        """
        picard MarkDuplicates -I {input.bam} -O {output.bam} -M {log}
        samtools index {output.bam}
        """