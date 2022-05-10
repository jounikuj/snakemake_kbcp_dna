# Tumor-only variant calling workflow
# (C) MSc Jouni Kujala, University of Eastern Finland

import os
import pandas as pd

# Define and validate configuration file.
configfile: "config/config.yaml"

# Identify input files from input directory.
SAMPLES, INDEX, = glob_wildcards("input/{sample}_S{index}_R1_001.fastq.gz")

# Target rule.
rule all:
    input:
        expand("output/{sample}/{sample}_results.tsv", sample=SAMPLES),
        expand("output/{sample}/fastqc/{sample}_bqsr_fastqc.html", sample=SAMPLES)

# Trim reads and align them to reference genome.
include: "workflow/rules/trimming.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/deduplication.smk"
include: "workflow/rules/bqsr.smk"

# Perform quality control.
include: "workflow/rules/qc.smk"

# Run variant calling algorithms.
include: "workflow/rules/mutect.smk"

# Annotate called variants with SnpEff and SnpSift.
include: "workflow/rules/annotation.smk"

# Generate report.
include: "workflow/rules/report.smk"