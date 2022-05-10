def get_fastq_files(wildcards):

    # Initialize dictionary for .fastq files.
    fastqs = {
        "r1": None,
        "r2": None
    }

    # Loop through all files in the input directory and identify files that match
    # 1) with the .fastq.gz extension and 2) given wildcard.
    for file in os.listdir("input/"):
        if (wildcards.sample in file) and (file.endswith(".fastq.gz")):

            # Check the direction of reads - is this R1 or R2 .fastq file?
            if "R1" in file:
                fastqs["r1"] = "input/" + file
            elif "R2" in file:
                fastqs["r2"] = "input/" + file

    # Confirm that both .fastq files were correctly identified.
    if (fastqs["r1"] is None) or (fastqs["r2"] is None):
        raise ValueError("Warning: workflow identified only one .fastq file for this sample!")

    # Return .fastq files.
    return fastqs
    
rule trim_reads:
    input:
        unpack(get_fastq_files)
    output:
        r1="output/{sample}/trimming/{sample}_R1_trimmed.fastq.gz",
        r2="output/{sample}/trimming/{sample}_R2_trimmed.fastq.gz"
    params:
        adapter_fwd=config["trimming"]["adapter_fwd"],
        adapter_rev=config["trimming"]["adapter_rev"],
        min_quality=config["trimming"]["min_quality"],
        min_length=config["trimming"]["min_length"]
    conda:
        "../envs/trimming.yaml"
    threads:
        4
    shell:
        """
        cutadapt \
        -j {threads} \
        -q {params.min_quality} \
        -m {params.min_length} \
        -a {params.adapter_fwd} \
        -A {params.adapter_rev} \
        -o {output.r1} \
        -p {output.r2} \
        {input.r1} {input.r2}
        """