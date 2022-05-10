import pandas as pd

def parse_oncokb_output(input_file, output_file):

    # Import OncoKB output file as pd.DataFrame.
    df = pd.read_table(input_file, sep="\t")

    # Convert chromosome to categorical variable.
    chrom = [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
        "chrM"
    ]
    
    df["Chromosome"] = pd.Categorical(df["Chromosome"], ordered=True, categories=chrom)

    # Generate temporary column for end position in case that there are multiple mutations
    # per one start position.
    df["End_Position"] = df["Start_Position"] + (df["Tumor_Seq_Allele1"].str.len() - 1)

    # Sort by genomic coordinates.
    df = df.sort_values(by=["Chromosome", "Start_Position"])

    # Select columns.
    cols = [
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele1",
        "VARIANT_IN_ONCOKB",
        "MUTATION_EFFECT",
        "ONCOGENIC"
    ]

    df = df[cols]

    # Convert start end end positions to 0-based format.
    for col in ["Start_Position", "End_Position"]:
        df[col] = df[col] - 1

    # Fill missing values.
    for col in ["VARIANT_IN_ONCOKB", "MUTATION_EFFECT", "ONCOGENIC"]:
        df[col] = df[col].fillna(".")
    
    # Save as tab-delimited file.
    df.to_csv(output_file, sep="\t", index=False, header=False)

# Run script.
parse_oncokb_output(snakemake.input[0], snakemake.output[0])