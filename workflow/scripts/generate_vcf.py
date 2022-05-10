# -*- coding: utf-8 -*-

from datetime import datetime
import pandas as pd

def generate_vcf(tsv, vcf):
    
    # Define header of VCF file.
    header = [
        "##fileformat=VCFv4.0",
        "##fileDate={}".format(datetime.today().strftime('%d-%m-%Y')),
        "##source=Tumor-only Variant Calling Workflow",
        "##reference=hg38",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>",
        "##FORMAT=<ID=AF,Number=1,Type=Float,Description='Variant allele frequency'>",
        "##FORMAT=<ID=AD,Number=1,Type=Integer,Description='Read Depths of reference and variant allele'>",
        "##FORMAT=<ID=TOOLS,Number=1,Type=String,Description='Supporting algorithms'>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    ]

    # Initialize VCF file and write VCF header.
    vcf = open(vcf, "x")
    for line in header:
        vcf.write(line + "\n")
    
    # Import variant information and write rows to VCF.
    variants = pd.read_table(tsv, sep="\t")
    
    # Loop through each variant.
    for index, row in variants.iterrows():
        
        # Extract information.
        chrom = str(row["Chr"])
        pos = str(row["Pos"])
        ref = str(row["Ref"])
        alt = str(row["Alt"])
        gt = str(row["Genotype"])
        af = str(row["AlleleFrequency"])
        ad = str(row["RefReads"]) + "," + str(row["AltReads"])
        tools = str(row["Tools"])
        
        # Generate variant line and write it to VCF.
        variant = "{}\t{}\t.\t{}\t{}\t.\tPASS\tGT:AF:AD:TOOLS\t{}:{}:{}:{}".format(chrom, pos, ref, alt, gt, af, ad, tools)
        vcf.write(variant + "\n")
    
    # Close VCF file.
    vcf.close()

# Run script.
generate_vcf(snakemake.input[0], snakemake.output[0])