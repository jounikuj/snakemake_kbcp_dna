# -*- coding: utf-8 -*-

import pandas as pd
import vcf

def generate_report(input_file, output_file):
    
    # Import and open VCF file.
    vcf_reader = vcf.Reader(open(input_file, "r"))
        
    # Initialize list to store parsed information.
    variants = []
        
    # Loop through variants and extract required data.
    for variant in vcf_reader:
                
        # Basic information.
        chrom = variant.CHROM
        start = variant.POS
        end = variant.POS + len(variant.ALT[0]) - 1
        ref = variant.REF
        alt = variant.ALT[0]
        
        # Gene name, nucleotide change and AA change.
        gene = variant.INFO["ANN"][0].split("|")[3]
        hgvs_c = variant.INFO["ANN"][0].split("|")[9]
        hgvs_p = variant.INFO["ANN"][0].split("|")[10]
        
        # Mutation effect.
        effect = variant.INFO["ANN"][0].split("|")[1]
        impact = variant.INFO["ANN"][0].split("|")[2]
        
        # Mutation source.
        if len(variant.FILTER) == 0:
            source = "SOMATIC"
        else:  
            if variant.FILTER[0] == "germline":
                source = "GERMLINE"
            else:
                source = "AMBIGUOUS"
        
        # Allele frequency.
        reference_reads = variant.samples[0]["AD"][0]
        variant_reads = variant.samples[0]["AD"][1]
        depth = variant.samples[0]["DP"]
        af = round(float(variant_reads / depth), 3)
        
        # Extract ClinVar annotation.
        if "CLNSIG" in variant.INFO.keys():
            clinvar = variant.INFO["CLNSIG"][0]
        else:
            clinvar = "."
        
        # Check if variant has OncoKB record.
        if "OncoKB_variant" in variant.INFO.keys():
            oncokb_variant = variant.INFO["OncoKB_variant"][0]
        else:
            oncokb_variant = "."
            
        # Check functional OncoKB annotation.
        if "OncoKB_function" in variant.INFO.keys():
            oncokb_function = variant.INFO["OncoKB_function"][0]
        else:
            oncokb_function = "."
            
        # Check oncogenic OncoKB annotation.
        if "OncoKB_oncogenic" in variant.INFO.keys():
            oncokb_oncogenic = variant.INFO["OncoKB_oncogenic"][0]
        else:
            oncokb_oncogenic = "."
            
        # Add to variants.
        data = [
            chrom,
            start,
            end,
            ref,
            alt,
            gene,
            effect,
            impact,
            source,
            hgvs_c,
            hgvs_p,
            reference_reads,
            variant_reads,
            af,
            depth,
            clinvar,
            oncokb_variant,
            oncokb_function,
            oncokb_oncogenic
        ]
        
        variants.append(data)
        
    # Convert to pd.DataFrame.
    cols = [
        "Chromosome",
        "Start_position",
        "End_position",
        "Reference_allele",
        "Variant_allele",
        "Gene",
        "Effect",
        "Predicted impact",
        "Source",
        "HGVS_c",
        "HGVS_p",
        "Reference_reads",
        "Variant_reads",
        "Allele_frequency",
        "Read_depth",
        "ClinVar",
        "OncoKB_match",
        "OncoKB_functionality",
        "OncoKB_oncogenicity"
    ]
    
    df = pd.DataFrame(variants, columns=cols)
        
    # Write to tab-delimited TSV.
    df.to_csv(output_file, sep="\t", index=False)
        
# Run script.
generate_report(snakemake.input[0], snakemake.output[0])