import pandas as pd
import vcf as vcf

def generate_oncokb_input(input_file, wildcards, output_file):

    # Import and open VCF file.
    vcf_reader = vcf.Reader(open(input_file, "r"))
    
    # Initialize list to store parsed information.
    variants = []
    
    # Loop through variants and extract required data.
    for variant in vcf_reader:
        
        # Basic information.
        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        alt = variant.ALT[0]
        
        # Gene name.
        gene = variant.INFO["ANN"][0].split("|")[3]
        
        # Aminoa acid change, will be converted to OncoKB compatible
        # format.
        aachange = variant.INFO["ANN"][0].split("|")[10]
        
        # Generate HGVSp_short column.
        aachanges = {
            "Gly":"G",
            "Ala":"A",
            "Val":"V",
            "Leu":"L",
            "Ile":"I",
            "Pro":"P",
            "Ser":"S",
            "Thr":"T",
            "Asn":"N",
            "Gln":"Q",
            "Cys":"C",
            "Met":"M",
            "Phe":"F",
            "Tyr":"Y",
            "Trp":"W",
            "Asp":"D",
            "Glu":"E",
            "His":"H",
            "Lys":"K",
            "Arg":"R"
        }
        
        aachange_short = aachange
        for long, short in aachanges.items():
            aachange_short = aachange_short.replace(long, short)

        # Save to list.
        variant = [
            gene,
            wildcards.sample,
            aachange,
            aachange_short,
            chrom,
            pos,
            ref,
            alt
        ]
        
        variants.append(variant)
        
    # Convert to pd.DataFrame.
    cols = [
        "Hugo_Symbol",
        "Tumor_Sample_Barcode",
        "HGVSp",
        "HGVSp_Short",
        "Chromosome",
        "Start_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele1"
    ]
    
    df = pd.DataFrame(variants, columns=cols)
    
    # Write output file.
    df.to_csv(output_file, sep="\t", index=False)

# Run script.
generate_oncokb_input(snakemake.input[0], snakemake.wildcards, snakemake.output[0])