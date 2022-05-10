# -*- coding: utf-8 -*-

import pandas as pd

def read_mutect_results(tsv):
    
    # Define column names.
    cols = [
        "Chr",
        "Pos",
        "Ref",
        "Alt",
        "Genotype",
        "AlleleFrequency",
        "RefReads",
        "AltReads"        
    ]    

    # Convert .tsv file to pd.DataFrame.
    df = pd.read_table(tsv, header=None, names=cols)
    df = df.fillna(".")
    
    # Set index.
    df = df.set_index(["Chr", "Pos", "Ref", "Alt"])
    
    # Return pd.DataFrame.
    return df

def read_varscan_results(tsv):

    # Define column names.
    cols = [
        "Chr",
        "Pos",
        "Ref",
        "Alt",
        "Genotype",
        "AlleleFrequency",
        "RefReads",
        "AltReads"        
    ]    

    # Convert .tsv file to pd.DataFrame.
    df = pd.read_table(tsv, header=None, names=cols)
    df = df.fillna(".")
    
    # VarScan2 reports allele frequencies as percentages.
    # Convert them to float numbers.
    df["AlleleFrequency"] = df["AlleleFrequency"].str.strip("%").astype("float") / 100
    
    # Set index.
    df = df.set_index(["Chr", "Pos", "Ref", "Alt"])
    
    # Return pd.DataFrame.
    return df

def read_strelka_results(tsv):

    # Define column names.
    cols = [
        "Chr",
        "Pos",
        "Ref",
        "Alt",
        "Genotype",
        "RefReads",
        "AltReads"        
    ]    

    # Convert .tsv file to pd.DataFrame.
    df = pd.read_table(tsv, header=None, names=cols)
    df = df.fillna(".")
    
    # Calculate allele frequency.
    df["AlleleFrequency"] = df["AltReads"] / (df["AltReads"] + df["RefReads"])
    
    # Set index.
    df = df.set_index(["Chr", "Pos", "Ref", "Alt"])
    
    # Return pd.DataFrame.
    return df

def merge_variant_calling_results(mutect, varscan, strelka, output):
    
    # Start with Mutect2 results.
    df = read_mutect_results(mutect)
    df["Tools"] = "mutect"
    
    # Continue with VarScan2 results.
    df_ = read_varscan_results(varscan)
    df_["Tools"] = "varscan"
    
    for index, row in df_.iterrows():
        if index in df.index.tolist():
            df.at[index, "Tools"] = df.at[index, "Tools"] + ",varscan"
        else:
            df = df.append(row)
    
    # Continue with Strelka2 results.
    df_ = read_strelka_results(strelka)
    df_["Tools"] = "strelka"
    
    for index, row in df_.iterrows():
        if index in df.index.tolist():
            df.at[index, "Tools"] = df.at[index, "Tools"] + ",strelka"
        else:
            df = df.append(row)

    # Round variant allele frequencies to same decimal accuracy.
    df["AlleleFrequency"] = round(df["AlleleFrequency"], 3)
    
    # Reset indexing.
    df = df.reset_index(drop=False)

    # Format chromosome to Ensembl compatible format.
    # Sort by position and chromosome.
    df["Chr"] = df["Chr"].str.lstrip("chr")
    df = df.sort_values(by=["Chr", "Pos"])

    # Save pd.DataFrame.
    df.to_csv(output, sep="\t", index=False)
    
# Run script.
merge_variant_calling_results(snakemake.input[0], snakemake.input[1], snakemake.input[2], snakemake.output[0])