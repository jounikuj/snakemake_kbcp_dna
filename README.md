
# Tumor-only variant calling workflow

 
This is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for the somatic variant valling from targeted tumor-only NGS data. Workflow covers the preprocessing of NGS reads (trimming, alignment, optional deduplication and base quality score recalibration) and calls somatic and germline variants with [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) algorithm. Called variants are annotated with the [SnpEff](http://pcingola.github.io/SnpEff/).

### Authors
MSc Jouni Kujala, University of Eastern Finland, Kuopio, Finland

  

### Installation and dependencies

Following dependencies must be installed before running this workflow:

* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [Mamba](https://github.com/mamba-org/mamba)
* [Conda](https://docs.conda.io/en/latest/)

We recommend you to first install Conda and then install Mamba and Snakemake with following commands. These steps install Mamba to your local system and generate a virtual environment for Snakemake workflow. Actual workflow runs are performed within this environment. Rest of the dependencies are handled automatically by Snakemake.

```
# Install Mamba package manager.
conda install -c conda-forge mamba

# Use Mamba to create virtual environment for Snakemake workflow.
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Next step is to obtain your working copy of this workflow. You can either download this repository as .zip file or clone repository with `git`. Following command can be used to clone this repository:
```
# Go to directory to which the repository should be cloned.
cd path/to/target/directory

# Clone repository.
git clone https://github.com/jounikuj/tumor_only_variant_calling_workflow.git
```
  

### Prepare your run

Preparing your workflow run is straightforwards as it only contains three steps.

1. Import your FASTQ files to `input/` directory. Naming convention `{sample_id}_S1_R1_001.fastq.gz` should be used used for FASTQ files where the `{sample_id}` is the sample identifier used by workflow.
2. Import your gene panel design as BED file to `config/` directory.
3. Check that the configuration file `config/config.yaml` is correctly set up. This file defines the used adapter sequences and reference files (reference genome etc.) and it is essential to configurate workflow before your analysis.


### Start workflow

Go to workflow directory and activate Snakemake's virtual environment. It is recommended to test Snakemake in its built-in dryrun mode as it will usually identify possible errors in the workflow configuration.

```
# Activate virtual environment.
conda activate snakemake

# Run Snakemake in dryrun mode. 
snakemake --use-conda --cores 1 --dryrun
```

If everything seems to be fine, proceed to actual workflow by dropping the `--dryrun` mode from the command. Parameter `--cores` specifies the number of available cores for the workflows. Provided value must be at least one, but you can significantly speed up the workflow by providing more cores for the analysis.

```
# Run workflow.
snakemake --use-conda --cores [int]
```

Workflow should now start. It will take some time to generate all required virtual environments especially if this is your first run. 

### Workflow description

1. Adapter and primer sequences are first trimmed from raw FASTQ reads with [cutadapt](https://cutadapt.readthedocs.io/en/stable/). By default, reads are quality-trimmed so that the bases with base quality scores less than Q10 are trimmed from 3' end of the read and reads shorter than 10 bp are discarded from further analysis.
2. Trimmed reads are aligned against human reference genome hg38 with [BWA-MEM](http://bio-bwa.sourceforge.net/bwa.shtml) algorithm.
3. Aligned reads are deduplicated with [Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-) tool to identify possible PCR artifacts and optical sequencing errors. Although deduplication is performed by default, it can be disabled from `config/config.yaml` if required by specific sequencing approaches (such as amplicon based target enrichment).
4. Base quality score recalibration (BQSR) is applied to detect and compensate systematic errors in the sequencing data that originate from sequencing machine.
5. Somatic and germline variants are called with [GATK4 Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2). Low-quality somatic variants are subsequently discarded by [GATK4 FilterMutectCalls](https://gatk.broadinstitute.org/hc/en-us/articles/360036856831-FilterMutectCalls) to obtain high-quality variant calls.
6. Variant calls are annotated with gene name, amino acid changes, functional consequence etc. beneficial information by using [SnpEff](http://pcingola.github.io/SnpEff/). Variant calls are further annotated by comparing variants to existing [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) and [OncoKB](https://www.oncokb.org/) database records.
9. Key results of variant calling workflow are parsed to Excel file.

### Bugs and feature requests
Found a bug from source code or have idea for improvements and new features? Open a new Issue or make a new Pull request through GitHub to help improve the development of this workflow!

### How to cite?
This workflow has not been published yet. If you decide to use this workflow in your research, simply refer to this GitHub repository.