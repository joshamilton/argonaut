# emilytrybulec/genomeassembly: Output

## Introduction

This document describes the output produced by the pipeline. All of the programs being run will have their own folder in the output directory specified by your params.yaml file.

The directories created after the pipeline has finished will depend on which options are selected in the configuration and which programs are run. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Busco](#busco) - Assembly quality checking for completeness
- [Centrifuge](#centrifuge) - Contaminant detection for long reads
- [Extract](#extract) - Extracting genome size estimation from kmerfreq and gce
- [Fastp](#fastp) - Adapter trimming for short reads
- [Fastqc](#fastqc) - Short read quality checking
- [Flye](#flye) - Flye assembly
- [GCE](#gce) - Genome size estimation using kmerfreq
- [Genomescope2](#genomescope2) - Short read [what does it do]
- [Gunzip](#gunzip) - Converting files from .gz to unzipped (no .gz)
- [Gzip](#gzip) - Converting files from unzipped to .gz
- [Jellyfish](#jellyfish) - Short read [what does it do]
- [Kmerfreq](#kmerfreq) - Genome size estimation using an inputted kmer size 
- [Kraken2](#kraken2) - Contaminant detection for short reads
- [MaSuRCA](#masurca) - Hybrid assembly
- [Medaka](#medaka) - Long read polishing of flye assembly
- [Merqury](#merqury) - Assembly quality checking for accuracy
- [Meryl](#meryl) - Building a database for merqury quality checking
- [Minimap2](#minimap2) - Assemblies with aligned reads
- [Nanoplot](#nanoplot) - Long read quality checking
- [Output](#output) - Final quality stats (quast, busco, merqury) of all assemblies produced
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
- [POLCA](#polca) - Short read polishing of flye assembly
- [Purge](#purge) - Purge haplotigs to reduce duplication
- [PycoQC](#pycoqc) - Assembly quality checking with sequencing summary
- [Quast](#quast) - Assembly quality checking for contiguity
- [Recentrifuge](#recentrifuge) - Contaminant filtering visualization for short and long reads
- [Samtools](#nanoplot) - Indexed Minimap2 alignments
- [Seqkit](#seqkit) - Filtering contaminants detected by Centrifuge from long reads for downstream processing


### busco

<details markdown="1">
<summary>Output files</summary>

- `busco/`
  - `canu_*/`: Canu assembly BUSCO output directory
  - `flye_*/`: Flye assembly BUSCO output directory
  - `masurca_*/`: MaSuRCA assembly BUSCO output directory
  - `polca_*/`: MaSuRCA assembly BUSCO output directory

</details>

[BUSCO](https://busco.ezlab.org/) gives completeness quality metrics about your assemblies. It provides information about whether your assembly matches with known protein sequences for similar species. Please be sure to specify the correct BUSCO lineage for your species. For further reading and documentation see the [BUSCO user guide](https://busco.ezlab.org/busco_userguide.html).

### centrifuge

<details markdown="1">
<summary>Output files</summary>

- `centrifuge/`
  - `centrifuge/`: Directory containting contaminant detection data
    -`filtered_*.report.txt/`
    -`filtered_*.results.txt/`
    -`filtered_*.unmapped.fastq.gz/`: Centrifuge filtered reads (not used for downstream processing, as filtering with seqkit is more accurate)
  - `kreport/`: Directory containing a kraken style report about the contaminants detected
    -`*.txt`

</details>

[Centrifuge](https://ccb.jhu.edu/software/centrifuge/) is a contaminant detection tool that generates a report summarizing which of your reads contain known bacterial, fungal, viral, etc. sequences. These reads are flagged and filtered from the reads used downstream, as they are likely contaminants. Please provide a centrifuge database containing contaminant sequences in params.yaml for the pipeline to run smoothly. For more information, please see the [Centrifuge manual](https://github.com/DaehwanKimLab/centrifuge/blob/master/MANUAL).
Centrifuge results are summarized in the Recentrifuge out directory.

### extract

<details markdown="1">
<summary>Output files</summary>

- `extract/`
  - `finalSize.txt`: Genome size estimation in megabases or gigabases
  - `standardSize.txt`: Genome size estimation in standard numerical format

</details>

The extract module uses awk and numfmt to isolate the genome size estimate for downstream use in assemblers like flye and canu, as well as quality checking with merqury.

### fastp

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `*1.fail.fastq.gz`
  - `*1.fastp.fastq.gz`
  - `*2.fail.fastq.gz`
  - `*2.fastp.fastq.gz`
  - `*fastp.html`
  - `*fastp.json`
  - `*fastp.log`
  - `*merged.fastq.gz`

</details>

[Fastp](https://github.com/OpenGene/fastp) does short read adapter trimming and provides quality metrics about your short reads. For further reading and documentation see the [Fastp usage documentation](https://github.com/OpenGene/fastp#simple-usage).

### fastqc

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `2_adaper_trim/`
    - `*1.fastqc.html`
    - `*1.fastqc.zip`
    - `*2.fastqc.html`
    - `*2.fastqc.zip`
  - `3_decontam/`
    - `*1.fastqc.html`
    - `*1.fastqc.zip`
    - `*2.fastqc.html`
    - `*2.fastqc.zip`
  - `*1.fastqc.html`
  - `*1.fastqc.zip`
  - `*2.fastqc.html`
  - `*2.fastqc.zip`

</details>

[Fastqc](https://github.com/s-andrews/FastQC) provides quality metrics about your short reads. It is run on the raw short reads, after adapter trimming with fastp, and after contaminant filtering with kraken2. For further reading and documentation see the [Fastqc help page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### flye

<details markdown="1">
<summary>Output files</summary>

- `flye/`
  - `assembly_graph.gfa`
  - `assembly_graph.gv`
  - `flye_*.assembly_info.txt`
  - `flye_*.fasta`
  - `flye_*.flye.log`
  - `flye_*.params.json`

</details>

[Flye](https://github.com/fenderglass/Flye) is a long read assembler. The fasta file is used downstream. For further reading and documentation see the [Flye usage documentation](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md).

### gce

<details markdown="1">
<summary>Output files</summary>

- `gce/`
  - `gce2.log`

</details>

[GCE](https://github.com/fanagislab/GCE) is a genomic character estimator. It uses the kmerfreq output to estimate genome size. 

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

### seqkit

<details markdown="1">
<summary>Output files</summary>

- `seqkit/`
  - `*_filtered.fastq`: Filtered long reads after Centrifuge contaminant detection

</details>

[seqkit](https://github.com/shenwei356/seqkit) is a tool for file manipulation and is used in the pipeline to filter flagged contaminants from long reads. For further reading and documentation see the [seqkit user guide](https://bioinf.shenwei.me/seqkit/usage/).