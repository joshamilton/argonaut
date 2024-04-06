[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/genomeassembly)

## Introduction


**Argonaut** performs **a**utomated **r**eads to **g**enome **o**perations for de **n**ovo **a**ssemblies; it is a bioinformatics pipeline that performs genome assembly on long and short read data. A fastq file and input information is fed to the pipeline, resulting in final assemblies with completeness, contiguity, and correctnesss quality checking at each step. The pipeline accepts short reads, long reads, or both. 

<img align="right" height="300" src="https://github.com/emilytrybulec/genomeassembly/assets/114685119/9b900dab-44cb-479e-9362-0c0d9dc00ae0">

## Table of Contents
- [Pipeline Summary](#Pipeline-Summary)
- [Quick Start](#Quick-Start)
- [Output Overview](#Pipeline-Output)
- [Credits](#Credits)
- [Contributions & Support](#Contributions-and-Support)
- [Citations](#Citations)
   
## Pipeline Summary

Illumina Short Read 
1. Read QC, Adaptor Trimming, Contaminant Filtering([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [`FastP`](https://github.com/OpenGene/fastp), [`GenomeScope`](http://qb.cshl.edu/genomescope/),[`Jellyfish`](https://github.com/gmarcais/Jellyfish),[`Kraken2`](https://ccb.jhu.edu/software/kraken2/), [`Recentrifuge`](https://github.com/khyox/recentrifuge),)
    
ONT Long Read
1. Read QC and Contaminant Filtering([`Nanoplot`](https://github.com/wdecoster/NanoPlot),[`KmerFreq`](https://github.com/fanagislab/kmerfreq), [`GCE`](https://github.com/fanagislab/GCE), [`Centrifuge`](https://ccb.jhu.edu/software/centrifuge/), [`Recentrifuge`](https://github.com/khyox/recentrifuge))
2. Length Filtering (optional)([`Bioawk`](https://github.com/lh3/bioawk), [`Nanoplot`](https://github.com/wdecoster/NanoPlot))
   
PacBio Hifi Long Read (CCS format)
1. Read QC, Adaptor Trimming, Contaminant Filtering([`Nanoplot`](https://github.com/wdecoster/NanoPlot),[`CutAdapt`](https://cutadapt.readthedocs.io/en/stable/),[`GenomeScope`](http://qb.cshl.edu/genomescope/),[`Jellyfish`](https://github.com/gmarcais/Jellyfish),[`Kraken2`](https://ccb.jhu.edu/software/kraken2/), [`Recentrifuge`](https://github.com/khyox/recentrifuge))
2. Length Filtering (optional)([`Bioawk`](https://github.com/lh3/bioawk), [`Nanoplot`](https://github.com/wdecoster/NanoPlot))

All reads are used for the following steps:  
  
3. Assembly &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;             4. Assembly QC**
- [`Flye`](https://github.com/fenderglass/Flye)  &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;- [`BUSCO`](https://busco.ezlab.org/)
- [`Canu`](https://github.com/marbl/canu) &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;&nbsp; - [`Quast`](https://quast.sourceforge.net/)
- [`Hifiasm`](https://github.com/chhylp123/hifiasm) &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;  - [`Minimap2`](https://github.com/lh3/minimap2)
- [`MaSuRCA`](https://github.com/alekseyzimin/masurca) &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;  - [`PycoQC`](https://github.com/a-slide/pycoQC)
- [`Redundans`](https://github.com/Gabaldonlab/redundans) &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;  - [`Merqury`](https://github.com/marbl/merqury)  
  
5. Polish  &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;&nbsp;          6. Polish QC**  
- [`Medaka`](https://github.com/nanoporetech/medaka)
- [`Racon`](https://github.com/isovic/racon)
- [`POLCA`](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007981)
  
7. Purge  &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;&nbsp;   8. Purge QC**  
- [`PurgeHaplotigs`](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/)
- [`Redundans`](https://github.com/Gabaldonlab/redundans)
  
9. Scaffolding &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;     10. Scaffold QC**  
- ([`RagTag`](https://github.com/malonge/RagTag))  

11. Assembly Visualization
- ([`Blobtools`](https://blobtoolkit.genomehubs.org/blobtools2/))  

Below is a figure detailing the major workflow steps involved in hybrid assembly.

<img align="center" width="700" alt="Argonaut Hybrid Workflow" src="https://github.com/emilytrybulec/argonaut/assets/114685119/54fd9e52-d657-4a29-99a6-953f736e1630">

    
Purge Haplotigs is the first step of manual curation, as it produces a histogram that needs to be analyzed for -l, -m, -h flags. The pipeline will stop at the purge step if purge is activated in your configuration and wait for manual input of parameters according to the histogram of your assembly, which can be found in your out directory.

## Quick Start

First, prepare a samplesheet with your input data as follows:

`ont_samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,single_end
maca_jans_ont,SRR11191910.fastq.gz,,TRUE
```

If more than one read input type is available, prepare a second (and third) samplesheet with your other input data as follows:

`illumina_samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,single_end
maca_jans_ill,SRR11191912_1.fastq.gz,SRR11191912_2.fastq.gz,FALSE
```

`pb_hifi_samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,single_end
maca_jans_pb,SRR11191909.fastq.gz,,TRUE
```

!!! PLEASE ADD "ont", "pb", AND/OR "ill" TO YOUR SAMPLES NAMES !!! Failure to do so will result in assemblers not recognizing your read type.

Additionally, the sample name inputted in your samplesheet will serve as the prefix for your output files. Please indicate which kind of read is being inputted in the sample name. Failure to do so may result in outputs being overwritten. 



After you have your samplesheet(s), create a params.yaml file to specify the paths to your samplesheet(s), contaminant databases, etc. Most likely, a config file will also need to be made to modify the default settings of the pipeline. Please look through the [nextflow.config](nextflow.config) file to browse the defaults and specify which you would like to change in your my_config file. More information is located in the [usage](docs/usage.md) section.

Now, you can run the pipeline using:


```bash
nextflow run emilytrybulec/argonaut \
  -r main \
  -params-file params.yaml \
  -c my_config \
  -profile singularity,xanadu \
```

## Pipeline output
All of the output from the programs run in the pipeline pipeline will be located in the out directory specified in params.yaml. The pipeline produces the following labeled directories depending on configurations:
* short read qc  
* long read qc  
* assembly  
* polish  
* purge  
* scaffold  
* assembly qc
* visualization
* output  
* pipeline info

Information about interpreting output is located in the [output](docs/output.md) section.

## Credits
emilytrybulec/genomeassembly was originally written by Emily Trybulec.

I thank the following people for their extensive assistance in the development of this pipeline:

University of Connecticut:  
<img align="right" height="210" src="https://github.com/emilytrybulec/argonaut/assets/114685119/7a3fd47c-0fbf-443c-a121-9fd8a3da9ba3">

* Biodiversity and Conservation Genomics Center  
     * Jill Wegrzyn  
     * Cynthia Webster  
     * Anthony He  
     * Laurel Humphrey  
     * Keertana Chagari  
     * Amanda Mueller  
     * Cristopher Guzman  
     * Harshita Akella
  
<img align="right" height="140" src="https://github.com/emilytrybulec/argonaut/assets/114685119/91c25e9f-f70b-481f-8aab-55d2d529eca4">

* Rachel O'Neill Lab  
     * Rachel O’Neill  
     * Michelle Neitzey  
     * Nicole Pauloski  
     * Vel Johnston
  
<img align="right" height="150" src="https://github.com/emilytrybulec/argonaut/assets/114685119/161c0c34-4f05-496d-9436-2d087ba5ccd1">  

* Computational Biology Core  
     * Noah Reid  
     * Gabe Barrett  

* nf-core Community  

* Zbigniew Trybulec

 
## Contributions and Support

Development of this pipeline was funded by the University of Connecticut Office of Undergraduate Research through the Summer Undergraduate Research Fund (SURF) Grant.

<img align="right" height="110" src="https://github.com/emilytrybulec/argonaut/assets/114685119/925bab1e-ba82-44d3-b640-3d7cf6c2028f">  
<img align="right" height="110" src="https://github.com/emilytrybulec/argonaut/assets/114685119/3276ed4c-e272-4ff5-93ba-4a7802bd78aa">  

The Biodiversity and Conservation Genomics Center is a part of the [Earth Biogenome Project](https://www.earthbiogenome.org/), working towards capturing the genetic diversity of life on Earth.

## Citations
Argonaut is currently unpublished. For now, please use the GitHub URL when referencing.

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
