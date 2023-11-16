[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/genomeassembly/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/genomeassembly)

## Introduction

<img align="right" height="250" src="https://github.com/emilytrybulec/genomeassembly/assets/114685119/9b900dab-44cb-479e-9362-0c0d9dc00ae0">

**Argonaut** is a bioinformatics pipeline that performs de novo genome assembly on long and short read data. A fastq file and input information is fed to the pipeline, resulting in final assemblies with quality checking at each step. The pipeline accepts short reads, long reads, or both. [detailed output]

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->


## Pipeline Summary

Long Read Assembly
1. Read QC and Contaminant Filtering([`Nanoplot`](https://github.com/wdecoster/NanoPlot),[`KmerFreq`](https://github.com/fanagislab/kmerfreq), [`GCE`](https://github.com/fanagislab/GCE), [`Centrifuge`](https://ccb.jhu.edu/software/centrifuge/), [`Recentrifuge`](https://github.com/khyox/recentrifuge))
2. Length Filtering (optional)([`Bioawk`](https://github.com/lh3/bioawk), [`Nanoplot`](https://github.com/wdecoster/NanoPlot))
3. Assembly([`Flye`](https://github.com/fenderglass/Flye), [`Canu`](https://github.com/marbl/canu))
4. Assembly QC([`BUSCO`](https://busco.ezlab.org/),[`Quast`](https://quast.sourceforge.net/),[`Minimap2`](https://github.com/lh3/minimap2),[`PycoQC`](https://github.com/a-slide/pycoQC),[`Merqury`](https://github.com/marbl/merqury))
5. Polish([`Medaka`](https://github.com/nanoporetech/medaka))
6. Polish QC([`BUSCO`](https://busco.ezlab.org/),[`Quast`](https://quast.sourceforge.net/),[`Minimap2`](https://github.com/lh3/minimap2),[`PycoQC`](https://github.com/a-slide/pycoQC),[`Merqury`](https://github.com/marbl/merqury))
7. Purge([`PurgeHaplotigs`](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/))
8. Purge QC([`BUSCO`](https://busco.ezlab.org/),[`Quast`](https://quast.sourceforge.net/),[`Minimap2`](https://github.com/lh3/minimap2),[`PycoQC`](https://github.com/a-slide/pycoQC),[`Merqury`](https://github.com/marbl/merqury))
9. Scaffolding (optional)([`RagTag`](https://github.com/malonge/RagTag))
10. Final QC([`BUSCO`](https://busco.ezlab.org/),[`Quast`](https://quast.sourceforge.net/),[`Minimap2`](https://github.com/lh3/minimap2),[`PycoQC`](https://github.com/a-slide/pycoQC),[`Merqury`](https://github.com/marbl/merqury))

Short Read and Hybrid Assembly
1. Read QC, Contaminant Filtering, Adaptor Trimming([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/),[`GenomeScope`](http://qb.cshl.edu/genomescope/),[`Jellyfish`](https://github.com/gmarcais/Jellyfish),[`Kraken2`](https://ccb.jhu.edu/software/kraken2/), [`Recentrifuge`](https://github.com/khyox/recentrifuge),[`FastP`](https://github.com/OpenGene/fastp))
2. Assembly([`MaSuRCA`](https://github.com/alekseyzimin/masurca), [`Redundans`](https://github.com/Gabaldonlab/redundans))
3. Polish([`POLCA`](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007981))
4. Purge([`Redundans`](https://github.com/Gabaldonlab/redundans))

Hybrid assembly is conducted within the long read assembly subworkflow, and short read assembly is conducted within the main workflow. Downstream quality checking of short read and hybrid assemblies is also conducted in the long read QC subworkflows.

First, prepare a samplesheet with your long read input data that looks as follows:

`longread_samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,single_end
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,,TRUE
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

If short reads are available, prepare a second samplesheet with your short read input data that looks as follows:

`shortread_samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,single_end
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,FALSE
```

Next, create a params.yaml file to specify the paths to your samplesheet(s), contaminant databases, etc. Most likely, a config file will also need to be made to modify the default settings of the pipeline. Please look through the [nextflow.config](nextflow.config) file to browse the defaults and specify which you would like to change in you my_config file. More information is located in the [usage](docs/usage.md) section.

Now, you can run the pipeline using:


```bash
nextflow run emilytrybulec/argonaut \
  -r main \
  -params-file params.yaml \
  -c my_config \
  -profile singularity,xanadu \
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;

## Pipeline output
All of the output from the programs run in the pipeline pipeline will be located in the out directory specified in params.yaml. Information about interpreting output is located in the [output](docs/output.md) section.

## Credits

emilytrybulec/genomeassembly was originally written by Emily Trybulec.

We thank the following people for their extensive assistance in the development of this pipeline:

University of Connecticut:  
  
* Biodiversity and Conservation Genomics Center  
     * Jill Wegrzyn  
     * Cynthia Webster  
     * Rachel O’Neill  
     * Michelle Neitzey  
  
* Computational Biology Core  
     * Noah Reid  
     * Gabe Barrett  

* nf-core Community  
  
* Zbigniew Trybulec


## Contributions and Support

Development of this pipeline was funded by the University of Connecticut Office of Undergraduate Research through the Summer Undergraduate Research Fund (SURF) Grant.

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations
Argonaut is currently unpublished. For now, please use the GitHub URL when referencing.

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
