[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/genomeassembly/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/genomeassembly)


## Introduction

**GenomeAssembly** is a bioinformatics pipeline that performs de novo genome assembly on long read data. A fastq file and input information is fed to the pipeline, resulting in a final assembly with quality checking at each step. [detailed output]

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->

1. Read QC([`Nanoplot`],[`Centrifuge`],[`Bioawk`],[`Kmerfreq`])
2. Assembly([`Flye`])
3. Assembly QC([`BUSCO`],[`Quast`],[`Minimap2`],[`PycoQC`])
4. Polishing([`Medaka`])
5. Polishing QC([`BUSCO`],[`Quast`],[`Minimap2`],[`PycoQC`])
6. Purge([`PurgeDups`])
7. Final QC([`BUSCO`],[`Quast`],[`PycoQC`])



If short reads are inputted as well, 

8. Short Read QC([`Fastp`],[`Kraken`],[`FastQC`],[`Genomescope`],[`Smudgeplot`])
9. Hybrid Assembly([`MaSuRCA`])
10. Alignment([`BWA`])
11. Polishing([`POLCA`])


<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run genomeassembly \
   -params-file genomeassembly/params.yaml \
   -profile <docker/singularity/test/.../institute> \
   --outdir <OUTDIR>
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;

## Pipeline output


## Credits

nf-core/genomeassembly was originally written by emilytrybulec.

We thank the following people for their extensive assistance in the development of this pipeline:

> University of Connecticut:

> Biodiversity and Conservation Genomics Center

     > Jill Wegrzyn

     > Cynthia Webster

     > Rachel O’Neill

     > Michelle Neitzey

> Computational Biology Core

     > Noah Reid

     > Gabe Barrett

> nf-core Community

> Zbigniew Trybulec


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/genomeassembly for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
