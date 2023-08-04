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


## Pipeline Summary

Long Read Assembly
1. Read QC and Contaminant Filtering([`Nanoplot`](https://github.com/wdecoster/NanoPlot),[`Centrifuge`](https://ccb.jhu.edu/software/centrifuge/))
2. Assembly([`Flye`](https://github.com/fenderglass/Flye),[`MaSuRCA`](https://github.com/alekseyzimin/masurca),[`Canu`](https://github.com/marbl/canu))
3. Assembly QC([`BUSCO`](https://busco.ezlab.org/),[`Quast`](https://quast.sourceforge.net/),[`Minimap2`](https://github.com/lh3/minimap2),[`PycoQC`](https://github.com/a-slide/pycoQC),[`Merqury`](https://github.com/marbl/merqury))
4. Polish([`Medaka`](https://github.com/nanoporetech/medaka))
5. Polish QC([`BUSCO`](https://busco.ezlab.org/),[`Quast`](https://quast.sourceforge.net/),[`Minimap2`](https://github.com/lh3/minimap2),[`PycoQC`](https://github.com/a-slide/pycoQC),[`Merqury`](https://github.com/marbl/merqury))
6. Purge([`PurgeDups`](https://github.com/dfguan/purge_dups))
7. Final QC([`BUSCO`](https://busco.ezlab.org/),[`Quast`](https://quast.sourceforge.net/),[`Minimap2`](https://github.com/lh3/minimap2),[`PycoQC`](https://github.com/a-slide/pycoQC),[`Merqury`](https://github.com/marbl/merqury))
8. Scaffolding([`RagTag`](https://github.com/malonge/RagTag))

Short Read Assembly (optional)
1. Read QC, Contaminant Filtering, Adaptor Trimming([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/),[`GenomeScope`](http://qb.cshl.edu/genomescope/),[`Kraken2`](https://ccb.jhu.edu/software/kraken2/),[`FastP`](https://github.com/OpenGene/fastp),[`Jellyfish`](https://github.com/gmarcais/Jellyfish))
2. Align([`BWA`](https://bio-bwa.sourceforge.net/bwa.shtml))
3. Polish([`POLCA`](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007981))


First, prepare a samplesheet with your long read input data that looks as follows:

`longread_samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,single_end
CONTROL_REP1,AEG588A1_S1_L002.fastq.gz,,TRUE
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

If short reads are available, prepare a second samplesheet with your short read input data that looks as follows:

`shortread_samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,single_end
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,FALSE
```

Specify your type of input, assemblers to run, and more in the nextflow.config file. Ensure that the file is correct, especially regarding the boolean statements about input and assembly. Construct a params.yaml file to specify paths to your samplesheets and databases. An example params.yaml is available for your reference.

Now, you can run the pipeline using:

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
<pre>
Biodiversity and Conservation Genomics Center

     Jill Wegrzyn

     Cynthia Webster

     Rachel O’Neill

     Michelle Neitzey

Computational Biology Core

     Noah Reid

     Gabe Barrett
</pre>

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
