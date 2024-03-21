process BLOBTOOLS_CONFIG {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(assembly) from assemblies
    path ont_fastq
    path pacbio_fastq
    path illumina_fastq

    output:
    tuple val(meta), path('*config.yaml'), emit: config

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def configContent = """
    assembly:
      accession: ${meta.id}
      file: $assembly
      level: scaffold
      prefix: ${meta.id}
    busco:
      download_dir: ${params.busco_lineage}
      lineages:
        - ${params.busco_lineage}
      basal_lineages:
        - ${params.busco_lineage}
    reads:
    """

    if (illumina_fastq) {
        configContent += """
          paired:
            - prefix: ${meta.id}
              platform: ILLUMINA
              file: $illumina_fastq
        """
    }

    if (ont_fastq) {
        configContent += """
          single:
            - prefix: ${meta.id}
              platform: ONT
              file: $ont_fastq
        """
    }

    if (pacbio_fastq) {
        configContent += """
          single:
            - prefix: ${meta.id}
              platform: PACBIO
              file: $pacbio_fastq
        """
    }

    // Write the combined configuration to the output file
    """
    echo '$configContent' > ${prefix}_config.yaml
    """
}
