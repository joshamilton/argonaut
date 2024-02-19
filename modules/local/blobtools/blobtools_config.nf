process BLOBTOOLS_CONFIG {
    tag "$meta.id"
    label 'process_medium'

    container 'chrishah/blobtools'

    input:
    tuple val(meta), path(assembly)
    tuple val(meta), path(ont_fastq)
    tuple val(meta), path(pacbio_fastq)
    tuple val(meta), path(illumina_fastq)

    output:
    tuple val(meta), path('*config.yaml'), emit: config
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (params.shortread == true){
    """
    echo "assembly:
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
        paired:
          - prefix: ${meta.id}
            platform: ILLUMINA
            file: $illumina_fastq
    version: 1" > ${prefix}_config.yaml  


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
        """
    }

    if (params.longread == true && params.ONT_lr == true){
    """
    echo "assembly:
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
        paired:
          - prefix: ${meta.id}
            platform: ONT
            file: $ont_fastq
    version: 1" > ${prefix}_config.yaml  


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
    }

    if (params.longread == true && params.PacBioHifi_lr == true){
    """
    echo "assembly:
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
        paired:
          - prefix: ${meta.id}
            platform: PACBIO
            file: $pacbio_fastq
    version: 1" > ${prefix}_config.yaml  


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtools: \$(blobtools --version)
    END_VERSIONS
    """
    }
}
