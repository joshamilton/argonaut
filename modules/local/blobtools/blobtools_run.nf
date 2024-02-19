process BLOBTOOLS_RUN {
    tag "$meta.id"
    label 'process_medium'

    container 'chrishah/blobtools'

    input:
    tuple val(meta), path(assembly)
    tuple val(meta), path(ont_fastq)
    tuple val(meta), path(pacbio_fastq)
    tuple val(meta), path(illumina_fastq)
    tuple val(meta), path(config)
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.png'), emit: png
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    blobtools create \\
        -i $assembly \\
        -b $bam \\
        --meta $config \\
        --taxdump ./blobtoolkit/taxdump


    blobtools view -i *.json

    blobtools plot -i *.json
    """
}
