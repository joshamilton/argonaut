process TAR {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    path("*.gz")        , emit: untar

    script: 
    """
    tar -xvzf $reads
    """
}
