process TAR {
    tag "$meta.id"
    label 'process_low'

    container 'nf-core/ubuntu:20.04'
        
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.gz")        , emit: untar

    script: 
    """
    tar -xvzf $reads
    """
}
