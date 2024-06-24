process TAR {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04':
        'ubuntu:20.04'}"
        
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.gz")        , emit: untar

    script: 
    """
    tar -xvzf $reads
    """
}
