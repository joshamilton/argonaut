process CAT {
    label 'process_low'
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"
        
    input:
    tuple val(meta), path(ont), path(hifi)

    output:
    path("*.gz")        , emit: cat_longreads

    script: 

    def prefix = task.ext.prefix ?: "longreads_combined"
    """
    cat "${ont}" "${hifi}"> ${prefix}.fastq.gz
    """
}
