process RACON {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/racon:1.4.20--h9a82719_1' :
        'biocontainers/racon:1.4.20--h9a82719_1' }"

    input:
    tuple val(meta), path(assembly), path(paf), path(reads)

    output:
    path('*_assembly_racon.fasta') , emit: improved_assembly
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    racon -t "$task.cpus" \\
        "${reads}" \\
        $auto_hybrid_mode \\
        "${paf}" \\
        $args \\
        "${assembly}" > \\
        ${prefix}_assembly_racon.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        racon: \$( racon --version 2>&1 | sed 's/^.*v//' )
    END_VERSIONS
    """
}
