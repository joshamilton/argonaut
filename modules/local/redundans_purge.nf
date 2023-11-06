process REDUNDANS_P {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/redundans:2.01--py310pl5321h43eeafb_0' :
        'biocontainers/redundans:2.01--py310pl5321h43eeafb_0' }"


    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(shortreads)

    output:
    tuple val(meta), path("*filled.fa")               , emit: assembly_fasta

    script:
    """
    redundans.py -v -f $fasta -i $shortreads -t $task.cpus --nogapclosing --noscaffolding -o .
    """
}
