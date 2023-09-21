process ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::minimap2=2.24 bioconda::samtools=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' }"

    input:
    tuple val(meta), path(contigs)
    tuple val(meta), path(subreads)

    output:           
    tuple val(meta), path("${subreads.baseName}.aligned.bam") , emit: aligned

    """
    samtools fastq $subreads | \
        minimap2 -t ${task.cpus} -ax map-pb --secondary=no $contigs - \
        | samtools sort -m 10G -o ${contigs.baseName}.aligned.bam
    samtools index ${contigs.baseName}.aligned.bam
    """
}

process HISTOGRAM {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_haplotigs:1.1.2--hdfd78af_0' :
        'biocontainers/purge_haplotigs:1.1.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(contigs)
    tuple val(meta), path(aligned_bam)

    output:
    tuple val(meta), path ("*aligned.bam.gencov"), emit: gencov
    tuple val(meta), path ("*aligned.bam.histogram.png"), emit: hist

    """
    purge_haplotigs hist -t $task.cpus -b $aligned_bam -g $contigs
    """
}

process PURGE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_haplotigs:1.1.2--hdfd78af_0' :
        'biocontainers/purge_haplotigs:1.1.2--hdfd78af_0' }"
        
    input:
    val low
    val mid
    val high
    tuple val(meta), path(assembly)
    path gencov

    output:
    tuple val(meta), path("curated.artefacts.fasta"), emit: purged

    """
    purge_haplotigs cov -in $gencov \
        -low $low -mid $mid -high $high
    purge_haplotigs purge -t $task.cpus -g $assembly -c coverage_stats.csv
    """
}