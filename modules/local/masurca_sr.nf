process MASURCA_SR {
    tag "$meta.id"
    label 'process_high_memory', 'error_ignore'

    container 'staphb/masurca'
    
    input:
    tuple val(meta), path(shortreads) //path_to/pe_R1.fa /path_to/pe_R2.fa

    output:
    path("CA.mr.99.17.15.0.02/masurca*")                , emit: fasta
    path ("CA.mr.99.17.15.0.02/versions.yml")                , emit: versions

    script:
    def VERSION = '4.1.0'
    def sr
    def prefix = task.ext.prefix ?: "masurca_sr_only_${meta.id}"
    """
    sr=\$(echo '${shortreads}' | sed -e "s/ /,/g")
    masurca -t $task.cpus -i \$sr

    cd CA.mr.99.17.15.0.02
    mv primary.genome.scf.fasta ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MaSuRCA: $VERSION
    END_VERSIONS
    """
}