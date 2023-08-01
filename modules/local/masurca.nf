process MASURCA {
    tag "$meta.id"
    label 'process_high'

    container 'staphb/masurca:4.1.0'

    input:
    tuple val(meta), path(longreads) //path_to/longreads.gz
    tuple val(meta), path(shortreads) //path_to/pe_R1.fa,/path_to/pe_R2.fa

    output:
    tuple val(meta), path("*.fasta")   , emit: fasta
    path "versions.yml"                , emit: versions

    script:
    def VERSION = '4.1.0'
    
    """
    /MaSuRCA-4.1.0/bin/masurca -t 32 -i $shortreads -r $longreads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MaSuRCA: $VERSION
    END_VERSIONS
    """
}