process POLCA {
tag "$meta.id"
    label 'process_medium'

    container 'staphb/masurca:4.1.0'

    input:
    tuple val(meta), path(longreads) //path_to/longreads.gz
    tuple val(meta), path(shortreads) //path_to/pe_R1.fa,/path_to/pe_R2.fa

    output:
    tuple val(meta), path("*.vcf")     , emit: vcf
    tuple val(meta), path("*.report")  , emit: report
    path "versions.yml"                , emit: versions

    script:
    def VERSION = '4.1.0'
    """
    /MaSuRCA-4.1.0/polca.sh -t 6 -a $longreads -r $shortreads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        POLCA/MaSuRCA: $VERSION
    END_VERSIONS
    """
}