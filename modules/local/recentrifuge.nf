process RECENTRIFUGE {
tag "$meta.id"
    label 'process_medium'

    container 'replikation/recentrifuge'

    input:
    tuple val(meta), path(reads) //path_to/longreads.gz

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