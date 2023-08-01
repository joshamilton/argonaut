process JELLYFISH_HIST {
    tag "$meta.id"
    label 'process_medium'

    container 'ohdzagenetics/jellyfish'

    input:
    tuple val(meta), path(shortkmer)

    output:
    tuple val(meta), path("*.histo")               , emit: shortkmer_hist

    script:
    def kmernum = '21'
    """
    cd /
    tar -xvf jellyfish.tar.gz
    
    jellyfish histo -o {$kmernum}mer_out.histo {$shortkmer}
    """
}