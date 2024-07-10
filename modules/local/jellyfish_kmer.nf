process JELLYFISH_KMER {
    tag "$meta.id"
    label 'process_high_memory'

    container 'biocontainers/jellyfish:2.2.6--0'

    input:
    tuple val(meta), path(shortreads)
    val kmernum

    output:
    tuple val(meta), path("*.jf")   , emit: shortkmer

    script:
    """
    jellyfish count -t 32 -C -m $kmernum -s 10000000 -o ${kmernum}_mer_out.jf $shortreads
    """
}