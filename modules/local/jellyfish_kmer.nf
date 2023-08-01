process JELLYFISH_KMER {
    tag "$meta.id"
    label 'process_high'

    container 'ohdzagenetics/jellyfish'

    input:
    tuple val(meta), path(shortreads)

    output:
    tuple val(meta), path("*.jf")   , emit: shortkmer

    script:
    def kmernum = '21'
    """
    cd /
    tar -xvf jellyfish.tar.gz
    
    jellyfish count -t 32 -C -m $kmernum -s 100G -o {$kmernum}mer_out.jf {$shortreads}
    """
}