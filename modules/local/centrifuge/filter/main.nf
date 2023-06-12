process CENTRIFUGE_FILTER {
    tag "$meta.id"
    label 'process_high'

    input:
    path reads
    path unclassified_seq

    output:
    path "filtered.fastq", emit: fastq

    script:
    """
    seqtk subseq ${reads} ${unclassified_seq} > filtered.fastq
    """

}