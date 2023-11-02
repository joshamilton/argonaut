process TOTAL_BASES_SR {
    label 'process_low'

    input:
    tuple val(meta), path(fastp_report)

    output:
    path("totalBasesSR.txt")        , emit: total_bases

    script: 
    """
    sed -n '7p;18p' < $fastpreport | grep -o -E "[0-9]+" > totalBasesSR.txt
    """
}