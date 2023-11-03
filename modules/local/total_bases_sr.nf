process TOTAL_BASES_SR {
    label 'process_low'
    tag "$meta.id"
    
    input:
    tuple val(meta), path(fastp_report)

    output:
    path("totalBasesSR_before*")        , emit: total_bases_before
    path("totalBasesSR_after*")        , emit: total_bases_after

    script: 
    """
    sed -n '7p' < $fastp_report | grep -o -E "[0-9]+" > totalBasesSR_before.txt
    sed -n '18p' < $fastp_report | grep -o -E "[0-9]+" > totalBasesSR_after.txt
    """
}