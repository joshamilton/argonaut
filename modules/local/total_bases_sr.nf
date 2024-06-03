process TOTAL_BASES_SR {
    label 'process_low'
    tag "$meta.id"
    
    input:
    tuple val(meta), path(fastp_report)

    output:
    tuple val(meta), path("totalBasesSR_before.txt")        , emit: total_bases_before
    tuple val(meta), path("totalBasesSR_filtered.txt")        , emit: total_bases_after

    script: 
    """
    sed -n '7p' < $fastp_report | grep -o -E "[0-9]+" > totalBasesSR_before.txt
    sed -n '18p' < $fastp_report | grep -o -E "[0-9]+" > totalBasesSR_filtered.txt

    head -1 totalBasesSR_before.txt | numfmt --grouping > totalBasesSR_before.txt
    head -1 totalBasesSR_filtered.txt | numfmt --grouping > totalBasesSR_filtered.txt
    """
}
