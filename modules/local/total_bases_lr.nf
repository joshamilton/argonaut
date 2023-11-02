process TOTAL_BASES_LR {
    label 'process_low'

    input:
    tuple val(meta), path(nanoplot_report)

    output:
    path("totalBasesLR.txt")        , emit: total_bases

    script: 
    """
    sed -n '9p' < $nanoplot_report | awk '{print "\""$3"\""}' | sed -e 's/^"//' -e 's/"$//' > totalBasesLR.txt
    """
}