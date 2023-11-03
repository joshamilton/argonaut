process TOTAL_BASES_LR {
    label 'process_low'
    tag "$meta.id"

    input:
    tuple val(meta), path(nanoplot_report)

    output:
    path("*totalBasesLR.txt")        , emit: total_bases

    script: 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sed -n '9p' < $nanoplot_report | awk '{print "\""\$3"\""}' | sed -e 's/^"//' -e 's/"\$//' | sed 's/,//g' > ${prefix}_totalBasesLR.txt
    """
}