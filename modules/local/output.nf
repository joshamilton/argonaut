process OUTPUT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(ch_quast_tsv)
    tuple val(meta), path(ch_busco)
    tuple val(meta), path(ch_merqury)

    output:
    tuple val(meta), path("assemblyStats.txt")                 , emit: assemblyStats
   
    script: 
    """
    echo -ne "quast output\n" > assemblyStats.txt
    less $ch_quast_tsv > assemblyStats.txt

    echo -ne "busco score\n" > assemblyStats.txt
    grep -A 17 'Results:' $ch_busco > assemblyStats.txt

    echo -ne "merqury score\n" > assemblyStats.txt
    awk '{ print \$4 }' $ch_merqury > assemblyStats.txt
    """
}