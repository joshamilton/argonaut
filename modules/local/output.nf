process OUTPUT {
    label 'process_low'

    input:
    tuple val(meta), path(ch_quast_tsv)
    tuple val(meta), path(ch_busco)
    tuple val(meta), path(ch_merqury)

    output:
    tuple val(meta), path("*.assemblyStats.txt")       , emit: assemblyStats
    tuple val(meta), path("all_assemblies.tsv")       , emit: all_assemblyStats
   
    script: 
    def prefix
    """
    prefix=\$(awk 'NR==1 {print \$2}' $ch_quast_tsv)
    echo -ne "quast output\n" >> \$prefix.assemblyStats.txt
    less $ch_quast_tsv >> \$prefix.assemblyStats.txt

    echo -ne "\n busco score\n" >> \$prefix.assemblyStats.txt
    grep -A 17 'Results:' $ch_busco >> \$prefix.assemblyStats.txt

    echo -ne "merqury quality score\n" >> \$prefix.assemblyStats.txt
    awk '{ print \$4 }' $ch_merqury >> \$prefix.assemblyStats.txt

    total_summ=\$(cat \$prefix.assemblyStats.txt)
    
    awk 'BEGIN{ FS = OFS = "\t" } { print $0, "\$prefix" : "\$total_summ" }'  > all_assemblies.tsv
    
    """
}
