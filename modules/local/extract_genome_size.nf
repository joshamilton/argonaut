process EXTRACT {
    tag "$meta.id"
    label 'process_low'

    container 'emilytrybulec/genassembly:kmer'

    input:
    tuple val(meta), path(gce2log)

    output:
    tuple val(meta), path("shortenedSizeFinal.txt")        , emit: genome_size_est
    tuple val(meta), path("standardSizeFinal.txt")     , emit: standard_fmt_est    

    script: 
    def number
    """
    tail -3 $gce2log | awk '{print \$5}' > scientificSize.txt
    less scientificSize.txt | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print \$1 * (10 ^ \$2)}' > standardSize.txt 
    
    head -1 standardSize.txt > standardSizeFinal.txt
    
    head -1 standardSize.txt | numfmt --to=si > shortenedSizeFinal.txt
    """
}