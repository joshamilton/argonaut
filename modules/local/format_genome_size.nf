process FORMAT {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04':
        'ubuntu:20.04'}"
        
    input:
    val genome_size

    output:
    path("shortenedSizeFinal.txt")        , emit: genome_size_est
    path("standardSize.txt")     , emit: standard_fmt_est    

    script: 
    """
    echo $genome_size > standardSize.txt 
    
    head -1 standardSize.txt | numfmt --to=si > shortenedSizeFinal.txt
    """
}