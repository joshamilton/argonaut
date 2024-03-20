process BLOBTOOLS_VIEW {
    tag "$meta.id"
    label 'process_medium'

    container 'genomehubs/blobtk'

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path('*.png') , emit: png

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """        
    blobtk plot -v snail -d $db -o ${prefix}_snail.png 

    blobtk plot -v cumulative -d $db -o ${prefix}_cumulative.png

    blobtk plot -v blob -d $db -o ${prefix}_blob.png 

    """
}
