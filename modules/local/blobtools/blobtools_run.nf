process BLOBTOOLS_RUN {
    tag "$meta.id"
    label 'process_high'

    container 'genomehubs/blobtoolkit:latest'

    input:
    tuple val(meta), path(assembly), path(busco_full_table_tsv)
    tuple val(meta), path(config)
    tuple val(meta), path(bam)
    path taxon_taxid
    path taxon_taxdump

    output:
    tuple val(meta), path('*/*.json'), emit: json
    tuple val(meta), path('db_*') , emit: db
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def taxid = taxon_taxid ? "--taxid ${taxon_taxid}" : ''
    def taxdump = taxon_taxdump ? "--taxdump ${taxon_taxdump}" : ''
    """
    blobtools create \\
        --fasta $assembly \\
        --meta $config \\
        $taxid \\
        $taxdump \\
        db_${prefix}


    blobtools add \\
        --busco $busco_full_table_tsv \\
        db_${prefix}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtools: \$(blobtools --version | sed -e "s/blobtoolkit v//g")
    END_VERSIONS
    """
}
