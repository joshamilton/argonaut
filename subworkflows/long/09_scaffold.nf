include { RAGTAG } from '../../modules/nf-core/purgedups/pbcstat/main' 

workflow SCAFFOLD {

    take:
        assemblies // channel: val(meta), path(polished assembly alignment to reads)

    main:

    ch_versions = Channel.empty()

        RAGTAG (assemblies)

    emit:
    
        assembly_polished_purged_scaffolded      = RAGTAG.out.scaffolded_assembly           
        
    versions = ch_versions                     // channel: [ versions.yml ]
}