include { RAGTAG } from '../../modules/nf-core/purgedups/pbcstat/main' 

workflow SCAFFOLD {

    take:
        assemblies // channel: val(meta), path(polished assembly alignment to reads)
        reference

    main:

    ch_versions = Channel.empty()

        //optional scaffolding with the same species or most closely related species available
        RAGTAG (assemblies, reference)

    emit:
    
        assembly_polished_purged_scaffolded      = RAGTAG.out.scaffolded_assembly           
        
    versions = ch_versions                     // channel: [ versions.yml ]
}