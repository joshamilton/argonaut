include { ALIGN } from '../../modules/local/purgehap' 
include { HISTOGRAM } from '../../modules/local/purgehap' 
include { PURGE } from '../../modules/local/purgehap' 

workflow PURGE {

    take:
        assembly
        subreads
        
    main:

    ch_versions = Channel.empty()

        if ( ! (params.low && params.mid && params.high )) {

        ALIGN(assembly, subreads)
        HISTOGRAM(assembly, ALIGN.out.aligned)
        assemblies_polished_purged      = Channel.empty()

        } else {
        
        PURGE(assembly, HISTOGRAM.out.hist)
        assemblies_polished_purged      = PURGE.out.purged           

        }

    emit:
        assemblies_polished_purged
        
    versions = ch_versions                     // channel: [ versions.yml ]
}