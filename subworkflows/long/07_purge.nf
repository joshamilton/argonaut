include { ALIGN } from '../../modules/local/purgehap' 
include { HISTOGRAM } from '../../modules/local/purgehap' 
include { PURGE } from '../../modules/local/purgehap' 

workflow HAPS {

    take:
        assembly
        reads
        
    main:
    ch_versions = Channel.empty()

        if (params.low == null && params.mid == null && params.high == null) {
        println "generating assembly histogram with purge haplotigs!"

        ALIGN(reads, assembly, params.bam_format, params.cigar_paf_format, params.cigar_bam)
        HISTOGRAM(assembly, ALIGN.out.bam)
        assemblies_polished_purged      = Channel.empty()

        } else if (params.low != null && params.mid != null && params.high != null){
        println "purging assemblies with purge haplotigs!"

        PURGE(params.low, params.mid, params.high, assembly, params.gencov)
        purged_assemblies      = PURGE.out.purged           
        
        purged_assemblies
                .map { file -> tuple(id: file.baseName, file)  }
                .set { assemblies_polished_purged }
        }

    emit:
        assemblies_polished_purged
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
