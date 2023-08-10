include { PURGEDUPS_PBCSTAT } from '../../modules/nf-core/purgedups/pbcstat/main' 
include { PURGEDUPS_CALCUTS } from '../../modules/nf-core/purgedups/calcuts/main' 
include { PURGEDUPS_PURGEDUPS } from '../../modules/nf-core/purgedups/purgedups/main'


workflow PURGE {

    take:
        paf_alignment // channel: val(meta), path(polished assembly alignment to reads)
        flye_assembly_polished
        fastq_filt
        
    main:

    ch_versions = Channel.empty()

        PURGEDUPS_PBCSTAT (paf_alignment)

        PURGEDUPS_CALCUTS (PURGEDUPS_PBCSTAT.out.stat)

        PURGEDUPS_PURGEDUPS (PURGEDUPS_PBCSTAT.out.basecov, PURGEDUPS_CALCUTS.out.cutoff, paf_alignment)

    emit:
    
        flye_assembly_polished_purged      = PURGEDUPS_PURGEDUPS.out.bed           
        
    versions = ch_versions                     // channel: [ versions.yml ]
}