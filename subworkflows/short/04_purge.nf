include { REDUNDANS_P } from '../../modules/local/redundans_purge'  

workflow PURGE2 {

    take:
        assembly
        shortreads
        
    main:

    ch_versions = Channel.empty() 

        REDUNDANS_P (assembly, shortreads)


    emit:
    
        sr_assembly_polished      = REDUNDANS_P.out.assembly_fasta          
        
    versions = ch_versions                     // channel: [ versions.yml ]
}