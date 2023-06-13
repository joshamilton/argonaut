include { MEDAKA } from '../../modules/nf-core/medaka/main'  

workflow READ_QC {

    take:
            reads
            assembly

    main:

    ch_versions = Channel.empty() 

            MEDAKA (reads, assembly)


    emit:
    
        flye_assembly_polished      = MEDAKA.out.assembly           
        quast_filtered_polished_out   = QUAST.out.html
        
    versions = ch_versions                     // channel: [ versions.yml ]
}