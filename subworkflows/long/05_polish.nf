include { MEDAKA } from '../../modules/nf-core/medaka/main'  

workflow POLISH {

    take:
        flye_assembly
        fastq_filt
        model
    main:

    ch_versions = Channel.empty() 

        MEDAKA (fastq_filt, flye_assembly, model)

    emit:
    
        flye_assembly_polished      = MEDAKA.out.assembly           
        
    versions = ch_versions                     // channel: [ versions.yml ]
}