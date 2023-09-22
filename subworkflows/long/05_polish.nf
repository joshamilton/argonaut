include { MEDAKA } from '../../modules/nf-core/medaka/main'  

workflow POLISH {

    take:
        flye_assembly
        fastq_filt
        model
    main:

    ch_versions = Channel.empty() 

        MEDAKA (fastq_filt, flye_assembly, model)
        medaka_assembly = MEDAKA.out.assembly

        medaka_assembly
                .map { file -> tuple([id: file.baseName], file)  }
                .set { flye_assembly_polished }
    emit:
    
        flye_assembly_polished      
        
    versions = ch_versions                     // channel: [ versions.yml ]
}