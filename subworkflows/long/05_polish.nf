include { MEDAKA } from '../../modules/nf-core/medaka/main'  
include { RACON } from '../../modules/nf-core/racon/main'  

workflow POLISH {

    take:
        assembly
        fastq_filt
        model
        paf
    main:

    ch_versions = Channel.empty() 

        RACON (fastq_filt, assembly, paf)

        MEDAKA (fastq_filt, RACON.out.improved_assembly, model)
        ch_versions = ch_versions.mix(MEDAKA.out.versions)
 
        polished_assembly = MEDAKA.out.assembly

        polished_assembly
                .map { file -> tuple(id: file.baseName, file)  }
                .set { assembly_polished }
    emit:
    
        assembly_polished      
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
