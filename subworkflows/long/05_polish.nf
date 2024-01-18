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

        if (params.racon_polish == true) {
            fastq_filt
            .join(assembly, by:0)
            .concat(paf)
            .view()
            .set { ch_racon }

            RACON (fastq_filt, assembly, paf)
        ch_versions = ch_versions.mix(RACON.out.versions)
        }
        
        if (params.medaka_polish == true && params.racon_polish == true) {
            MEDAKA (fastq_filt, RACON.out.improved_assembly, model)
        ch_versions = ch_versions.mix(MEDAKA.out.versions)
 
        polished_assembly = MEDAKA.out.assembly
        } else {
            MEDAKA (fastq_filt, assembly, model)
        ch_versions = ch_versions.mix(MEDAKA.out.versions)
 
        polished_assembly = MEDAKA.out.assembly
        }
        

        polished_assembly
                .map { file -> tuple(id: file.baseName, file)  }
                .set { assembly_polished }
    emit:
    
        assembly_polished      
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
