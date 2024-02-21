include { MEDAKA } from '../../modules/nf-core/medaka/main'  
include { RACON } from '../../modules/nf-core/racon/main'  
include { GUNZIP } from '../../modules/nf-core/gunzip/main' 

workflow POLISH {

    take:
        assembly
        fastq_filt
        model
        paf
        ch_racon
    main:

    ch_versions = Channel.empty() 

        if (params.racon_polish == true) {
            println "polishing assemblies with racon!"
            RACON (ch_racon)
            
            no_meta_polished_assembly = RACON.out.improved_assembly

            no_meta_polished_assembly
                .map { file -> tuple(id: file.baseName.replaceAll(/\\.fasta$/, ''), file) }
                .set { racon_polished_assembly }

            GUNZIP (racon_polished_assembly)
            assembly_polished = GUNZIP.out.gunzip

        ch_versions = ch_versions.mix(RACON.out.versions)
        }
        
        if (params.medaka_polish == true && params.racon_polish == true) {
            println "polishing racon-polished assemblies with medaka!"
            MEDAKA (fastq_filt, racon_polished_assembly, model)
        ch_versions = ch_versions.mix(MEDAKA.out.versions)
 
            assembly_polished = MEDAKA.out.assembly
        
            assembly_polished
                .map { file -> tuple(id: file.baseName, file)  }
                .set { medaka_assembly_polished }

            assembly_polished
                .concat(no_meta_polished_assembly)
                .set{no_meta_polished_assembly}
            
        } else if (params.medaka_polish == true && params.racon_polish == false){

            println "polishing assemblies with medaka!"
            MEDAKA (fastq_filt, assembly, model)
        ch_versions = ch_versions.mix(MEDAKA.out.versions)
 
        no_meta_polished_assembly = MEDAKA.out.assembly

        no_meta_polished_assembly
                .map { file -> tuple(id: file.baseName, file)  }
                .set { assembly_polished }
        } 
        
    emit:
    
        no_meta_polished_assembly      
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
