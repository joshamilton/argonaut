include { MEDAKA } from '../../modules/nf-core/medaka/main'  

workflow POLISH {

    take:
        flye_assembly
        fastq_filt
    main:

    ch_versions = Channel.empty() 

        ch_assembly_reads = Channel.empty() 
        ch_assembly_reads.concat(fastq_filt, flye_assembly.map{it[1]}).view()

        MEDAKA (ch_assembly_reads)


    emit:
        flye_assembly_polished      = MEDAKA.out.assembly           
        
    versions = ch_versions                     // channel: [ versions.yml ]
}