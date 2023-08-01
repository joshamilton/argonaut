include { BWAMEM2_INDEX } from '../../modules/nf-core/bwamem2/index/main'  
include { BWAMEM2_MEM } from '../../modules/nf-core/bwamem2/mem/main' 

//align short reads against assemblies
workflow ALIGN {

    take:
        filt_shortreads
        flye_assembly // channel: [ val(meta), path(assembly.fasta) ] ideally all assemblies

    main:

    ch_versions = Channel.empty() 

    if ( params.flye == true ) {
        BWAMEM2_INDEX(flye_assembly)
        BWAMEM2_MEM(filt_shortreads, BWAMEM2_INDEX.out.index, params.bwa_sort_bam)
    }

    emit:
        ch_flye_aligned   =   BWAMEM2_MEM.out.bam
        
    versions = ch_versions                     // channel: [ versions.yml ]
}