include { FLYE } from '../../modules/nf-core/flye/main' 
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'

workflow ASSEMBLY {

    take:
        reads           // channel: [ val(meta), [ fastq ] ]
        ch_fastq_out 

    main:

    ch_versions = Channel.empty() 

        NANOPLOT (reads)
        FLYE(reads, ch_fastq_out, params.mode)

        //add option to add manual assembly as input (params.manual_assembly)


    emit:
    // TODO nf-core: edit emitted channels
        flye_assembly      = FLYE.out.fasta           
        nanoplot_filtered_out   = NANOPLOT.out.html
        
    versions = ch_versions                     // channel: [ versions.yml ]
}