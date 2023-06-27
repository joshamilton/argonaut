include { GUNZIP } from '../../modules/nf-core/gunzip/main' 
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'
include { FLYE } from '../../modules/nf-core/flye/main'
include { MASURCA } from '../../modules/local/masurca'

workflow ASSEMBLY {

    take:
        fastq_out           // channel: [ val(meta), decontam fastq ] ] 

    main:
        
    ch_versions = Channel.empty()
        
        GUNZIP (fastq_out) 

        fastq_assembly = { GUNZIP.out.gunzip }

        //makes sure input is filtered reads
        NANOPLOT (fastq_assembly)

        //if statement to run assembly and create channels for each resulting assembly
        if ( params.flye == true ) {
            println "assembling with flye!"
            FLYE(fastq_assembly, params.mode)
            flye_assembly      = FLYE.out.fasta    
        }
        if ( params.shortread == true ) {
            println "assembling with maSuRCA!"
            //put maSuRCA command and channel here
        }
        if ( params.ex_assembly == true ) {
            println "inputting existing assembly!"
          //  existing_assembly   = channel from something
        }

    emit:
        flye_assembly  
        fastq_assembly 
        nanoplot_filtered_out   = NANOPLOT.out.html
        
    versions = ch_versions                     // channel: [ versions.yml ]
}