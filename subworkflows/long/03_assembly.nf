include { GUNZIP } from '../../modules/nf-core/gunzip/main'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'
include { FLYE } from '../../modules/nf-core/flye/main' 
include { MASURCA } from '../../modules/local/masurca'
include { CANU } from '../../modules/nf-core/canu/main' 

workflow ASSEMBLY {

    take:
        longreads   // channel: [ val(meta), [ decontaminated and/or length filtered fastq ] ]
        shortreads
        genome_size_est

    main:
    ch_versions = Channel.empty() 
    
        //unzips centrifuge fastq output
        GUNZIP (longreads) 

        //makes sure long read input is filtered reads
        NANOPLOT (GUNZIP.out.gunzip)

        //if statement to run assembly and create channels for each resulting assembly
        if ( params.flye == true ) {
            println "assembling with flye!"
            FLYE(GUNZIP.out.gunzip, params.flye_mode, genome_size_est)
            flye_assembly      = FLYE.out.fasta   
        }
        if ( params.canu == true ) {
            println "assembling with canu!"
            CANU(GUNZIP.out.gunzip, params.canu_mode, )
        }
        if ( params.shortread == true ) {
            println "assembling with maSuRCA!"
            MASURCA(GUNZIP.out.gunzip, shortreads)
            //put maSuRCA command and channel here
        }
        if ( params.ex_assembly == true ) {
            println "inputting existing assembly!"
          //  existing_assembly   = channel from something
        }

    emit:
        flye_assembly  
        longreads_unzipped     = GUNZIP.out.gunzip     
        nanoplot_filtered_out   = NANOPLOT.out.html
        
    versions = ch_versions                     // channel: [ versions.yml ]
}