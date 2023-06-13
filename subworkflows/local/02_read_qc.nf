// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { CENTRIFUGE_CENTRIFUGE } from '../../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_FILTER } from '../../modules/local/centrifuge/filter/main'
include { CENTRIFUGE_KREPORT } from '../../modules/nf-core/centrifuge/kreport/main'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'
include { BIOAWK } from '../../modules/nf-core/bioawk/main'

workflow READ_QC {

    take:
  
        reads  // channel: [ val(meta), [ reads ] ]
        ch_fastq
        db 
           
    main:
    
    ch_versions = Channel.empty()

        ch_db = Channel.fromPath(params.db)

        NANOPLOT(reads)

        // if a centrifuge database is provided, run centrifuge and filter out all classified results
        if( ch_db ){
             CENTRIFUGE_CENTRIFUGE        ( reads, ch_db, params.save_unaligned, params.save_aligned, params.sam_format )
             //CENTRIFUGE_FILTER            ( reads, CENTRIFUGE_CENTRIFUGE.out.unmapped.fastq.gz )
            // ideally we would kill pipeline execution here if percent retained is too low. 
                // not sure how to do that. 
             CENTRIFUGE_KREPORT           ( reads, CENTRIFUGE_CENTRIFUGE.out.report )
        }


    // if a centrifuge database was provided, filtered fastq file, otherwise, original fastq file
    ch_fastq_out = { ch_db ? CENTRIFUGE_CENTRIFUGE.out.unmapped.fastq.gz : ch_fastq }
 
    // add if statement for if min_read_length exists for length filter

    emit:
 
        ch_fastq_out            // channel: [ val(meta), [ decontam fastq ] ]
        nanoplot_reads_out   = NANOPLOT.out.html
        centrifuge_out       = CENTRIFUGE_KREPORT.out.kreport
        
    versions = ch_versions                     // channel: [ versions.yml ]
}

