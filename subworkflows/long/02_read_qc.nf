include { CENTRIFUGE_CENTRIFUGE } from '../../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_FILTER } from '../../modules/local/centrifuge/filter/main'
include { CENTRIFUGE_KREPORT } from '../../modules/nf-core/centrifuge/kreport/main'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'
include { BIOAWK } from '../../modules/nf-core/bioawk/main'

workflow READ_QC {

    take:
        reads  // channel: [ val(meta), [ reads ] ]
        ch_db
             
    main:
    
    ch_versions = Channel.empty()

        NANOPLOT(reads)

        // if a centrifuge database is provided, run centrifuge and filter out all classified results
        if( ch_db ){
             CENTRIFUGE_CENTRIFUGE        ( reads, ch_db, params.save_unaligned, params.save_aligned, params.sam_format )
             CENTRIFUGE_KREPORT           ( CENTRIFUGE_CENTRIFUGE.out.results, ch_db )
        }

        if(params.min_readlength > 0){
            BIOAWK (reads)
        }
    
    emit:
        CENTRIFUGE_CENTRIFUGE.out.fastq_unmapped            // channel: [ [ val(meta), decontam fastq ] ]
        nanoplot_reads_out   = NANOPLOT.out.html
        centrifuge_out       = CENTRIFUGE_KREPORT.out.kreport
        
    versions = ch_versions                     // channel: [ versions.yml ]
}


