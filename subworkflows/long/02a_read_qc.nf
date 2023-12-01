include { GUNZIP } from '../../modules/nf-core/gunzip/main'
include { CENTRIFUGE_CENTRIFUGE } from '../../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT } from '../../modules/nf-core/centrifuge/kreport/main'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'
include { KMER_FREQ } from '../../modules/local/kmerfreq'
include { GCE } from '../../modules/local/gce'
include { RECENTRIFUGE_C } from '../../modules/local/recentrifuge/centrifuge'
include { SEQKIT_GREP } from '../../modules/local/seqkit/grep/main'
include { TAR } from '../../modules/local/tar'

workflow READ_QC {

    take:
  
        reads  // channel: [ val(meta), [ reads ] ]
        ch_db 
           
    main:
    
    ch_versions = Channel.empty()

        if ( params.tar == true ) { 
            TAR(reads)
            reads = TAR.out.untar }

        NANOPLOT(reads)

        if (params.manual_genome_size == null){
            KMER_FREQ(reads)

            GCE(KMER_FREQ.out.kmerstat, KMER_FREQ.out.kmernum)
            gce_genome_size      = GCE.out.gce2log
        } else{
            gce_genome_size      = Channel.empty() 
        }

        // if a centrifuge database is provided, run centrifuge and filter out all classified results
        if( ch_db ){
             CENTRIFUGE_CENTRIFUGE        ( reads, ch_db, params.save_unaligned, params.save_aligned, params.sam_format )
             CENTRIFUGE_KREPORT           ( CENTRIFUGE_CENTRIFUGE.out.results, ch_db )
        }

        SEQKIT_GREP(CENTRIFUGE_CENTRIFUGE.out.results, reads)

        RECENTRIFUGE_C(CENTRIFUGE_CENTRIFUGE.out.results, params.rcf_db)

        fastq_filt           = SEQKIT_GREP.out.filter

        fastq_filt
            .map { file -> tuple([id:file.baseName, single_end:true], file)  }
            .set { filtered_fastq }

    emit:
        filtered_fastq    // channel: [ val(meta), path(decontaminated fastq) ]
        nanoplot_reads_out   = NANOPLOT.out.html
        centrifuge_out       = CENTRIFUGE_KREPORT.out.kreport
        gce_genome_size
        nanoplot_report_txt  = NANOPLOT.out.txt
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
