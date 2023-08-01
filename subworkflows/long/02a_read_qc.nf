include { CENTRIFUGE_CENTRIFUGE } from '../../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT } from '../../modules/nf-core/centrifuge/kreport/main'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'
include { KMER_FREQ } from '../../modules/local/kmerfreq'
include { GCE } from '../../modules/local/gce'
include { EXTRACT } from '../../modules/local/extract_genome_size'


workflow READ_QC {

    take:
  
        reads  // channel: [ val(meta), [ reads ] ]
        ch_db 
           
    main:
    
    ch_versions = Channel.empty()

        NANOPLOT(reads)

        KMER_FREQ(reads)

        GCE(KMER_FREQ.out.kmerstat, KMER_FREQ.out.kmernum)

        EXTRACT(GCE.out.gce2log)

        // if a centrifuge database is provided, run centrifuge and filter out all classified results
        if( ch_db ){
             CENTRIFUGE_CENTRIFUGE        ( reads, ch_db, params.save_unaligned, params.save_aligned, params.sam_format )
             CENTRIFUGE_KREPORT           ( CENTRIFUGE_CENTRIFUGE.out.results, ch_db )
        }

        //fastq_filt
         //   .map { file -> tuple(file.baseName, file)  }
          //  .set { filtered_fastq }

        //filtered_fastq.view()

    emit:
        fastq_filt           = CENTRIFUGE_CENTRIFUGE.out.fastq_unmapped // channel: [ val(meta), path(decontaminated fastq) ]
        nanoplot_reads_out   = NANOPLOT.out.html
        centrifuge_out       = CENTRIFUGE_KREPORT.out.kreport
        est_genome_size      = EXTRACT.out.genome_size_est
        
    versions = ch_versions                     // channel: [ versions.yml ]
}