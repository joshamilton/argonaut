include { GZIP } from '../../modules/local/gzip'
include { FASTP } from '../../modules/nf-core/fastp/main'
include { FASTQC } from '../../modules/nf-core/fastqc/main'
include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { FASTQC_2 } from '../../modules/local/fastqc2/main'
include { FASTQC_3 } from '../../modules/local/fastqc3/main'
include { GUNZIP } from '../../modules/nf-core/gunzip/main'
include { GENOMESCOPE2 } from '../../modules/nf-core/genomescope2/main'
include { JELLYFISH_KMER } from '../../modules/local/jellyfish_kmer'
include { JELLYFISH_HIST } from '../../modules/local/jellyfish_hist'
include { RECENTRIFUGE_KR } from '../../modules/local/recentrifuge/kraken'

workflow READ_QC2 {

    take:
        shortreads  // channel: [ val(meta), [ reads ] ]
        ch_db 
           
    main:
    
    ch_versions = Channel.empty()

        //zipping files for downstream processing
        GZIP(shortreads)

        //qc raw short reads
        FASTQC(shortreads)

        //adapter trimming
        FASTP(GZIP.out.gzip, params.adapter_fasta, params.save_trimmed_fail, params.save_merged)

        //qc adapter trimmed short reads
        FASTQC_2(FASTP.out.reads)

        //decontamination of trimmed short reads
        KRAKEN2_KRAKEN2(FASTP.out.reads, ch_db, params.save_output_fastqs, params.save_reads_assignment)

        //summarizing and visualizing decontam
        RECENTRIFUGE_KR(KRAKEN2_KRAKEN2.out.classified_reads_assignment, params.rcf_db)

        //unzip decontam short reads
        GUNZIP(KRAKEN2_KRAKEN2.out.unclassified_reads_fastq)

        //qc decontaminated short reads
        FASTQC_3(KRAKEN2_KRAKEN2.out.unclassified_reads_fastq)

        JELLYFISH_KMER(GUNZIP.out.gunzip, params.kmer_num)
        JELLYFISH_HIST(JELLYFISH_KMER.out.shortkmer)
        
        GENOMESCOPE2(JELLYFISH_HIST.out.shortkmer_hist, params.kmer_num)

    emit:
        filt_shortreads = KRAKEN2_KRAKEN2.out.unclassified_reads_fastq   // channel: [ val(meta), [ decontaminated and adaptor trimmed short reads ] ]
        unzip_filt_shortreads = GUNZIP.out.gunzip

    versions = ch_versions                     // channel: [ versions.yml ]
}