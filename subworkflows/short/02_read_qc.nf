include { FASTP } from '../../modules/nf-core/fastp/main'
include { FASTQC } from '../../modules/nf-core/fastqc/main'
include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { GENOMESCOPE2 } from '../../modules/nf-core/genomescope2/main'
include { SMUDGEPLOT } from '../../modules/local/smudgeplot'

workflow READ_QC {

    take:
        shortreads  // channel: [ val(meta), [ reads ] ]
        ch_db 
           
    main:
    
    ch_versions = Channel.empty()

        FASTP(reads)

        KRAKEN2_KRAKEN2()
        shortreads_filt = KRAKEN2_KRAKEN2.out.fastq_unmapped

        

    emit:
        fastq_filt         // channel: [ val(meta), [ decontaminated fastq (and decontam+length filtered fastq) ] ]
        nanoplot_reads_out   = NANOPLOT.out.html
        centrifuge_out       = CENTRIFUGE_KREPORT.out.kreport
        
    versions = ch_versions                     // channel: [ versions.yml ]
}