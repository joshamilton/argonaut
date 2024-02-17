include { PBBAM_PBMERGE } from '../../modules/local/pb_bam'
include { BAM2FASTX } from '../../modules/local/bam2fastx'
include { CUTADAPT } from '../../modules/local/cutadapt'
include { KMER_FREQ } from '../../modules/local/kmerfreq'
include { GCE } from '../../modules/local/gce'
include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { RECENTRIFUGE_KR } from '../../modules/local/recentrifuge/kraken'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'
include { TOTAL_BASES_LR } from '../../modules/local/total_bases_lr' 

workflow READ_QC3 {

    take:
  
        input_pacbio  // channel: [ val(meta), [ reads ] ]
        ch_db 
           
    main:
    
    ch_versions = Channel.empty()

        //PBBAM_PBMERGE(input_pacbio)
		//ch_versions = ch_versions.mix(PBBAM_PBMERGE.out.versions)
		//final_pacBio_bam	= PBBAM_PBMERGE.out.bam
		//final_pacBio_bam_index	= PBBAM_PBMERGE.out.pbi

        //BAM2FASTX (final_pacBio_bam.join(final_pacBio_bam_index))
        //ch_versions = ch_versions.mix(BAM2FASTX.out.versions)
        //bam2fastx_output	= BAM2FASTX.out.reads

        NANOPLOT(input_pacbio)
        TOTAL_BASES_LR (NANOPLOT.out.txt)

	    CUTADAPT (input_pacbio)
	    ch_versions = ch_versions.mix(CUTADAPT.out.versions)

        //decontamination of trimmed short reads
        KRAKEN2_KRAKEN2(CUTADAPT.out.reads, ch_db, params.save_output_fastqs, params.save_reads_assignment)

        //summarizing and visualizing decontam
        RECENTRIFUGE_KR(KRAKEN2_KRAKEN2.out.classified_reads_assignment, params.rcf_db)

        filt_pbhifi = KRAKEN2_KRAKEN2.out.unclassified_reads_fastq   

        KMER_FREQ(filt_pbhifi)

        GCE(KMER_FREQ.out.kmerstat, KMER_FREQ.out.kmernum)
        gce_genome_size      = GCE.out.gce2log

    emit:
        filt_pbhifi    // channel: [ val(meta), path(decontaminated fastq) ]
        nanoplot_reads_out   = NANOPLOT.out.html
        nanoplot_report_txt  = NANOPLOT.out.txt
        base_count           = TOTAL_BASES_LR.out.total_bases
        gce_genome_size
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
