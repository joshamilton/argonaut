include { PBBAM_PBMERGE } from '../../modules/local/pb_bam'
include { BAM2FASTX } from '../../modules/local/bam2fastx'
include { CUTADAPT } from '../../modules/local/cutadapt'
include { MERYL_COUNT } from '../../modules/nf-core/meryl/count/main'
include { MERYL_HISTOGRAM } from '../../modules/nf-core/meryl/histogram'
include { GENOMESCOPE2 } from '../../modules/nf-core/genomescope2/main'
include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { RECENTRIFUGE_KR } from '../../modules/local/recentrifuge/kraken'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'

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
	    CUTADAPT (input_pacbio)
	    ch_versions = ch_versions.mix(CUTADAPT.out.versions)

        MERYL_COUNT (CUTADAPT.out.reads, params.kmer_num)

        MERYL_COUNT.out.meryl_db
                .map { file -> tuple(id: file.baseName, file)  }
                .set { hifi_meryldb }

        MERYL_HISTOGRAM (hifi_meryldb)

        GENOMESCOPE2(MERYL_HISTOGRAM.out.hist)

        //decontamination of trimmed short reads
        KRAKEN2_KRAKEN2(input_pacbio, ch_db, params.save_output_fastqs, params.save_reads_assignment)

        //summarizing and visualizing decontam
        RECENTRIFUGE_KR(KRAKEN2_KRAKEN2.out.classified_reads_assignment, params.rcf_db)

        filt_pbhifi = KRAKEN2_KRAKEN2.out.unclassified_reads_fastq   

    emit:
        filt_pbhifi    // channel: [ val(meta), path(decontaminated fastq) ]
        nanoplot_reads_out   = NANOPLOT.out.html
        nanoplot_report_txt  = NANOPLOT.out.txt
        genome_size_est = GENOMESCOPE2.out.summary
        
    versions = ch_versions                     // channel: [ versions.yml ]
}