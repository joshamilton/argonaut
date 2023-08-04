include { QUAST } from '../../modules/local/quast'  
include { BUSCO } from '../../modules/nf-core/busco/main' 
include { PYCOQC } from '../../modules/nf-core/pycoqc/main'  
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main' 
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'  
include { MERYL_COUNT } from '../../modules/nf-core/meryl/count/main' 
include { MERQURY } from '../../modules/nf-core/merqury/main' 
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main' 


workflow QC_2 {

    take:
        polished_flye_assembly // channel: [ val(meta), path(flye assembly.fasta) ]
        fastq_filt // channel: [ val(meta), path(filtered reads) ]
        summarytxt // channel from params.summarytxt

    main:

    ch_versions = Channel.empty() 

    if ( params.flye == true ) {
        // build index
        MINIMAP2_INDEX(polished_flye_assembly)
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        ch_index = MINIMAP2_INDEX.out.index

        // align reads
        MINIMAP2_ALIGN(fastq_filt, polished_flye_assembly.map{it[1]}, params.bam_format, params.cigar_paf_format, params.cigar_bam)
        ch_align_bam = MINIMAP2_ALIGN.out.bam
        ch_align_paf = MINIMAP2_ALIGN.out.paf
        
        // run quast
        QUAST(
            polished_flye_assembly.map{it -> it[1]}.collect() // this has to be aggregated because of how QUAST makes the output directory for reporting stats
        )
        ch_quast = QUAST.out.results
        ch_versions = ch_versions.mix(QUAST.out.versions)

        // run BUSCO
        BUSCO(polished_flye_assembly, params.busco_lineage, [], [])
        ch_busco = BUSCO.out.batch_summary
        ch_versions = ch_versions.mix(BUSCO.out.versions)

        SAMTOOLS_INDEX (MINIMAP2_ALIGN.out.bam)
    
        // create summary txt channel with meta id and run pycoQC
        ch_summarytxt = summarytxt.map { file -> tuple(file.baseName, file) }
        
        PYCOQC (
            ch_summarytxt, MINIMAP2_ALIGN.out.bam, SAMTOOLS_INDEX.out.bai
        )
        ch_versions = ch_versions.mix(PYCOQC.out.versions)

        MERYL_COUNT ( fastq_filt )

        MERQURY (
            polished_flye_assembly, MERYL_COUNT.out.meryl_db
        )
    }

    emit:
        ch_polished_index
        ch_polished_align_bam
        ch_polished_align_paf
        ch_polished_quast
        ch_polished_busco
        
    versions = ch_versions                     // channel: [ versions.yml ]
}