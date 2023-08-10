include { QUAST } from '../../modules/local/quast'  
include { BUSCO } from '../../modules/nf-core/busco/main' 
include { PYCOQC } from '../../modules/nf-core/pycoqc/main'  
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main' 
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'  
include { MERYL_COUNT } from '../../modules/nf-core/meryl/count/main' 
include { MERQURY } from '../../modules/nf-core/merqury/main' 
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main' 
include { OUTPUT } from '../../modules/local/output' 

workflow QC_3 {

    take:
        assemblies // channel: [ val(meta), path(assemblies) ]
        fastq_filt // channel: [ val(meta), path(filtered reads) ]
        summarytxt // channel from params.summarytxt
        ch_quast
        ch_busco
        ch_merqury

    main:

    ch_versions = Channel.empty() 

        // build index
        MINIMAP2_INDEX(assemblies)
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        ch_index = MINIMAP2_INDEX.out.index

        // align reads
        MINIMAP2_ALIGN(fastq_filt, assemblies.map{it[1]}, params.bam_format, params.cigar_paf_format, params.cigar_bam)
        ch_align_bam = MINIMAP2_ALIGN.out.bam
        ch_align_paf = MINIMAP2_ALIGN.out.paf
        
        // run quast
        QUAST(
            assemblies.map{it -> it[1]}.collect() // this has to be aggregated because of how QUAST makes the output directory for reporting stats
        )
        ch_quast
            .concat(QUAST.out.results)
            .set { ch_quast }
        ch_versions = ch_versions.mix(QUAST.out.versions)

        // run BUSCO
        BUSCO(assemblies, params.busco_lineage, [], [])
        ch_busco
            .concat(BUSCO.out.short_summaries_txt)
            .set { ch_busco }
        ch_versions = ch_versions.mix(BUSCO.out.versions)

        SAMTOOLS_INDEX (MINIMAP2_ALIGN.out.bam)
    
        // create summary txt channel with meta id and run pycoQC
        ch_summarytxt = summarytxt.map { file -> tuple(file.baseName, file) }
        
        PYCOQC (
            ch_summarytxt, MINIMAP2_ALIGN.out.bam, SAMTOOLS_INDEX.out.bai
        )
        ch_versions = ch_versions.mix(PYCOQC.out.versions)

        if ( params.shortread == true ) {
            MERYL_COUNT ( shortreads ) }
        else {
            MERYL_COUNT ( fastq_filt )
        }

        MERQURY (
            assemblies, MERYL_COUNT.out.meryl_db
        )
        ch_merqury
            .concat(MERQURY.out.assembly_qv)
            .set { ch_merqury }
        ch_versions = ch_versions.mix(MERQURY.out.versions)

        OUTPUT (ch_quast, ch_busco, ch_merqury)

    emit:
        ch_index
        ch_align_bam
        ch_align_paf
        ch_quast
        ch_busco
        ch_merqury
        assembly_stats  =   OUTPUT.out.assemblyStats
        
    versions = ch_versions                     // channel: [ versions.yml ]
}