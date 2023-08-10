include { QUAST } from '../../modules/local/quast'  
include { BUSCO } from '../../modules/nf-core/busco/main' 
include { PYCOQC } from '../../modules/nf-core/pycoqc/main'  
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main' 
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'   
include { OUTPUT } from '../../modules/local/output' 

workflow QC_3 {

    take:
    
        flye_assembly_polished // channel: [ val(meta), [ bam ] ]
        fastq_filt // channel: [ val(meta), path(filtered reads) ]
        ch_summarytxt // channel from params.summarytxt
        ch_quast
        ch_busco
        ch_merqury

    main:

    ch_versions = Channel.empty()     
        // build index
        MINIMAP2_INDEX(flye_assembly_polished)
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        ch_index = MINIMAP2_INDEX.out.index

        // match index to the corresponding reads
        //ch_mapping = ch_index.join(reads)

        // align reads to index
        MINIMAP2_ALIGN(fastq_filt, [], params.bam_format, params.cigar_paf_format, params.cigar_bam)
        ch_align2_bam = MINIMAP2_ALIGN.out.bam
        ch_align2_paf = MINIMAP2_ALIGN.out.paf
        
        //run quast
        QUAST(
            FLYE.out.fasta.map{it -> it[1]}.collect() // this has to be aggregated because of how QUAST makes the output directory for reporting stats
        )
        ch_quast = QUAST.out
        ch_versions = ch_versions.mix(QUAST.out.versions)

        //run BUSCO
        BUSCO(assembly, params.busco_db)
        ch_busco = BUSCO.out
        ch_versions = ch_versions.mix(BUSCO.out.versions)
    
        PYCOQC(ch_summarytxt)
        ch_versions = ch_versions.mix(PYCOQC.out.versions)

        OUTPUT (ch_quast, ch_busco, ch_merqury)

    emit:
        assembly_stats  =   OUTPUT.out.assemblyStats
        
    versions = ch_versions                     // channel: [ versions.yml ]
}