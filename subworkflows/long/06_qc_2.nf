include { QUAST } from '../../modules/local/quast'  
include { BUSCO } from '../../modules/nf-core/busco/main' 
include { PYCOQC } from '../../modules/nf-core/pycoqc/main'  
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main' 
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'    

workflow QC_2 {

    take:
    
        flye_assembly_polished // channel: [ val(meta), path(polished assembly)]
        fastq_filt // channel: [ val(meta), path(filtered reads) ]
        ch_summarytxt // channel from params.summarytxt

    main:

    ch_versions = Channel.empty()     
        // build index
        MINIMAP2_INDEX(flye_assembly_polished)
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        ch_index = MINIMAP2_INDEX.out.index

        // match index to the corresponding reads
        ch_mapping = ch_index.join(fastq_filt)

        // align reads to index
        MINIMAP2_ALIGN(ch_mapping)
        ch_align_bam = MINIMAP2_ALIGN.out.sorted_indexed_bam
        
        //run quast
        QUAST(
            FLYE.out.fasta.map{it -> it[1]}.collect() // this has to be aggregated because of how QUAST makes the output directory for reporting stats
        )
        ch_quast.concat(QUAST.out.results)
        ch_versions = ch_versions.mix(QUAST.out.versions)

        //run BUSCO
        BUSCO(assembly, params.busco_db)
        ch_busco.concat(BUSCO.out.batch_summary)
        ch_versions = ch_versions.mix(BUSCO.out.versions)
    
        //run pycoQC
       // PYCOQC(ch_summarytxt)
       // ch_versions = ch_versions.mix(PYCOQC.out.versions)

    emit:
        ch_index
        ch_align_bam
        ch_quast
        ch_busco
        
    versions = ch_versions                     // channel: [ versions.yml ]
}