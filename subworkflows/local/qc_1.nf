include { QUAST } from '../../modules/nf-core/quast/main'  
include { BUSCO } from '../../modules/nf-core/busco/main' 
include { PYCOQC } from '../../modules/nf-core/pycoqc/main'  
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main' 
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'   

workflow QC_1 {

    take:
    
        assembly // channel: [ val(meta), [ bam ] ]
        ch_fastq_out

    main:

    ch_versions = Channel.empty() 

        // build index
        MINIMAP2_INDEX(assembly)
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        ch_index = MINIMAP2_INDEX.out.index

        // match index to the corresponding reads
        ch_mapping = ch_index.join(reads)

        // align reads to index
        MINIMAP2_ALIGN(ch_mapping)
        ch_align_bam = MINIMAP2_ALIGN.out.sorted_indexed_bam
        
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
    
        PYCOQC(assembly)
        ch_versions = ch_versions.mix(PYCOQC.out.versions)
  

    emit:

        ch_index
        ch_align_bam
        ch_quast
        ch_busco
        
    versions = ch_versions                     // channel: [ versions.yml ]
}