include { BWAMEM2_INDEX } from '../../modules/nf-core/bwamem2/index/main'  
include { BWAMEM2_MEM } from '../../modules/nf-core/bwamem2/mem/main' 

//align short reads against assemblies

workflow ALIGN {

    take:
        flye_assembly // channel: [ val(meta), path(flye assembly.fasta) ]
        fastq_filt // channel: [ val(meta), path(filtered reads) ]
        summarytxt // channel from params.summarytxt

    main:

    ch_versions = Channel.empty() 

    if ( params.flye == true ) {
        // build index
        MINIMAP2_INDEX(flye_assembly)
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        ch_index = MINIMAP2_INDEX.out.index

        // match index to the corresponding reads
       // ch_mapping = ch_index.join(fastq_filt)

        // align reads
        MINIMAP2_ALIGN(fastq_filt, flye_assembly.map{it[1]}, params.bam_format, params.cigar_paf_format, params.cigar_bam)
        ch_align_bam = MINIMAP2_ALIGN.out.bam
        ch_align_paf = MINIMAP2_ALIGN.out.paf
        
        // run quast
        QUAST(
            flye_assembly.map{it -> it[1]}.collect() // this has to be aggregated because of how QUAST makes the output directory for reporting stats
        )
        ch_quast = QUAST.out.results
        ch_versions = ch_versions.mix(QUAST.out.versions)

        // run BUSCO
        BUSCO(flye_assembly, params.busco_lineage, [], [])
        ch_busco = BUSCO.out.batch_summary
        ch_versions = ch_versions.mix(BUSCO.out.versions)
    
        // create summary txt channel with meta id and run pycoQC
        ch_summarytxt = Channel.empty() 
        ch_summarytxt.concat(fastq_filt.map{it[0]}, summarytxt).view()  //channel: val(meta), summarytxt
        PYCOQC(ch_summarytxt)
        ch_versions = ch_versions.mix(PYCOQC.out.versions)
    }

    emit:
        ch_index
        ch_align_bam
        ch_quast
        ch_busco
        
    versions = ch_versions                     // channel: [ versions.yml ]
}