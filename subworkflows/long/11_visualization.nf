include { GZIP } from '../../modules/local/gzip' 
include { BLOBTOOLS_ADD } from '../../modules/local/blobtools/blobtools_add' 
include { BLOBTOOLS_CONFIG_1LINEAGE } from '../../modules/local/blobtools/blobtools_config' 
include { BLOBTOOLS_CREATE } from '../../modules/local/blobtools/blobtools_create' 
include { BLOBTOOLS_PIPELINE } from '../../modules/local/blobtools/blobtools_pipeline' 
include { BLOBTOOLS_VIEW } from '../../modules/local/blobtools/blobtools_view' 

workflow VISUALIZE {

    take:
        assemblies // channel: [ val(meta), path(assemblies) ]
        lr_fastq // channel: [ val(meta), path(filtered reads) ]
        sr_fastq

    main:
        ch_versions = Channel.empty() 

        GZIP(assemblies)

		BLOBTOOLS_CONFIG_1LINEAGE(GZIP.out.gzip, lr_fastq)
		blobtools_config=BLOBTOOLS_CONFIG_1LINEAGE.out.config
		
		BLOBTOOLS_PIPELINE(blobtools_config, GZIP.out.gzip)
	    BLOBTOOLS_CREATE(assemblies, blobtools_config)
	    BLOBTOOLS_ADD(BLOBTOOLS_PIPELINE.out.blast_out, BLOBTOOLS_PIPELINE.out.diamond_proteome_out, BLOBTOOLS_PIPELINE.out.diamond_busco_out, BLOBTOOLS_PIPELINE.out.assembly_minimap_bam, BLOBTOOLS_PIPELINE.out.hic_minimap_bam, BLOBTOOLS_PIPELINE.out.lineage1_full_table_tsv, BLOBTOOLS_CREATE.out.blobtools_folder)
	    BLOBTOOLS_VIEW(BLOBTOOLS_ADD.out.blobtools_folder)

	    ch_versions = ch_versions.mix(BLOBTOOLS_PIPELINE.out.versions)
		
    emit:
        BLOBTOOLS_VIEW.out.png
        
    versions = ch_versions          
}
