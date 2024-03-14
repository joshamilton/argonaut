include { BLOBTOOLS_RUN } from '../../modules/local/blobtools/blobtools_run' 
include { BLOBTOOLS_CONFIG } from '../../modules/local/blobtools/blobtools_config' 
include { BLOBTOOLS_VIEW } from '../../modules/local/blobtools/blobtools_view' 

workflow VISUALIZE {

    take:
        assemblies // channel: [ val(meta), path(assemblies) ]
        ont_fastq // channel: [ val(meta), path(filtered reads) ]
        pb_fastq
        sr_fastq
        bam
        busco_table

    main:
        ch_versions = Channel.empty() 

        if (params.PacBioHifi_lr == true && params.ONT_lr == false && params.shortread == false){
            BLOBTOOLS_CONFIG(assemblies, [], pb_fastq, [])
        } else if (params.PacBioHifi_lr == true && params.ONT_lr == true && params.shortread == false){
            BLOBTOOLS_CONFIG(assemblies, ont_fastq, pb_fastq, [])
        } else if (params.PacBioHifi_lr == false && params.ONT_lr == true && params.shortread == false){
            BLOBTOOLS_CONFIG(assemblies, ont_fastq, [], [])
        } else if (params.PacBioHifi_lr == false && params.ONT_lr == true && params.shortread == true){
            BLOBTOOLS_CONFIG(assemblies, ont_fastq, [], sr_fastq)
        } else if (params.PacBioHifi_lr == true && params.ONT_lr == false && params.shortread == true){
            BLOBTOOLS_CONFIG(assemblies, [], pb_fastq, sr_fastq)
        } else if (params.PacBioHifi_lr == false && params.ONT_lr == false && params.shortread == true){
            BLOBTOOLS_CONFIG(assemblies, [], [], sr_fastq)
        } else if (params.PacBioHifi_lr == true && params.ONT_lr == true && params.shortread == true){
            BLOBTOOLS_CONFIG(assemblies, ont_fastq, pb_fastq, sr_fastq)
        } 
		blobtools_config=BLOBTOOLS_CONFIG.out.config
		
        assemblies
            .join(busco_table)
            .set{assembly_busco_combo}

        if (params.taxon_taxid && params.taxon_taxdump){
            BLOBTOOLS_RUN(assembly_busco_combo, blobtools_config, bam, params.taxon_taxid, params.taxon_taxdump)
        } else {
            BLOBTOOLS_RUN(assembly_busco_combo, blobtools_config, bam, [], [])
        }

	    ch_versions = ch_versions.mix(BLOBTOOLS_RUN.out.versions)
		
        BLOBTOOLS_VIEW(BLOBTOOLS_RUN.out.db)

    emit:
        BLOBTOOLS_VIEW.out.png
        
    versions = ch_versions          
}
