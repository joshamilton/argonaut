include { BLOBTOOLS_RUN } from '../../modules/local/blobtools/blobtools_run' 
include { BLOBTOOLS_CONFIG } from '../../modules/local/blobtools/blobtools_config' 
include { BLOBTOOLS_VIEW } from '../../modules/local/blobtools/blobtools_view' 
include { BLOBTOOLS_BLAST } from '../../modules/local/blobtools/blobtools_blast' 

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

        if (params.PacBioHifi_lr == true){
            assemblies
                .concat(pb_fastq)
                .set{assembly_pb}
            if (params.ONT_lr == true){
                ont_fastq.view()
                 assemblies
                    .concat(ont_fastq)
                    .set{assembly_ont}
                if (params.shortread == true){
                    assemblies
                        .concat(sr_fastq)
                        .set{assembly_sr}
                    BLOBTOOLS_CONFIG(assembly_ont, assembly_pb, assembly_sr)
                } else if (params.shortread == false){
                    BLOBTOOLS_CONFIG(assembly_ont, assembly_pb, [])
                }
            } else if (params.ONT_lr == false && params.shortread == false){
                BLOBTOOLS_CONFIG([], assembly_pb, [])
            } else if (params.ONT_lr == false && params.shortread == true){
                assemblies
                    .concat(sr_fastq)
                    .set{assembly_sr}
                BLOBTOOLS_CONFIG([], assembly_pb, assembly_sr)
            }
        } else if (params.PacBioHifi_lr == false){
            if (params.ONT_lr == true){
                assemblies
                    .concat(ont_fastq)
                    .set{assembly_ont}
                if (params.shortread == true){
                    assemblies
                        .concat(sr_fastq)
                        .set{assembly_sr}
                    BLOBTOOLS_CONFIG(assembly_ont, [], assembly_sr)
                } else if (params.shortread == false){
                    BLOBTOOLS_CONFIG(assembly_ont, [], [])
                }
            } else if (params.ONT_lr == false && params.shortread == true) {
                assemblies
                    .concat(sr_fastq)
                    .set{assembly_sr}
                BLOBTOOLS_CONFIG([], [], assembly_sr)
            }
        }

		blobtools_config=BLOBTOOLS_CONFIG.out.config
		
        assemblies
            .join(busco_table)
            .set{assembly_busco_combo}

        assembly_busco_combo
            .join(bam)
            .set{assembly_busco_bam_combo}

        if (params.taxon_taxid && params.taxon_taxdump){
            if (params.blast_db){
                BLOBTOOLS_BLAST(assemblies)
                BLOBTOOLS_RUN(assembly_busco_bam_combo, blobtools_config, BLOBTOOLS_BLAST.out.hits, params.taxon_taxid, params.taxon_taxdump)}
            else{
                BLOBTOOLS_RUN(assembly_busco_bam_combo, blobtools_config, [], params.taxon_taxid, params.taxon_taxdump)}
        } else {
            if (params.blast_db){
                BLOBTOOLS_BLAST(assemblies)
                BLOBTOOLS_RUN(assembly_busco_bam_combo, blobtools_config, BLOBTOOLS_BLAST.out.hits, [], [])}
            else{
                BLOBTOOLS_RUN(assembly_busco_bam_combo, blobtools_config, [], [], [])
            }
        }

	    ch_versions = ch_versions.mix(BLOBTOOLS_RUN.out.versions)
		
        BLOBTOOLS_VIEW(BLOBTOOLS_RUN.out.db)

    emit:
        BLOBTOOLS_VIEW.out.png
        
    versions = ch_versions          
}
