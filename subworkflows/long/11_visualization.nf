include { BLOBTOOLS_RUN } from '../../modules/local/blobtools/blobtools_run' 
include { BLOBTOOLS_CONFIG } from '../../modules/local/blobtools/blobtools_config' 

workflow VISUALIZE {

    take:
        assemblies // channel: [ val(meta), path(assemblies) ]
        ont_fastq // channel: [ val(meta), path(filtered reads) ]
        pb_fastq
        sr_fastq
        bam

    main:
        ch_versions = Channel.empty() 
        println "visualizing assemblies with blobtools!"

		BLOBTOOLS_CONFIG(assemblies, ont_fastq, pb_fastq, sr_fastq)
		blobtools_config=BLOBTOOLS_CONFIG.out.config
		
		BLOBTOOLS_RUN(assemblies, ont_fastq, pb_fastq, sr_fastq, blobtools_config, bam)
	    
	    ch_versions = ch_versions.mix(BLOBTOOLS_RUN.out.versions)
		
    emit:
        BLOBTOOLS_RUN.out.png
        
    versions = ch_versions          
}
