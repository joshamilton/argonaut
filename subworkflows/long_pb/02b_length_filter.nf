include { BIOAWK } from '../../modules/nf-core/bioawk/main'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'

workflow LENGTH_FILT3 {

    take:
  
        decontam_reads  // channel: [ val(meta), [ decontam reads ] ] 
           
    main:
    
    ch_versions = Channel.empty()

        // if statement for if min_read_length exists for length filter
        if(params.min_readlength > 0){
            BIOAWK(decontam_reads)
            NANOPLOT(BIOAWK.out.output)

            longreads = BIOAWK.out.output  // channel: [ val(meta), path(decontam+length filtered fastq) ]
        }
        else{
            longreads = decontam_reads  // channel: [ val(meta), path(decontaminated fastq) ]
        }

    emit:
        longreads
        
    versions = ch_versions                     // channel: [ versions.yml ]
}