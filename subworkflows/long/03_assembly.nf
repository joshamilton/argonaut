//include { GUNZIP } from '../../modules/nf-core/gunzip/main'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'
include { FLYE } from '../../modules/nf-core/flye/main' 
include { MASURCA } from '../../modules/local/masurca'
include { CANU } from '../../modules/nf-core/canu/main' 

workflow ASSEMBLY {

    take:
        longreads   // channel: [ val(meta), [ decontaminated and/or length filtered fastq ] ]
        shortreads
        genome_size_est

    main:
    ch_versions = Channel.empty() 

        //makes sure long read input is filtered reads
        NANOPLOT (longreads)

        assemblies = Channel.empty() 

        //if statement to run assembly and create channels for each resulting assembly
        if ( params.flye == true ) {
            println "assembling with flye!"
            FLYE(longreads, params.flye_mode, genome_size_est)
            flye_assembly      = FLYE.out.fasta   

            flye_assembly
                .map { file -> tuple([id: file.baseName], file)  }
                .set { f_assembly }       
        }
        else {
            f_assembly = Channel.empty() 
        }
        if ( params.canu == true ) {
            println "assembling with canu!"
            CANU(longreads, params.canu_mode, genome_size_est)
            canu_assembly      = CANU.out.assembly   

            canu_assembly
                .map { file -> tuple([id: file.baseName], file)  }
                .set { c_assembly }
        }
        else {
            c_assembly = Channel.empty() 
        }
        if ( params.masurca == true ) {
            println "assembling with maSuRCA!"
            MASURCA(longreads, shortreads)
            masurca_assembly    = MASURCA.out.fasta

            masurca_assembly
                .map { file -> tuple([id: file.baseName], file)  }
                .set { m_assembly }
        }
        else {
            m_assembly = Channel.empty() 
        }
        if ( params.ex_assembly == true ) {
            println "inputting existing assembly!"
            existing_assembly = Channel.fromPath(params.existing_assembly)

            existing_assembly
                .map { file -> tuple([id: file.baseName], file)  }
                .set { ex_assembly }
        }
        else {
            ex_assembly = Channel.empty() 
        }

        assemblies
            .concat(f_assembly, c_assembly, m_assembly, ex_assembly)
            .set { all_assemblies }
        
        all_assemblies.view()

    emit:
        all_assemblies  
        longreads   
        nanoplot_filtered_out   = NANOPLOT.out.html
        f_assembly
        
    versions = ch_versions                     // channel: [ versions.yml ]
}