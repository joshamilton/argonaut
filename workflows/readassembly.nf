/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowGenomeassembly.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.centrifuge_db ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.shortinput) { ch_shortinput = file(params.shortinput) }
if (params.centrifuge_db) { ch_db = file(params.centrifuge_db) } else { exit 1, 'Centrifuge database not specified!' }
//if (params.summary_txt) {ch_sequencing_summary = file(params.sequencing_summary) } else { ch_sequencing_summary = []}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
//ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
//ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// MODULES
include { OUTPUT } from '../modules/local/output' 

// SUBWORKFLOWS
include { INPUT_CHECK } from '../subworkflows/long/01_input_check'
include { READ_QC } from '../subworkflows/long/02a_read_qc'
include { LENGTH_FILT } from '../subworkflows/long/02b_length_filter'
include { ASSEMBLY } from '../subworkflows/long/03_assembly'
include { QC_1 } from '../subworkflows/long/04_qc_1'
include { POLISH } from '../subworkflows/long/05_polish'
include { QC_2 } from '../subworkflows/long/06_qc_2'
include { HAPS } from '../subworkflows/long/07_purge'
include { QC_3 } from '../subworkflows/long/08_qc_3'
include { SCAFFOLD } from '../subworkflows/long/09_scaffold'
include { QC_4 } from '../subworkflows/long/10_qc_4'

include { INPUT_CHECK2 } from '../subworkflows/short/01_input_check'
include { READ_QC2 } from '../subworkflows/short/02_read_qc'
include { POLISH2 } from '../subworkflows/short/03_polish'

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MASURCA } from '../modules/local/masurca'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMEASSEMBLY {

    ch_versions = Channel.empty()
    
    ch_data = INPUT_CHECK ( ch_input )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    ch_centrifuge_db = Channel.fromPath(params.centrifuge_db)
    
    //decontamination and quality checking of long reads
    READ_QC (ch_data.reads, ch_centrifuge_db)
    ch_versions = ch_versions.mix(READ_QC.out.versions)

    if (params.manual_genome_size){
        genome_size = params.manual_genome_size
    } else {
        genome_size = READ_QC.out[3]
    }
    
    LENGTH_FILT (READ_QC.out[0])
    ch_versions = ch_versions.mix(LENGTH_FILT.out.versions)

    
    if ( params.shortread == true ) {
        ch_shortdata = INPUT_CHECK2 ( ch_shortinput )
    ch_versions = ch_versions.mix(INPUT_CHECK2.out.versions)

        //adaptor trimming and decontamination of short reads if available
        ch_kraken_db = Channel.fromPath(params.kraken_db)
        READ_QC2 (ch_shortdata.reads, ch_kraken_db)
    ch_versions = ch_versions.mix(READ_QC2.out.versions)

        //assembly inputting long & short reads
        ASSEMBLY (LENGTH_FILT.out[0], READ_QC2.out[0], genome_size)
        all_assemblies   = ASSEMBLY.out[0]
    ch_versions = ch_versions.mix(ASSEMBLY.out.versions)
    }
    else {
        ch_shortdata = Channel.empty() 
        //assembly of decontam and length filtered (if specified) long reads
        ASSEMBLY (LENGTH_FILT.out[0], [], genome_size)
        all_assemblies   = ASSEMBLY.out[0]
    ch_versions = ch_versions.mix(ASSEMBLY.out.versions)   
    }

    //short read only assembly
    if ( params.shortread == true && params.longread == false){
        MASURCA (LENGTH_FILT.out[0], READ_QC2.out[0])
        MASURCA.out[0]
            .concat(all_assemblies)
            .set { all_assemblies }
    } 
    
   ch_summtxt = Channel.fromPath(params.summary_txt)

    if ( params.shortread == true ) {
        QC_1 (all_assemblies, LENGTH_FILT.out[0], ch_summtxt, READ_QC2.out[0], READ_QC.out[4])
    ch_versions = ch_versions.mix(QC_1.out.versions)
    } else {
        QC_1 (all_assemblies, LENGTH_FILT.out[0], ch_summtxt, [], READ_QC.out[4])
    ch_versions = ch_versions.mix(QC_1.out.versions)
    }

    //polish flye assembly with medaka
    if ( params.flye == true ) {
        POLISH (ASSEMBLY.out[3], ASSEMBLY.out[1], params.model)
        lr_polish   = POLISH.out[0]

        lr_polish
                .map { file -> tuple([id: file.baseName], file)  }
                .set { medaka_polish }      
        
    ch_versions = ch_versions.mix(POLISH.out.versions)
    }

    //align flye assembly to short reads and polish with POLCA if short reads are available
    if ( params.flye == true && params.shortread == true) {
        POLISH2 (ASSEMBLY.out[3], READ_QC2.out[1])
        sr_polish   = POLISH2.out[0]

        sr_polish
                .map { file -> tuple([id: file.baseName], file)  }
                .set { polca_polish }   
        
        //combine polished flye assemblies w other assemblies
        polca_polish
            .concat(all_assemblies, medaka_polish)
            .set{ polished_assemblies }

    ch_versions = ch_versions.mix(POLISH2.out.versions)
    } else {
        medaka_polish
            .concat(all_assemblies)
            .set{ polished_assemblies }
    }

    polished_assemblies.view()

    if ( params.shortread == true ) {
        QC_2 (polished_assemblies, ASSEMBLY.out[1], ch_summtxt, QC_1.out[3], QC_1.out[4], QC_1.out[5], READ_QC2.out[0], QC_1.out[2], READ_QC.out[4], QC_1.out[7])
    ch_versions = ch_versions.mix(QC_2.out.versions)
    } else {
        QC_2 (polished_assemblies, ASSEMBLY.out[1], ch_summtxt, QC_1.out[3], QC_1.out[4], QC_1.out[5], [], QC_1.out[2], READ_QC.out[4], QC_1.out[7])
    ch_versions = ch_versions.mix(QC_2.out.versions)
    }

    HAPS (polished_assemblies, LENGTH_FILT.out[0])
    ch_versions = ch_versions.mix(HAPS.out.versions)

    if ( params.shortread == true ) {
        QC_3 (HAPS.out[0], ASSEMBLY.out[1], ch_summtxt, QC_2.out[3], QC_2.out[4], QC_2.out[5], READ_QC2.out[0], READ_QC.out[4])
    ch_versions = ch_versions.mix(QC_3.out.versions)
    } else {
        QC_3 (HAPS.out[0], ASSEMBLY.out[1], ch_summtxt, QC_2.out[3], QC_2.out[4], QC_2.out[5], [], READ_QC.out[4])  
    }
        
    if ( params.ragtag_scaffold == true ) {
    ch_reference = Channel.fromPath(params.ragtag_reference)
        SCAFFOLD (HAPS.out[0], ch_reference)
    ch_versions = ch_versions.mix(SCAFFOLD.out.versions)
    }

    if (params.ragtag_scaffold == true && params.shortread == true) {
        QC_4 (SCAFFOLD.out[0], ASSEMBLY.out[1], ch_summtxt, QC_3.out[1], QC_3.out[2], QC_3.out[3], READ_QC2.out[0], READ_QC.out[4]) 
    } else if (params.ragtag_scaffold == true && params.shortread != true){
        QC_4 (SCAFFOLD.out[0], ASSEMBLY.out[1], ch_summtxt, QC_3.out[1], QC_3.out[2], QC_3.out[3], [], READ_QC.out[4])  
    }

    if ( params.ragtag_scaffold == true ) {
        ch_quast = QC_4.out[1]
        ch_busco = QC_4.out[2]
        ch_merqury = QC_4.out[3]
    } else {
        ch_quast = QC_3.out[1]
        ch_busco = QC_3.out[2]
        ch_merqury = QC_3.out[3]
    }

    OUTPUT (ch_quast, ch_busco, ch_merqury)

    assembly_stats  =   OUTPUT.out.assemblyStats

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
   // workflow_summary    = WorkflowGenomeassembly.paramsSummaryMultiqc(workflow, summary_params)
   // ch_workflow_summary = Channel.value(workflow_summary)

    //methods_description    = WorkflowGenomeassembly.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    //ch_methods_description = Channel.value(methods_description)

    //ch_multiqc_files = Channel.empty()
    //ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    //ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
   // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
   // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

 //   MULTIQC (
 //       ch_multiqc_files.collect(),
  //      ch_multiqc_config.toList(),
 //       ch_multiqc_custom_config.toList(),
  //      ch_multiqc_logo.toList()
 //  )
 //   multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//workflow.onComplete {
 //   if (params.email || params.email_on_fail) {
 //      NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
  //  }
 //   NfcoreTemplate.summary(workflow, params, log)
 //   if (params.hook_url) {
 //       NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
  //  }
//}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/