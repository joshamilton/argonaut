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
def checkPathParamList = [ params.longinput, params.centrifuge_db ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.longinput) { ch_longinput = file(params.longinput) } else { exit 1, 'Input samplesheet not specified!' }
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

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules

// MODULES
//include { CENTRIFUGE_FILTER } from '../modules/local/centrifuge/filter/main'
// SUBWORKFLOWS
include { INPUT_CHECK } from '../subworkflows/long/01_input_check'
include { READ_QC } from '../subworkflows/long/02a_read_qc'
include { LENGTH_FILT } from '../subworkflows/long/02b_length_filter'
include { ASSEMBLY } from '../subworkflows/long/03_assembly'
include { QC_1 } from '../subworkflows/long/04_qc_1'
include { POLISH } from '../subworkflows/long/05_polish'
//include { QC_2 } from '../subworkflows/long/06_qc_2'
//include { PURGE } from '../subworkflows/long/07_purge'
//include { QC_3 } from '../subworkflows/long/08_qc_3'

include { INPUT_CHECK2 } from '../subworkflows/short/01_input_check'
include { READ_QC2 } from '../subworkflows/short/02_read_qc'
include { ALIGN } from '../subworkflows/short/03_align'
//include { POLISH2 } from '../subworkflows/short/04_polish'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

//include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
//def multiqc_report = []

workflow GENOMEASSEMBLY {

    ch_versions = Channel.empty()
    
    ch_data = INPUT_CHECK ( ch_longinput )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    ch_centrifuge_db = Channel.fromPath(params.centrifuge_db)

    //decontamination and quality checking of long reads
    READ_QC (
        ch_data.reads, ch_centrifuge_db
    )
    ch_versions = ch_versions.mix(READ_QC.out.versions)

    LENGTH_FILT (
        READ_QC.out[0]
    )

    //adaptor trimming and decontamination of short reads if available
    if ( params.shortread == true ) {

        ch_shortdata = INPUT_CHECK2 ( ch_shortinput )
        ch_versions = ch_versions.mix(INPUT_CHECK2.out.versions)

        ch_kraken_db = Channel.fromPath(params.kraken_db)

        READ_QC2 (ch_shortdata.reads, ch_kraken_db)
        ch_versions = ch_versions.mix(READ_QC2.out.versions)

        //assembly inputting everything + shortreads
        ASSEMBLY (
        LENGTH_FILT.out[0], READ_QC2.out[1], READ_QC.out[3]
        )
    ch_versions = ch_versions.mix(ASSEMBLY.out.versions)
    }
    else {
        ch_shortdata = Channel.empty() 

        //assembly of decontam fastq and length filtered fastq (if specified)
        ASSEMBLY (
        LENGTH_FILT.out[0], [], READ_QC.out[3]
        )
    }
    
   ch_summtxt = Channel.fromPath(params.summary_txt)

    if ( params.shortread == true ) {

        QC_1 (
            ASSEMBLY.out[0], ASSEMBLY.out[1], ch_summtxt, READ_QC2.out[0]
        )
        ch_versions = ch_versions.mix(QC_1.out.versions)
    }
    else {
        QC_1 (
            ASSEMBLY.out[0], ASSEMBLY.out[1], ch_summtxt, []
        )
        ch_versions = ch_versions.mix(QC_1.out.versions)
    }

    ALIGN(
        ASSEMBLY.out[0]
    )
    
    POLISH (
        ASSEMBLY.out[0], ASSEMBLY.out[1]
    )

    //QC_2 (
    //    POLISH.out[0], ASSEMBLY.out[1], QC_1.out[0]
    //)

   // PURGE (
   //     POLISH.out[0], ASSEMBLY.out[1]
   // )

   // if ( params.shortread == true ) {
   //     POLISH2 (
   //     )



    //
    // MODULE: Run FastQC
    //
   // FASTQC (
   //     INPUT_CHECK.out.reads
   // )
   // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

   // CUSTOM_DUMPSOFTWAREVERSIONS (
   //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
   // )

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