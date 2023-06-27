/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowShortreadassembly.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.shortinput ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.shortinput) { ch_shortinput = file(params.shortinput) } else { exit 1, 'Short read input samplesheet not specified!' }
//if (params.fastq) { ch_fastq = file(params.fastq) } else { exit 1, 'Input reads in fastq format not specified!' }
if (params.db) { ch_db = file(params.db) } else { exit 1, 'Centrifuge database not specified!' }


include { INPUT_CHECK } from '../subworkflows/short/01_input_check'
include { READ_QC } from '../subworkflows/short/02_read_qc'
include { ALIGN } from '../subworkflows/short/03_align'
include { POLISH } from '../subworkflows/short/04_polish'

workflow SHORTREAD {

    ch_versions = Channel.empty()

    ch_data = INPUT_CHECK ( ch_shortinput )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    ch_db = Channel.fromPath(params.kraken_db)

    READ_QC (ch_data.reads, ch_db)

    QC_1()

    POLISH()

}