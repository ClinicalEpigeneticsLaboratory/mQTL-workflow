#!/usr/bin/env nextflow

// Workflow Params
// CLI params - experiment specific options
params.methylation_idats_dir = 'EPIC/'
params.results_dir = 'results/'
params.sample_sheet = 'sample_sheet.csv'
params.array_position = 'Sentrix_Info'
params.sample_name = 'Sample_Name'

log.info """\
            mQTL v1.0 Worklow two [Methylation arrays]
            ==============
            
            Config:
            ==============
            Number of CPUs: ${params.CPUs}

            Input:
            ==============
            Methylation idats: ${params.methylation_idats_dir}
            Sample sheet: ${params.sample_sheet}
            Array position column: ${params.array_position}
            Sample name column: ${params.sample_name}

            Output:
            ==============
            results directory: ${params.results_dir}

""".stripIndent()

idats = file( params.methylation_idats_dir, checkIfExists: true )
sample_sheet = file( params.sample_sheet, checkIfExists: true )

process prepareIDATs {
    input:
    path idats

    output: 
    path prepared_idats_directory, type: 'dir'
    
    script:
    """
    
    prepare_idats.py $idats ${sample_sheet} ${params.array_position} prepared_idats_directory

    """
}

workflow {

    methylation_idats_dir_ch = prepareIDATs( idats )
}
