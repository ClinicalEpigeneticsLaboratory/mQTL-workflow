#!/usr/bin/env nextflow

// Workflow Params
// CLI params - experiment specific options
params.methylation_idats_dir = 'EPIC/'
params.results_dir = 'results/'
params.sample_sheet = 'sample_sheet.csv'
params.array_position = 'Sentrix_Info'
params.sample_name = 'Sample_Name'

// CLI default params
params.correction = true
params.deconvolution_method = 'RPC'


log.info """\
            mQTL v1.0 Worklow two [Methylation arrays]
            ==============
            
            Config:
            ==============
            Number of CPUs: ${params.CPUs}
            Cell fraction correction: ${params.correction}
            Deconvolution method: ${params.deconvolution_method}

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
    path 'prepared_idats_directory', type: 'dir'
    
    script:
    """
    
    prepare_idats.py $idats ${sample_sheet} ${params.array_position} prepared_idats_directory

    """
}

process processIDATs {

    input:
    path idats

    output:
    path 'mynorm.parquet'

    script:
    """
    
    sesame.R $idats ${params.CPUs}

    """
}

process renameSamples {
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: 'mynorm.parquet'

    input:
    path mynorm

    output:
    path 'mynorm.parquet'

    script:
    """

    rename_samples.py $mynorm $sample_sheet ${params.array_position} ${params.sample_name}

    """
}

process CellFractionCorrection {
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: 'mynorm_corrected.parquet'
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: '*.csv'

    input:
    path mynorm

    output:
    path '*.parquet', emit: mynorm_corrected
    path '*.csv'

    script: 
    """

    correction.R $mynorm ${params.deconvolution_method}
    
    """
}

workflow {

    methylation_idats_dir_ch = prepareIDATs( idats )
    mynorm = processIDATs( methylation_idats_dir_ch )
    mynorm = renameSamples( mynorm )

    if ( params.correction == true ) {
        mynorm = CellFractionCorrection( mynorm )
    }
    else if (params.correction == false ) {
        println('Cell fraction correction will not be applied')
    }
    else {
        println('--correction param should be true or false')
    }
}
