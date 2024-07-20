#!/usr/bin/env nextflow

// Workflow Params
params.sample_sheet = file('sample_sheet.csv')
params.methylation_idats_dir = file('EPIC/')
params.results_dir = file('results/')

params.array_position = 'Sentrix_Info'
params.sample_name = 'Sample_Name'

log.info """\
            mQTL $version Worklow two [Methylation arrays]
            ==============
            
            Config:
            ==============
            Cell fraction correction [--correction <boolean: true/false>]: ${params.correction} [default: false]
            Deconvolution method [--deconvolution_method <str: CP/RPC/CBS>]: ${params.deconvolution_method} [default: RPC]
            Collapse methylation readings to the cg prefix [--collapse_prefix <boolean: true/false>]: ${params.collapse_prefix} [default: true]
            Number of CPUs [--CPUs <int>]: ${params.CPUs} [default: 10]

            Input:
            ==============
            Methylation idats [--methylation_idats_dir <path>]: ${params.methylation_idats_dir}
            Sample sheet [--sample_sheet <path>]: ${params.sample_sheet}
            Array position column [--array_position <str>]: ${params.array_position}
            Sample name column [--sample_name <str>]: ${params.sample_name}

            Output:
            ==============
            results directory [--results_dir <path>]: ${params.results_dir}

""".stripIndent()

process validateParams {
    script:
    """
    #!/usr/bin/python3
    import os
    import multiprocessing
    import pandas as pd 

    # --CPUs flag check
    cpus_available = multiprocessing.cpu_count() - 1
    assert ${params.CPUs} <= cpus_available, f"Exceeded number of available cpus --> {${params.CPUs}} / {cpus_available}" 

    # Files flags check
    files = ["${params.sample_sheet}", "${params.methylation_idats_dir}"]
    for file in files:
        assert os.path.exists(file), f"File does not exists {file}"

    # --array_position and --sample_name flags check
    sample_sheet_fields = pd.read_csv("${params.sample_sheet}").columns
    assert "${params.array_position}" in sample_sheet_fields, f"${params.array_position} not present in sample sheet columns: {sample_sheet_fields}"
    assert "${params.sample_name}" in sample_sheet_fields, f"${params.sample_name} not present in sample sheet columns: {sample_sheet_fields}"

    # Additional flags check
    assert "${params.correction}" in ["true", "false"], f"--correction should be either true or false, is ${params.correction}"
    assert "${params.collapse_prefix}" in ["true", "false"], f"--collapse_prefix should be either true or false, is ${params.collapse_prefix}"
    assert "${params.deconvolution_method}" in ["RPC", "CP", "CBS"], f"--deconvolution_method should be CP, RPC or CBS is ${params.deconvolution_method}"

    """
}

process prepareIDATs {
    input:
    path idats

    output: 
    path 'prepared_idats_directory', type: 'dir'
    
    script:
    """
    
    prepare_idats.py $idats ${params.sample_sheet} ${params.array_position} prepared_idats_directory

    """
}

process processIDATs {
    errorStrategy 'retry'
    maxRetries 5
    // From time to time MulticoreParam fails to start 
    // It is not elegant solution but in most cases restarting solves the problem

    input:
    path idats

    output:
    path 'mynorm.parquet'

    script:
    """
    
    sesame.R $idats ${params.CPUs} ${params.collapse_prefix}

    """
}

process renameSamples {
    publishDir "$params.results_dir/flow_two", mode: 'copy', overwrite: true, pattern: 'mynorm.parquet'

    input:
    path mynorm

    output:
    path 'mynorm.parquet'

    script:
    """

    rename_samples.py $mynorm ${params.sample_sheet} ${params.array_position} ${params.sample_name}

    """
}

process CellFractionCorrection {
    publishDir "$params.results_dir/flow_two", mode: 'copy', overwrite: true, pattern: 'mynorm_corrected.parquet'
    publishDir "$params.results_dir/flow_two", mode: 'copy', overwrite: true, pattern: '*.csv'

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
    validateParams()

    prepared_idats = prepareIDATs( params.methylation_idats_dir )
    mynorm = processIDATs( prepared_idats )
    mynorm = renameSamples( mynorm )

    if ( params.correction == true ) {
        mynorm = CellFractionCorrection( mynorm )
    }
}
