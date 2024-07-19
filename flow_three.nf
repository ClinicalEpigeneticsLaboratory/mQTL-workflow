#!/usr/bin/env nextflow

// Workflow Params
// CLI params - experiment specific options
params.results_dir = 'results/'
params.sample_sheet = 'sample_sheet.csv'
params.sample_name = 'Sample_Name'

// CLI default params
params.correction = true
params.deconvolution_method = 'RPC'

// Input files: mynorm as well as genotype table should be available from params.results_dir  

genotype_table = file("$params.results_dir/genotype_table.parquet", checkIfExists: true)

if ( file("$params.results_dir/mynorm_corrected.parquet").exists() )  {
    mynorm = file( "$params.results_dir/mynorm_corrected.parquet" )
}
else if ( file("$params.results_dir/mynorm.parquet").exists() ){
    mynorm = file( "$params.results_dir/mynorm.parquet" )
}
else {
    error("Mynorm file does not exists in ${params.results_dir} directory!")
}


log.info """\
            mQTL $version Worklow three [mQTL]
            ==============
            
            Config:
            ==============
            Number of CPUs: ${params.CPUs}

            Input:
            ==============
            Mynorm: $mynorm
            Genotype table: $genotype_table

            Output:
            ==============
            results directory: ${params.results_dir}

""".stripIndent()

workflow {

}
