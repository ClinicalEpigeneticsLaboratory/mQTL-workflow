#!/usr/bin/env nextflow

// Workflow Params
// CLI params - experiment specific options
params.results_dir = 'results/'
params.csv_methylation_manifest = 'resources/infinium-methylationepic-v-1-0-b5-manifest-file.csv'
params.csv_gsa_manifest = 'resources/GSA-24v3-0_A1.csv'

// CLI default params
params.correction = true
params.deconvolution_method = 'RPC'
params.distance = 50000
params.alpha = 0.05
params.slope = 0.05
params.nCPU = 20

// Input files: manifests, mynorm as well as genotype table should be available

csv_methylation_manifest = file("$params.csv_methylation_manifest", checkIfExists: true)
csv_gsa_manifest  = file("$params.csv_gsa_manifest", checkIfExists: true)

genotype_table = file("$params.results_dir/genotype_table.parquet", checkIfExists: true)

if ( file("$params.results_dir/mynorm_corrected.parquet").exists() )  {
    mynorm = file( "$params.results_dir/mynorm_corrected.parquet" )
}
else if ( file("$params.results_dir/mynorm.parquet").exists() ){
    mynorm = file( "$params.results_dir/mynorm.parquet" )
}
else {
    error("mynorm file does not exists in ${params.results_dir} directory!")
}


log.info """\
            mQTL $version Worklow three [mQTL]
            ==============
            
            Config:
            ==============
            Number of CPUs: ${params.CPUs}
            GSA manifest: ${csv_gsa_manifest}
            Methylation manifest: ${csv_methylation_manifest}
            Alpha: ${params.alpha}
            Slope: ${params.slope}
            Distance: ${params.distance}
            
            Input:
            ==============
            Mynorm: $mynorm
            Genotype table: $genotype_table

            Output:
            ==============
            results directory: ${params.results_dir}

""".stripIndent()

process identifyMQTLs {
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: 'mQTL.parquet'

    input:
    path genotype_table
    path mynorm

    output:
    path 'mQTL.parquet'

    script:
    """
    mqtl.py $genotype_table $mynorm ${csv_gsa_manifest} ${csv_methylation_manifest} ${params.distance} ${params.alpha} ${params.slope} ${params.nCPU}
    """
}

workflow {
    identifyMQTLs(genotype_table, mynorm)
}
