#!/usr/bin/env nextflow
include { loadMynorm; loadGenotypeFrame; loadBCF } from './modules/utils.nf'

// Workflow Params
params.results_dir = ''
params.csv_gsa_manifest = ''
params.csv_methylation_manifest = ''

params.vep_annotations = false
params.genome_assembly = "GRCh37"
params.distance = 250000
params.alpha = 0.05
params.slope = 0.05

params.clump_p1 = 0.0001       
params.clump_p2 = 0.01      
params.clump_r2 = 0.50      
params.clump_kb = 250     

log.info """\
            mQTL $version Worklow three [mQTL]
            ==============
            
            Config:
            ==============
            GSA manifest [--csv_gsa_manifest <path>]: ${params.csv_gsa_manifest}
            Methylation manifest [--csv_methylation_manifest <path>]: ${params.csv_methylation_manifest}
            Genome assembly [--genome_assembly <str>]: ${params.genome_assembly} [default: GRCh37 (important if VEP is ON)]
            VEP annotations [--vep_annotations <boolean: true/false>]: ${params.vep_annotations} [default: false (very slow)]
            Alpha [--alpha <float>]: ${params.alpha} [default: 0.05]
            Slope [--slope <float>]: ${params.slope} [default: 0.05]
            Distance [--distance <int>]: ${params.distance} [default: 250000]
            Number of CPUs [--CPUs <int>]: ${params.CPUs} [default: 10]

            PLINK config:
            ==============
            P1: [--clump_p1 <float>]: ${params.clump_p1} [default: 0.0001]
            P2: [--clump_p2 <float>]: ${params.clump_p2} [default: 0.01]
            R1: [--clump_r2 <float>]: ${params.clump_r2} [default: 0.5]
            KB: [--clump_kb <int>]: ${params.clump_kb} [default: 250]

            Input:
            ==============
            Results from flow_one and flow_two are expected to be in [--results_dir <path>]: ${params.results_dir}
            
            Output:
            ==============
            results directory [--results_dir <path>]: ${params.results_dir}

""".stripIndent()

process validateParams {
    input: 
    path csv_gsa_manifest
    path csv_methylation_manifest
    path mynorm
    path genotype_table
    path bcf

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
    files = ["$csv_gsa_manifest", "$csv_methylation_manifest", "$mynorm", "$genotype_table", "$bcf"]
    for file in files:
        assert os.path.exists(file), f"File does not exists {file}"
    
    # Assembly flag check
    assert "${params.genome_assembly}" in ["GRCh37", "GRCh38"], f"--genome_assembly should be either GRCh37 or GRCh38, is ${params.genome_assembly}"

    # VEP flag check 
    assert "${params.vep_annotations}" in ["true", "false"], f"--vep_annotations should be either true or false, is ${params.vep_annotations}"

    """
}

process identifyMQTLs {
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: 'mQTL.parquet'

    input:
    path genotype_table
    path mynorm
    path csv_gsa_manifest
    path csv_methylation_manifest

    output:
    path 'mQTL.parquet'

    script:
    """

    mqtl.py $genotype_table $mynorm $csv_gsa_manifest $csv_methylation_manifest ${params.distance} ${params.CPUs}
    
    """
}

process filterMQTLs {
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: 'filtered_mQTL.parquet'
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: '*.txt'

    input:
    path mqtls

    output:
    path 'filtered_mQTL.parquet', emit: mqtl_filtered
    path 'cpg_list.txt', emit: cpg_list
    path 'cpg_pval.txt', emit: cpg_pval_list
    path 'rs_list.txt', emit: rs_list
    path 'rs_pval.txt', emit: rs_pval_list

    script:
    """
    #!/usr/bin/python3
    
    import pandas as pd 

    mqtl = pd.read_parquet("$mqtls")
    mqtl[["cpg", "FDR"]].rename({"cpg": "CpG", "FDR": "P"}, axis=1).groupby("CpG").min().to_csv("cpg_pval.txt", sep="\t")
    mqtl[["rs", "FDR"]].rename({"rs": "SNP", "FDR": "P"}, axis=1).groupby("SNP").min().to_csv("rs_pval.txt", sep="\t")

    filtered = mqtl[(mqtl.slope >= float(${params.slope})) & (mqtl.FDR <= float(${params.alpha}))]
    filtered.to_parquet("filtered_mQTL.parquet")

    filtered.cpg.drop_duplicates().to_csv("cpg_list.txt", index=False, header=None, sep="\t")
    filtered.rs.drop_duplicates().to_csv("rs_list.txt", index=False, header=None, sep="\t")

    """
}

process anotateSNPs {
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: 'report*'

    input:
    path rs_list 

    output:
    path 'vep_report*'

    script:
    """

    vep -i $rs_list --format id -o vep_report --assembly ${params.genome_assembly} --database --biotype --variant_class --sift b --polyphen b

    """
}

process clumping {
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: 'plink.clumped', optional: true

    input: 
    path bcf
    path rs

    output:
    path 'plink.clumped', optional: true

    script: 
    """

    plink --bcf $bcf --clump $rs --clump-p1 ${params.clump_p1} --clump-p2 ${params.clump_p2} --clump-r2 ${params.clump_r2} --clump-kb ${params.clump_kb}  

    """
}

workflow {
    // Params
    csv_gsa_manifest = file(params.csv_gsa_manifest)
    csv_methylation_manifest = file(params.csv_methylation_manifest)

    mynorm = loadMynorm( params.results_dir )
    genotype_table = loadGenotypeFrame( params.results_dir )
    bcf = loadBCF( params.results_dir )

    // Validate params
    validateParams( csv_gsa_manifest, csv_methylation_manifest, mynorm, genotype_table, bcf )

    // Workflow
    mqtls = identifyMQTLs( genotype_table, mynorm, csv_gsa_manifest, csv_methylation_manifest )
    filtered_data = filterMQTLs( mqtls )

    if ( params.vep_annotations ){
        anotateSNPs( filtered_data.rs_list )
    }

    clumping( bcf, filtered_data.rs_pval_list )
}
