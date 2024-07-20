#!/usr/bin/env nextflow
include { loadMynorm; loadGenotypeFrame; loadBCF } from './modules/utils.nf'

// Workflow Params
params.results_dir = file('results/')
params.csv_gsa_manifest = file('resources/GSA-24v3-0_A1.csv')
params.csv_methylation_manifest = file('resources/infinium-methylationepic-v-1-0-b5-manifest-file.csv')

params.vep_annotations = false
params.genome_assembly = "GRCh37"
params.distance = 50000
params.alpha = 0.05
params.slope = 0.05

params.clump_p1 = 0.0001       
params.clump_p2 = 0.01      
params.clump_r2 = 0.50      
params.clump_kb = 250     

mynorm = loadMynorm( params.results_dir )
genotype_table = loadGenotypeFrame( params.results_dir )
bcf = loadBCF( params.results_dir )

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
            Distance [--distance <int>]: ${params.distance} [default: 50000]
            Number of CPUs [--CPUs <int>]: ${params.CPUs} [default: 10]

            PLINK config:
            ==============
            P1: [--clump_p1 <float>]: ${params.clump_p1} [default: 0.0001]
            P2: [--clump_p2 <float>]: ${params.clump_p2} [default: 0.01]
            R1: [--clump_r2 <float>]: ${params.clump_r2} [default: 0.5]
            KB: [--clump_kb <int>]: ${params.clump_kb} [default: 250]

            Input:
            ==============
            Methylation frame: $mynorm
            Genotype table: $genotype_table
            BCF file: $bcf
            
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
    files = ["${params.csv_gsa_manifest}", "${params.csv_methylation_manifest}", "$mynorm", "$genotype_table", "$bcf"]
    for file in files:
        assert os.path.exists(file), f"File does not exists {file}"

    # Assembly flag check
    assert "${params.genome_assembly}" in ["GRCh37", "GRCh38"], f"--genome_assembly should be either GRCh37 or GRCh38, is ${params.genome_assembly}"

    """
}

process identifyMQTLs {
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: 'mQTL.parquet'

    input:
    path genotype_table
    path mynorm

    output:
    path 'mQTL.parquet'

    script:
    """

    mqtl.py $genotype_table $mynorm ${params.csv_gsa_manifest} ${params.csv_methylation_manifest} ${params.distance} ${params.CPUs}
    
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
    path 'rs_list.txt', emit: rs_list
    path 'rs_pval.txt', emit: rs_pval_list

    script:
    """

    #!/usr/bin/python3
    import pandas as pd 

    mqtl = pd.read_parquet("$mqtls")
    filtered = mqtl[(mqtl.slope.abs() >= float(${params.slope})) & (mqtl.FDR <= float(${params.alpha}))]

    filtered.to_parquet("filtered_mQTL.parquet")
    filtered.cpg.drop_duplicates().to_csv("cpg_list.txt", index=False, header=None, sep="\t")

    filtered.rs.drop_duplicates().to_csv("rs_list.txt", index=False, header=None, sep="\t")
    filtered[["rs", "FDR"]].rename({"rs": "SNP", "FDR": "P"}, axis=1).drop_duplicates().to_csv("rs_pval.txt", index=False, sep="\t")

    """
}

process anotateSNPs {
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: 'report*'

    input:
    path rs_list 

    output:
    path 'report*'

    script:
    """

    vep -i $rs_list --format id -o report --assembly ${params.genome_assembly} --database

    """
}

process clumping {
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: 'plink.clumped', optional: true

    input: 
    path rs

    output:
    path 'plink.clumped', optional: true

    script: 
    """

    plink --bcf $bcf --clump $rs --clump-p1 ${params.clump_p1} --clump-p2 ${params.clump_p2} --clump-r2 ${params.clump_r2} --clump-kb ${params.clump_kb}  

    """
}

workflow {
    validateParams()

    mqtls = identifyMQTLs(genotype_table, mynorm)
    filtered_data = filterMQTLs( mqtls )

    if ( params.vep_annotations ){
        anotateSNPs( filtered_data.rs_list )
    }

    clumping( filtered_data.rs_pval_list )
}
