#!/usr/bin/env nextflow
include { loadMynorm; loadGenotypeFrame; loadVCF } from './modules/utils.nf'
include { adjustCPUs } from './modules/utils.nf'

// Workflow Params
params.results_dir = ''
params.csv_gsa_manifest = ''
params.csv_methylation_manifest = ''
params.conversion_file = ''
params.sample_sheet = ''
params.sample_name = 'Sample_Name'
params.covars = ''

params.vep_annotations = false
params.genome_assembly = "GRCh37"

params.min_distance = 0
params.max_distance = 100000
params.alpha = 0.05
params.slope = 0.05

params.clump_p1 = 0.0001
params.clump_p2 = 0.01
params.clump_r2 = 0.50
params.clump_kb = 250

// HOMER
homer_cpus = adjustCPUs( params.CPUs )

log.info """\
            mQTL $version Worklow three [mQTL]
            ==============
            
            Config:
            ==============
            GSA manifest [--csv_gsa_manifest <path>]: ${params.csv_gsa_manifest}
            Methylation manifest [--csv_methylation_manifest <path>]: ${params.csv_methylation_manifest}
            Conversion file [--conversion_file <path>]: ${params.conversion_file}
            Genome assembly [--genome_assembly <str>]: ${params.genome_assembly} [default: GRCh37]
            Sample sheet [--sample_sheet <path>]: ${params.sample_sheet}
            Sample name column [--sample_name <str>]: ${params.sample_name} [default: Sample_Name]
            Number of CPUs [--CPUs <int>]: ${params.CPUs} [default: 10]

            mQTL config:
            ==============
            Alpha [--alpha <float>]: ${params.alpha} [default: 0.05]
            Slope [--slope <float>]: ${params.slope} [default: 0.05]
            Minimum distance [--min_distance <int>]: ${params.min_distance} [default: 0]
            Maximum distance [--max_distance <int>]: ${params.max_distance} [default: 100000]
            Model covariates [--covars list[<str>]]: ${params.covars}

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
    path conversion_file
    path mynorm
    path genotype_table
    path vcf
    path sample_sheet

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
    files = ["$csv_gsa_manifest", "$csv_methylation_manifest", "$conversion_file", "$mynorm", "$genotype_table", "$vcf", "$sample_sheet"]
    for file in files:
        assert os.path.exists(file), f"File does not exists {file}"
    
    # Assembly flag check
    assert "${params.genome_assembly}" in ["GRCh37", "GRCh38"], f"--genome_assembly should be either GRCh37 or GRCh38, is ${params.genome_assembly}"

    # VEP flag check 
    assert "${params.vep_annotations}" in ["true", "false"], f"--vep_annotations should be either true or false, is ${params.vep_annotations}"

    # Distance
    assert int(${params.min_distance}) >= 0, "Min distance should be >= 0"
    assert int(${params.max_distance}) > int(${params.min_distance}), "Max distance <= min distance between SNP and CpG"

    # Sample sheet
    ss = pd.read_csv("$sample_sheet")
    assert "${params.sample_name}" in ss, "Sample Name not in sample sheet"

    if "${params.covars}":
        covars = [covar.strip() for covar in "${params.covars}".split(",")]
        for covar in covars:
            assert covar in ss.columns, f"{covar} not available in sample sheet"

        try:
            ss[covars].astype(float)
        except ValueError:
            raise Exception("Sample sheet should contain only numerical values!")
    """
}

process identifyMQTLs {
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: 'mQTL.parquet'

    input:
    path genotype_table
    path mynorm
    path csv_gsa_manifest
    path csv_methylation_manifest
    path conversion_file
    path sample_sheet

    output:
    path 'mQTL.parquet'

    script:
    """

    mqtl.py $genotype_table $mynorm $csv_gsa_manifest \
    $csv_methylation_manifest $conversion_file $sample_sheet \
    ${params.sample_name} ${params.covars} ${params.min_distance} \
    ${params.max_distance} ${params.CPUs}
    
    """
}

process filterMQTLs {
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: 'filtered_mQTL.parquet'
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: '*pval.txt'

    input:
    path mqtls

    output:
    path 'filtered_mQTL.parquet', emit: mqtl_filtered
    path 'cpg_pval.txt', emit: cpg_pval_list
    path 'rs_pval.txt', emit: rs_pval_list
    path 'snp_list.txt', emit: snp_list

    script:
    """
    #!/usr/bin/python3
    
    import pandas as pd 

    mqtl = pd.read_parquet("$mqtls")
    mqtl[["cpg", "FDR"]].rename({"cpg": "CpG", "FDR": "P"}, axis=1).groupby("CpG").min().to_csv("cpg_pval.txt", sep="\t")
    mqtl[["RsID", "FDR"]].rename({"RsID": "SNP", "FDR": "P"}, axis=1).groupby("SNP").min().to_csv("rs_pval.txt", sep="\t")

    filtered = mqtl[(mqtl["|slope|"] >= float(${params.slope})) & (mqtl.FDR <= float(${params.alpha}))]
    filtered.to_parquet("filtered_mQTL.parquet")

    # Raw SNPs [Illumina probes] not RSids! To make this list consistent with VCF file!
    filtered.snp.drop_duplicates().to_csv("snp_list.txt", index=False, header=None, sep="\t")

    """
}

process anotateSNPs {
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: 'vep_report*'
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: 'mQTL.vcf.gz'

    input:
    path vcf
    path snp_list

    output:
    path 'vep_report*'
    path 'mQTL.vcf.gz'

    script:
    """

    bcftools view -i 'ID=@$snp_list' --threads ${params.CPUs} -O z -o mQTL.vcf.gz $vcf
    vep -i mQTL.vcf.gz -o vep_report --assembly ${params.genome_assembly} --cache --offline --fork ${params.CPUs} --dir_cache "/usr/local/.vep"

    """
}

process exportBEDfiles {
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: '*.bed'

    input:
    path mqtl

    output:
    tuple val('SNP'), path('snpInput.bed'), path('snpBg.bed'), emit: snp
    tuple val('CpG'), path('cpgInput.bed'), path('cpgBg.bed'), emit: cpg

    script:
    """
    #!/usr/bin/python3
    
    import pandas as pd 

    mqtl = pd.read_parquet("$mqtl")
    mqtl[["cpg pos", "snp pos"]] = mqtl[["cpg pos", "snp pos"]].astype(int)

    cpg = mqtl[["CHR", "cpg", "cpg pos"]]
    cpg["start"] = cpg["cpg pos"] - 1
    cpg[["CHR", "start", "cpg pos", "cpg"]].drop_duplicates().sort_values(["CHR", "start"]).to_csv(f"cpgBg.bed", index=False, header=None, sep="\t")

    snp = mqtl[["CHR", "snp", "snp pos"]]
    snp["start"] = snp["snp pos"] - 1
    snp[["CHR", "start", "snp pos", "snp"]].drop_duplicates().sort_values(["CHR", "start"]).to_csv(f"snpBg.bed", index=False, header=None, sep="\t")

    mqtl = mqtl[(mqtl["|slope|"] >= float(${params.slope})) & (mqtl.FDR <= float(${params.alpha}))]

    cpg = mqtl[["CHR", "cpg", "cpg pos"]]
    cpg["start"] = cpg["cpg pos"] - 1
    cpg[["CHR", "start", "cpg pos", "cpg"]].drop_duplicates().sort_values(["CHR", "start"]).to_csv(f"cpgInput.bed", index=False, header=None, sep="\t")

    snp = mqtl[["CHR", "snp", "snp pos"]]
    snp["start"] = snp["snp pos"] - 1
    snp[["CHR", "start", "snp pos", "snp"]].drop_duplicates().sort_values(["CHR", "start"]).to_csv(f"snpInput.bed", index=False, header=None, sep="\t")

    """
}


process runHomer {
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: '*_homer.*'

    input:
    tuple val(type), path(input), path(bg)

    output:
    path '*_homer.*'

    script:
    """
    if [[ "${params.genome_assembly}" == "GRCh37" ]]; then
        findMotifsGenome.pl $input hg19 "homer/$type" -bg $bg -mask -nomotif -h -p ${homer_cpus}
    else
        findMotifsGenome.pl $input hg38 "homer/$type" -bg $bg -mask -nomotif -h -p ${homer_cpus}
    fi

    cp homer/$type/knownResults.html $type"_homer.html"
    cp homer/$type/knownResults.txt $type"_homer.txt"
    """
}


process clumping {
    publishDir "$params.results_dir/flow_three", mode: 'copy', overwrite: true, pattern: 'plink.clumped', optional: true

    input: 
    path vcf
    path rs

    output:
    path 'plink.clumped', optional: true

    script: 
    """

    plink --vcf $vcf --clump $rs --clump-p1 ${params.clump_p1} --clump-p2 ${params.clump_p2} --clump-r2 ${params.clump_r2} --clump-kb ${params.clump_kb}

    """
}

workflow {
    // Params
    csv_gsa_manifest = file(params.csv_gsa_manifest)
    csv_methylation_manifest = file(params.csv_methylation_manifest)
    conversion_file = file(params.conversion_file)
    sample_sheet = file(params.sample_sheet)

    mynorm = loadMynorm( params.results_dir )
    genotype_table = loadGenotypeFrame( params.results_dir )
    vcf = loadVCF( params.results_dir )

    // Validate params
    validateParams( csv_gsa_manifest, csv_methylation_manifest, conversion_file, mynorm, genotype_table, vcf, sample_sheet )

    // Workflow
    mqtls = identifyMQTLs( genotype_table, mynorm, csv_gsa_manifest, csv_methylation_manifest, conversion_file, sample_sheet )
    filtered_data = filterMQTLs( mqtls )
    
    beds = exportBEDfiles( mqtls )
    beds.snp.mix(beds.cpg).set { combinedChannel }

    runHomer( combinedChannel )
    anotateSNPs( vcf, filtered_data.snp_list )
    clumping( vcf, filtered_data.rs_pval_list )
}
