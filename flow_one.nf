#!/usr/bin/env nextflow

// Workflow Params
params.bmp_manifest = file('resources/GSA-24v3-0_A1.bpm')
params.csv_manifest = file('resources/GSA-24v3-0_A1.csv')
params.cluster_file = file('resources/GSA-24v3-0_A1_ClusterFile.egt')
params.reference_fa = file('resources/GRCh37_genome.fa')

params.sample_sheet = file('sample_sheet.csv')
params.gsa_idats_dir = file('GSA/')
params.results_dir = file('results/')

params.array_position = 'ArrayPicker'
params.sample_name = 'Sample_Name'

log.info """\
            mQTL ${version} Worklow one [SNPs arrays]
            ==============
            
            Config:
            ==============
            Reference genome [--reference_fa <path>]: ${params.reference_fa}
            BMP manifest [--bmp_manifest <path>]: ${params.bmp_manifest}
            CSV manifest [--csv_manifest <path>]: ${params.csv_manifest}
            Cluster file [--cluster_file <path>]: ${params.cluster_file}
            Number of CPUs [--CPUs <int>]: ${params.CPUs} [default: 10]

            PLINK params:
            ==============
            MAF [--MAF <float>]: ${params.MAF} [default: 0.01]
            HWE [--HWE <float>]: ${params.HWE} [default: 1e-50]
            GENO [--GENO <float>]: ${params.GENO} [default: 0.1]

            Input:
            ==============
            GSA idats [--gsa_idats_dir <path>]: ${params.gsa_idats_dir}
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
    files = ["${params.reference_fa}", "${params.csv_manifest}", "${params.cluster_file}", "${params.reference_fa}", "${params.sample_sheet}"]
    for file in files:
        assert os.path.exists(file), f"File does not exists {file}"

    # --array_position and --sample_name flags check
    sample_sheet_fields = pd.read_csv("${params.sample_sheet}").columns
    assert "${params.array_position}" in sample_sheet_fields, f"${params.array_position} not present in sample sheet columns: {sample_sheet_fields}"
    assert "${params.sample_name}" in sample_sheet_fields, f"${params.sample_name} not present in sample sheet columns: {sample_sheet_fields}"

    """

}

process prepareIDATs {
    input:
    path gsa

    output: 
    path 'prepared_idats_directory', type: 'dir'
    
    script:
    """
    
    prepare_idats.py $gsa ${params.sample_sheet} ${params.array_position} prepared_idats_directory

    """
}

process callGenotypes {
    input: 
    path idats_dir

    output:
    path 'gtc_files_dir', type: 'dir'

    script:
    """
    array-analysis-cli genotype call --bpm-manifest ${params.bmp_manifest} --cluster-file ${params.cluster_file} --idat-folder $idats_dir --num-threads ${params.CPUs} --output-folder gtc_files_dir

    """
}

process GtcToVcf {
    input:
    path gtc_dir

    output:
    path 'vcf_files_dir', type: 'dir'

    script:
    """

    array-analysis-cli genotype gtc-to-vcf --csv-manifest ${params.csv_manifest} --bpm-manifest ${params.bmp_manifest} --gtc-folder $gtc_dir --genome-fasta-file ${params.reference_fa} --output-folder vcf_files_dir

    """
}

process indexVCF {
    input:
    path vcf_files_dir

    output:
    path vcf_files_dir, type: 'dir'

    script:
    """
    cd $vcf_files_dir

    ls *.vcf | xargs -n1 bgzip -@ ${params.CPUs}
    ls *.vcf.gz | xargs -n1 tabix -f -p vcf

    """
}

process mergeVCF {
    input:
    path vcf_files_dir

    output:
    path 'merged.bcf'

    script: 
    """
    
    bcftools merge --threads ${params.CPUs} -o merged.bcf -O b $vcf_files_dir/*.vcf.gz

    """
}

process filterBCF {
    publishDir "$params.results_dir/flow_one", mode: 'copy', overwrite: true, pattern: 'filtered_merged.bcf'

    input:
    path bcf_file

    output:
    path 'filtered_merged.bcf'

    script: 
    """
    
    bcftools view -m 2 -M 2 $bcf_file -O b -o filtered_merged.bcf

    """
}

process createGenotypeFrame {
    publishDir "$params.results_dir/flow_one", mode: 'copy', overwrite: true, pattern: 'genotype_table.traw'
    publishDir "$params.results_dir/flow_one", mode: 'copy', overwrite: true, pattern: '*.frq'

    input:
    path filtered_merged_bcf_file

    output: 
    path 'genotype_table.traw', emit: frame
    path '*.frq'

    script: 
    """

    plink --bcf $filtered_merged_bcf_file --freq --maf ${params.MAF} --hwe ${params.HWE} --geno ${params.GENO} --autosome --recode A-transpose -out genotype_table

    """
}

process refineGenotypeFrame {
    publishDir "$params.results_dir/flow_one", mode: 'copy', overwrite: true, pattern: 'genotype_table.parquet'

    input:
    path genotype_table

    output:
    path 'genotype_table.parquet'

    script:
    
    """
   
    refine_genotype_frame.py $genotype_table ${params.sample_sheet} ${params.array_position} ${params.sample_name}
    
    """   
}

workflow {
    validateParams()

    prepared_idats = prepareIDATs( params.gsa_idats_dir )
    gtc_dir = callGenotypes( prepared_idats )
    vcf_dir = GtcToVcf( gtc_dir )
    
    indexed_vcf_dir = indexVCF( vcf_dir )
    merged_bcf = mergeVCF( indexed_vcf_dir )
    filtered_bcf = filterBCF( merged_bcf )
    genotypes = createGenotypeFrame( filtered_bcf )

    refineGenotypeFrame( genotypes.frame )
}
