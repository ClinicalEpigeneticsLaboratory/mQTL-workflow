#!/usr/bin/env nextflow

// Workflow Params
params.bmp_manifest = ''
params.csv_manifest = ''
params.cluster_file = ''
params.reference_fa = ''

params.sample_sheet = ''
params.gsa_idats_dir = ''
params.results_dir = ''

params.array_position = 'ArrayPicker'
params.sample_name = 'Sample_Name'

log.info """\
            mQTL ${version} Worklow one [SNPs arrays]
            ==============
            
            Config:
            ==============
            Ref. genome directory [should include *.fa and *.fai files] [--reference_fa <path>]: ${params.reference_fa}
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
    input:
    path reference_fa
    path reference_fa_index
    path bmp_manifest
    path csv_manifest
    path cluster_file
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
    files = ["$reference_fa", "$reference_fa_index", "$bmp_manifest", "$csv_manifest", "$cluster_file", "$sample_sheet"]
    for file in files:
        assert os.path.exists(file), f"File does not exists {file} {os.getcwd()}"    
    
    # --array_position and --sample_name flags check
    sample_sheet_fields = pd.read_csv("$sample_sheet").columns
    assert "${params.array_position}" in sample_sheet_fields, f"${params.array_position} not present in sample sheet columns: {sample_sheet_fields}"
    assert "${params.sample_name}" in sample_sheet_fields, f"${params.sample_name} not present in sample sheet columns: {sample_sheet_fields}"

    """

}

process prepareIDATs {
    input:
    path gsa
    path sample_sheet

    output: 
    path 'prepared_idats_directory', type: 'dir'
    
    script:
    """
    
    prepare_idats.py $gsa $sample_sheet ${params.array_position} prepared_idats_directory

    """
}

process callGenotypes {
    input: 
    path bmp_manifest
    path cluster_file
    path idats_dir    

    output:
    path 'gtc_files_dir', type: 'dir'

    script:
    """
    array-analysis-cli genotype call --bpm-manifest $bmp_manifest --cluster-file $cluster_file --idat-folder $idats_dir --num-threads ${params.CPUs} --output-folder gtc_files_dir

    """
}

process GtcToVcf {
    input:
    path csv_manifest
    path bmp_manifest
    path reference_fa
    path reference_fa_index
    path gtc_dir


    output:
    path 'vcf_files_dir', type: 'dir'

    script:
    """

    array-analysis-cli genotype gtc-to-vcf --csv-manifest $csv_manifest --bpm-manifest $bmp_manifest --genome-fasta-file $reference_fa --gtc-folder $gtc_dir --output-folder vcf_files_dir

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

    input:
    path filtered_merged_bcf_file

    output: 
    path 'genotype_table.traw', emit: frame

    script: 
    """

    plink --bcf $filtered_merged_bcf_file --maf ${params.MAF} --hwe ${params.HWE} --geno ${params.GENO} --autosome --recode A-transpose -out genotype_table

    """
}

process refineGenotypeFrame {
    publishDir "$params.results_dir/flow_one", mode: 'copy', overwrite: true, pattern: 'genotype_table.parquet'

    input:
    path genotype_table
    path sample_sheet

    output:
    path 'genotype_table.parquet'

    script:
    
    """
   
    refine_genotype_frame.py $genotype_table $sample_sheet ${params.array_position} ${params.sample_name}
    
    """   
}

workflow {
    // Params
    reference_fa = file("$params.reference_fa/*.fa")
    reference_fa_index = file("$params.reference_fa/*.fai")
    
    bmp_manifest = file(params.bmp_manifest)
    csv_manifest = file(params.csv_manifest)
    cluster_file = file(params.cluster_file)
    sample_sheet = file(params.sample_sheet)

    // Params validation
    validateParams( reference_fa, reference_fa_index, bmp_manifest, csv_manifest, cluster_file, sample_sheet )

    // Workflow
    idats_dir = file(params.gsa_idats_dir)
    prepared_idats = prepareIDATs( idats_dir,  sample_sheet )

    gtc_dir = callGenotypes( bmp_manifest, cluster_file, prepared_idats )
    vcf_dir = GtcToVcf( csv_manifest, bmp_manifest, reference_fa, reference_fa_index, gtc_dir )
    
    indexed_vcf_dir = indexVCF( vcf_dir )
    merged_bcf = mergeVCF( indexed_vcf_dir )
    filtered_bcf = filterBCF( merged_bcf )
    genotypes = createGenotypeFrame( filtered_bcf )

    refineGenotypeFrame( genotypes.frame, sample_sheet )
}
