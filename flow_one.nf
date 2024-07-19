#!/usr/bin/env nextflow

// Workflow Params
// CLI params - array specific options
params.bmp_manifest = 'resources/GSA-24v3-0_A1.bpm'
params.csv_manifest = 'resources/GSA-24v3-0_A1.csv'
params.cluster_file = 'resources/GSA-24v3-0_A1_ClusterFile.egt'
params.reference_fa = 'resources/GRCh37_genome.fa'

// CLI params - experiment specific options
params.gsa_idats_dir = 'GSA/'
params.results_dir = 'results/'
params.sample_sheet = 'sample_sheet.csv'
params.array_position = 'ArrayPicker'
params.sample_name = 'Sample_Name'

log.info """\
            mQTL ${version} Worklow one [SNPs arrays]
            ==============
            
            Config:
            ==============
            Reference genome: ${params.reference_fa}
            BMP manifest: ${params.bmp_manifest}
            CSV manifest: ${params.csv_manifest}
            Cluster file: ${params.cluster_file}
            Number of CPUs: ${params.CPUs}

            Input:
            ==============
            GSA idats: ${params.gsa_idats_dir}
            Sample sheet: ${params.sample_sheet}
            Array position column: ${params.array_position}
            Sample name column: ${params.sample_name}

            Output:
            ==============
            results directory: ${params.results_dir}

""".stripIndent()

idats = file( params.gsa_idats_dir, checkIfExists: true )
sample_sheet = file( params.sample_sheet, checkIfExists: true )
bmp_manifest = file( params.bmp_manifest, checkIfExists: true )
csv_manifest = file( params.csv_manifest, checkIfExists: true )
cluster_file = file( params.cluster_file, checkIfExists: true )
reference_fa = file( params.reference_fa, checkIfExists: true )


process prepareIDATs {
    input:
    path gsa

    output: 
    path prepared_idats_directory, type: 'dir'
    
    script:
    """
    
    prepare_idats.py $gsa ${sample_sheet} ${params.array_position} prepared_idats_directory

    """
}

process callGenotypes {
    input: 
    path idats_dir

    output:
    path gtc_files_dir, type: 'dir'

    script:
    """
    array-analysis-cli genotype call --bpm-manifest ${bmp_manifest} --cluster-file ${cluster_file} --idat-folder $idats_dir --num-threads ${params.CPUs} --output-folder gtc_files_dir

    """
}


process GtcToVcf {
    input:
    path gtc_dir

    output:
    path vcf_files_dir, type: 'dir'

    script:
    """

    array-analysis-cli genotype gtc-to-vcf --csv-manifest ${csv_manifest} --bpm-manifest ${bmp_manifest} --gtc-folder $gtc_dir --genome-fasta-file ${reference_fa} --output-folder vcf_files_dir

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

process createGenotypeFrame {
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: 'genotype_table.traw'
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: '*.frq'

    input:
    path 'merged.bcf'

    output: 
    path 'genotype_table.traw', emit: frame
    path '*.frq'

    script: 
    """

    plink --bcf merged.bcf --freq --maf 0.05 --hwe 1e-50 --geno 0.1 --autosome --recode A-transpose -out genotype_table

    """
}

process refineGenotypeFrame {
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: 'genotype_table.parquet'

    input:
    path genotype_table

    output:
    path 'genotype_table.parquet'

    script:
    
    """
   
    refine_genotype_frame.py $genotype_table ${sample_sheet} ${params.array_position} ${params.sample_name}
    
    """   
}

workflow {

    gsa_idats_dir = prepareIDATs( idats )
    gtc_dir = callGenotypes( gsa_idats_dir )
    vcf_dir = GtcToVcf( gtc_dir )
    
    indexed_vcf_dir = indexVCF( vcf_dir )
    merged_bcf = mergeVCF( indexed_vcf_dir )
    traw = createGenotypeFrame( merged_bcf )

    refineGenotypeFrame( traw.frame )
}
