#!/usr/bin/env nextflow

params.bmp_manifest = '~/projects/mQTL/GSA/resources/GSA-24v3-0_A1.bpm'
params.csv_manifest = '~/projects/mQTL/GSA/resources/GSA-24v3-0_A1.csv'
params.cluster_file = '~/projects/mQTL/GSA/resources/GSA-24v3-0_A1_ClusterFile.egt'
params.reference_fa = '~/projects/mQTL/GSA/resources/GRCh37_genome.fa'
params.CPUs = '20'

params.gsa_idats_dir = 'GSA/'
params.epic_idats_dir = 'EPIC/'
params.results_dir = 'results/'


log.info """\
            mQTL v1.0 PIPELINE
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
            EPIC idats: ${params.epic_idats_dir}
            
            Output:
            ==============
            results directory: ${params.results_dir}

""".stripIndent()

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
    path vcf_files_dir, type: 'dir'

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
    publishDir params.results_dir, mode: 'copy', overwrite: true

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
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: '*.frq'
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: '*.traw'

    input:
    path 'merged.bcf'

    output: 
    path 'genotype.*'

    script: 
    """
    plink --bcf merged.bcf --freq --maf 0.05 --hwe 1e-50 --geno 0.1 --autosome --recode A-transpose -out genotype

    """
}

workflow {
    gsa_idats_dir_ch = Channel.fromPath( params.gsa_idats_dir )
    gtc_dir_ch = callGenotypes( gsa_idats_dir_ch )
    vcf_dir_ch = GtcToVcf( gtc_dir_ch )
    
    indexed_vcf_dir_ch = indexVCF( vcf_dir_ch )
    merged_bcf_ch = mergeVCF( indexed_vcf_dir_ch )
    createGenotypeFrame( merged_bcf_ch )
}
