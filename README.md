# mQTL-workflow
This Nextflow project is designed to identify methylation quantitative trait loci (mQTL), which are interactions between the genome and the methylome. The project consists of three main workflows:

1. `flow_one` - Preprocess SNP arrays data (Illumina).
2. `flow_two` - Preprocess methylation arrays data (Illumina).
3. `flow_three` - Identify mQTLs.

### Getting Started

#### Prerequisites
Ensure you have the following software installed:

- Java >11, <=22
- Nextflow
- Docker

#### Docker Container
This projects contains two docker containers for Python3.10 as well as R4.4 environment used as execution environment for Nextflow pipelines.
To build images, use the following commands:

```
docker build -t python-image -f containers/python_env .
docker build -t r-image -f containers/R_env .
```

To validate builds:

```
docker run -it python-image
# Python 3.10.14

docker run -it r-image
# R version 4.4.1 (2024-06-14) -- "Race for Your Life"
```

#### Workflow one: SNP Array Data Preprocessing
The first workflow preprocesses SNP array data from Illumina microarrays. It involves the following steps:

```
0. Sanity check
1. Convert idat to gtc: Convert raw idat files to gtc format.
2. Convert gtc to vcf: Convert gtc files to VCF format.
3. Index and merge VCF files: Index and merge VCF files into a BCF file.
4. Filtering and formatting:
    1. Use PLINK to filter the merged BCF file.
    2. Apply filters: minor allele frequency (--maf 0.05), Hardy-Weinberg equilibrium (--hwe 1e-50), genotype missingness (--geno 0.1), and autosome restriction.
    3. Recode to a transposed format.
5. Exporting results files (merged.bcf, *.frq, *.traw as well as recoded traw file - *.parquet).
```

Workflow command:

```
nextflow run flow_one.nf --bmp_manifest <path> --csv_manifest <path> --cluster_file <path> --reference_fa <path> --CPUs <int> --gsa_idats_dir <path> --results_dir <path> --sample_sheet <path> --array_position <str> --sample_name <str> -with-docker python-image
```

- bmp_manifest: path to GSA array specific BMP manifest file
- csv_manifest: path to GSA array specific CSV manifest file
- cluster_file: path to GSA array specific cluster file
- reference_fa: path to reference genome
- CPUs: number of CPUs to use, should be psoitive integer 
- gsa_idats_dir: directory comprising GSA idats to analyse
- results_dir: results directory
- sample_sheet: sample sheet contianing at least info about `sample name` and `array position`
- array_position: sample sheet column name containing information about array position (sentrix ID and sentrix position) e.g. `205723740073_R03C02`
- sample_name: sample sheet column name containing sample name

**NOTE** Only IDATs common between `gsa_idats_dir/` and  `sample_sheet` are going to be analysed!
**NOTE** BMP manifest, CSV manifest as well as cluster file are product specific files therefore should be downloaded from Illumina Product Files page e.g. [GSA v3.0](https://emea.support.illumina.com/downloads/infinium-global-screening-array-v3-0-product-files.html)
**NOTE** Reference genome should be hg19 for EPICv1 and hg38 for EPICv2, and downloaded from [iGenomes](https://emea.support.illumina.com/sequencing/sequencing_software/igenome.html)

Workflow config example:

```
mQTL v1.0 Worklow one [SNPs arrays]
==============

Config:
==============
Reference genome: resources/GRCh37_genome.fa
BMP manifest: resources/GSA-24v3-0_A1.bpm
CSV manifest: resources/GSA-24v3-0_A1.csv
Cluster file: resources/GSA-24v3-0_A1_ClusterFile.egt
Number of CPUs: 20

Input:
==============
GSA idats: GSA/
Sample sheet: sample_sheet.csv
Array position column: ArrayPicker
Sample name column: Sample_Name

Output:
==============
results directory: results/
```

#### Workflow two: Methylation Array Data Preprocessing
Note: This workflow is not yet implemented.

This workflow will preprocess methylation array data from Illumina. Planned steps include:

```
0. Sanity check
1. Sesame preprocessing: use the Sesame R package for preprocessing the methylation data.
2. Cell fraction correction: use linear models to adjust for cellular composition in the samples.
```

#### Workflow three: mQTL Identification
Note: This workflow is not yet implemented.

This workflow will analyze data from flow_one and flow_two to identify mQTLs. Planned steps include:

```
1. mQTL analysis: Integrate genetic and methylation data to identify mQTLs.
2. Clumping: Perform clumping to group SNPs based on linkage disequilibrium.
3. Enrichment analysis: Conduct enrichment analysis to interpret the biological significance of the identified mQTLs.
```

### Contributing
Contributions to the project are welcome. Please fork the repository, make your changes, and submit a pull request.

### License
This project is licensed under the MIT License.

### Contact
For questions or issues, please contact [Jan Bińkowski] at [jan.binkowski@pum.edu.pl].
